__all__ = ["IsrTaskLSST", "IsrTaskLSSTConfig"]

import numpy
import math

from . import isrFunctions
from . import isrQa
from . import linearize
from .defects import Defects

from contextlib import contextmanager
from lsst.afw.cameraGeom import NullLinearityType
import lsst.pex.config as pexConfig
import lsst.afw.math as afwMath
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as cT
from lsst.daf.butler import DimensionGraph
from lsst.meas.algorithms.detection import SourceDetectionTask

from .overscan import OverscanCorrectionTask
from .assembleCcdTask import AssembleCcdTask
from .deferredCharge import DeferredChargeTask
from .crosstalk import CrosstalkTask
from .masking import MaskingTask
from .isrStatistics import IsrStatisticsTask
from .isr import maskNans


def crosstalkSourceLookup(datasetType, registry, quantumDataId, collections):
    """Lookup function to identify crosstalkSource entries.

    This should return an empty list under most circumstances.  Only
    when inter-chip crosstalk has been identified should this be
    populated.

    Parameters
    ----------
    datasetType : `str`
        Dataset to lookup.
    registry : `lsst.daf.butler.Registry`
        Butler registry to query.
    quantumDataId : `lsst.daf.butler.ExpandedDataCoordinate`
        Data id to transform to identify crosstalkSources.  The
        ``detector`` entry will be stripped.
    collections : `lsst.daf.butler.CollectionSearch`
        Collections to search through.

    Returns
    -------
    results : `list` [`lsst.daf.butler.DatasetRef`]
        List of datasets that match the query that will be used as
        crosstalkSources.
    """
    newDataId = quantumDataId.subset(DimensionGraph(registry.dimensions, names=["instrument", "exposure"]))
    results = set(registry.queryDatasets(datasetType, collections=collections, dataId=newDataId,
                                         findFirst=True))
    # In some contexts, calling `.expanded()` to expand all data IDs in the
    # query results can be a lot faster because it vectorizes lookups.  But in
    # this case, expandDataId shouldn't need to hit the database at all in the
    # steady state, because only the detector record is unknown and those are
    # cached in the registry.
    return [ref.expanded(registry.expandDataId(ref.dataId, records=newDataId.records)) for ref in results]


class IsrTaskLSSTConnections(pipeBase.PipelineTaskConnections,
                             dimensions={"instrument", "exposure", "detector"},
                             defaultTemplates={}):
    ccdExposure = cT.Input(
        name="raw",
        doc="Input exposure to process.",
        storageClass="Exposure",
        dimensions=["instrument", "exposure", "detector"],
    )
    camera = cT.PrerequisiteInput(
        name="camera",
        storageClass="Camera",
        doc="Input camera to construct complete exposures.",
        dimensions=["instrument"],
        isCalibration=True,
    )
    dnlLUT = cT.PrerequisiteInput(
        name="dnlLUT",
        doc="Look-up table for differential non-linearity.",
        storageClass="IsrCalib",
        dimensions=["instrument", "exposure", "detector"],
        isCalibration=True,
        # TODO DM 36636
    )
    bias = cT.PrerequisiteInput(
        name="bias",
        doc="Input bias calibration.",
        storageClass="ExposureF",
        dimensions=["instrument", "detector"],
        isCalibration=True,
    )
    deferredChargeCalib = cT.PrerequisiteInput(
        name="cpCtiCalib",
        doc="Deferred charge/CTI correction dataset.",
        storageClass="IsrCalib",
        dimensions=["instrument", "detector"],
        isCalibration=True,
    )
    linearizer = cT.PrerequisiteInput(
        name='linearizer',
        storageClass="Linearizer",
        doc="Linearity correction calibration.",
        dimensions=["instrument", "detector"],
        isCalibration=True,
        minimum=0,  # can fall back to cameraGeom
    )
    ptc = cT.PrerequisiteInput(
        name="ptc",
        doc="Input Photon Transfer Curve dataset",
        storageClass="PhotonTransferCurveDataset",
        dimensions=["instrument", "detector"],
        isCalibration=True,
    )
    crosstalk = cT.PrerequisiteInput(
        name="crosstalk",
        doc="Input crosstalk object",
        storageClass="CrosstalkCalib",
        dimensions=["instrument", "detector"],
        isCalibration=True,
        minimum=0,  # can fall back to cameraGeom
    )
    crosstalkSources = cT.PrerequisiteInput(
        name="isrOverscanCorrected",
        doc="Overscan corrected input images.",
        storageClass="Exposure",
        dimensions=["instrument", "exposure", "detector"],
        deferLoad=True,
        multiple=True,
        lookupFunction=crosstalkSourceLookup,
        minimum=0,  # not needed for all instruments, no config to control this
    )
    defects = cT.PrerequisiteInput(
        name='defects',
        doc="Input defect tables.",
        storageClass="Defects",
        dimensions=["instrument", "detector"],
        isCalibration=True,
    )
    bfKernel = cT.PrerequisiteInput(
        name='brighterFatterKernel',
        doc="Newer complete kernel + gain solutions.",
        storageClass="BrighterFatterKernel",
        dimensions=["instrument", "detector"],
        isCalibration=True,
        minimum=0,  # can use either bfKernel or newBFKernel
    )
    dark = cT.PrerequisiteInput(
        name='dark',
        doc="Input dark calibration.",
        storageClass="ExposureF",
        dimensions=["instrument", "detector"],
        isCalibration=True,
    )
    outputExposure = cT.Output(
        name='postISRCCD',
        doc="Output ISR processed exposure.",
        storageClass="Exposure",
        dimensions=["instrument", "exposure", "detector"],
    )
    preInterpExposure = cT.Output(
        name='preInterpISRCCD',
        doc="Output ISR processed exposure, with pixels left uninterpolated.",
        storageClass="ExposureF",
        dimensions=["instrument", "exposure", "detector"],
    )
    outputStatistics = cT.Output(
        name="isrStatistics",
        doc="Output of additional statistics table.",
        storageClass="StructuredDataDict",
        dimensions=["instrument", "exposure", "detector"],
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if config.doDiffNonLinearCorrection is not True:
            self.prerequisiteInputs.remove("dnlLUT")
        if config.doBias is not True:
            self.prerequisiteInputs.remove("bias")
        if config.doDeferredCharge is not True:
            self.prerequisiteInputs.remove("deferredChargeCalib")
        if config.doLinearize is not True:
            self.prerequisiteInputs.remove("linearizer")
        if config.usePtcGains is not True and config.usePtcReadNoise is not True:
            self.prerequisiteInputs.remove("ptc")
        if config.doCrosstalk is not True:
            self.prerequisiteInputs.remove("crosstalkSources")
            self.prerequisiteInputs.remove("crosstalk")
        if config.doDefect is not True:
            self.prerequisiteInputs.remove("defects")
        if config.doBrighterFatter is not True:
            self.prerequisiteInputs.remove("bfKernel")
        if config.doDark is not True:
            self.prerequisiteInputs.remove("dark")

        if config.doWrite is not True:
            self.outputs.remove("outputExposure")
            self.outputs.remove("preInterpExposure")

        if config.doSaveInterpPixels is not True:
            self.outputs.remove("preInterpExposure")

        if config.doCalculateStatistics is not True:
            self.outputs.remove("outputStatistics")


class IsrTaskLSSTConfig(pipeBase.PipelineTaskConfig,
                        pipelineConnections=IsrTaskLSSTConnections):
    """Configuration parameters for IsrTaskLSST.

    Items are grouped in the order in which they are executed by the task.
    """
    datasetType = pexConfig.Field(
        dtype=str,
        doc="Dataset type for input data; users will typically leave this alone, "
        "but camera-specific ISR tasks will override it.",
        default="raw",
    )
    expectWcs = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Expect input science images to have a WCS (set False for e.g. spectrographs)."
    )
    qa = pexConfig.ConfigField(
        dtype=isrQa.IsrQaConfig,
        doc="QA related configuration options.",
    )
    doHeaderProvenance = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Write calibration identifiers into output exposure header.",
    )

    # Differential non-linearity correction.
    doDiffNonLinearCorrection = pexConfig.Field(
        dtype=bool,
        doc="Do differential non-linearity correction?",
        default=True,
    )

    doOverscan = pexConfig.Field(
        dtype=bool,
        doc="Do overscan subtraction?",
        default=True,
    )
    overscan = pexConfig.ConfigurableField(
        target=OverscanCorrectionTask,
        doc="Overscan subtraction task for image segments.",
    )

    # Amplifier to CCD assembly configuration.
    doAssembleCcd = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Assemble amp-level exposures into a ccd-level exposure?"
    )
    assembleCcd = pexConfig.ConfigurableField(
        target=AssembleCcdTask,
        doc="CCD assembly task.",
    )

    # Combine snaps.
    doSnapCombine = pexConfig.Field(
        dtype=bool,
        doc="Combine snaps?",
        default=False,
    )

    # Bias subtraction.
    doBias = pexConfig.Field(
        dtype=bool,
        doc="Apply bias frame correction?",
        default=True,
    )

    # Deferred charge correction.
    doDeferredCharge = pexConfig.Field(
        dtype=bool,
        doc="Apply deferred charge correction?",
        default=False,
    )
    deferredChargeCorrection = pexConfig.ConfigurableField(
        target=DeferredChargeTask,
        doc="Deferred charge correction task.",
    )

    # Linearization.
    doLinearize = pexConfig.Field(
        dtype=bool,
        doc="Correct for nonlinearity of the detector's response?",
        default=True,
    )

    # Normalize gain.
    doGainNormalize = pexConfig.Field(
        dtype=bool,
        doc="Normalize by the gain.",
        default=True,
    )

    # Variance construction.
    doVariance = pexConfig.Field(
        dtype=bool,
        doc="Calculate variance?",
        default=True
    )
    gain = pexConfig.Field(
        dtype=float,
        doc="The gain to use if no Detector is present in the Exposure (ignored if NaN).",
        default=float("NaN"),
    )
    usePtcGains = pexConfig.Field(
        dtype=bool,
        doc="Use the gain values from the Photon Transfer Curve?",
        default=False,
    )
    readNoise = pexConfig.Field(
        dtype=float,
        doc="The read noise to use if no Detector is present in the Exposure.",
        default=0.0,
    )
    doEmpiricalReadNoise = pexConfig.Field(
        dtype=bool,
        doc="Calculate empirical read noise instead of value from AmpInfo data?",
        default=False,
    )
    usePtcReadNoise = pexConfig.Field(
        dtype=bool,
        doc="Use readnoise values from the Photon Transfer Curve?",
        default=False,
    )
    maskNegativeVariance = pexConfig.Field(
        dtype=bool,
        doc="Mask pixels that claim a negative variance.  This likely indicates a failure "
        "in the measurement of the overscan at an edge due to the data falling off faster "
        "than the overscan model can account for it.",
        default=True,
    )
    negativeVarianceMaskName = pexConfig.Field(
        dtype=str,
        doc="Mask plane to use to mark pixels with negative variance, if `maskNegativeVariance` is True.",
        default="BAD",
    )
    saturatedMaskName = pexConfig.Field(
        dtype=str,
        doc="Name of mask plane to use in saturation detection and interpolation.",
        default="SAT",
    )
    suspectMaskName = pexConfig.Field(
        dtype=str,
        doc="Name of mask plane to use for suspect pixels.",
        default="SUSPECT",
    )

    # Crosstalk.
    doCrosstalk = pexConfig.Field(
        dtype=bool,
        doc="Apply intra-CCD crosstalk correction?",
        default=False,
    )
    doCrosstalkBeforeAssemble = pexConfig.Field(
        dtype=bool,
        doc="Apply crosstalk correction before CCD assembly?",
        default=False,
    )
    crosstalk = pexConfig.ConfigurableField(
        target=CrosstalkTask,
        doc="Intra-CCD crosstalk correction.",
    )

    # Masking options.
    doDefect = pexConfig.Field(
        dtype=bool,
        doc="Apply correction for CCD defects, e.g. hot pixels?",
        default=True,
    )
    doNanMasking = pexConfig.Field(
        dtype=bool,
        doc="Mask non-finite (NAN, inf) pixels.",
        default=True,
    )
    doWidenSaturationTrails = pexConfig.Field(
        dtype=bool,
        doc="Widen bleed trails based on their width.",
        default=True,
    )
    doCameraSpecificMasking = pexConfig.Field(
        dtype=bool,
        doc="Mask camera-specific bad regions?",
        default=False,
    )
    masking = pexConfig.ConfigurableField(
        target=MaskingTask,
        doc="Masking task."
    )

    # Interpolation options.
    doInterpolate = pexConfig.Field(
        dtype=bool,
        doc="Interpolate masked pixels?",
        default=True,
    )
    maskListToInterpolate = pexConfig.ListField(
        dtype=str,
        doc="List of mask planes that should be interpolated.",
        default=['SAT', 'BAD'],
    )
    doSaveInterpPixels = pexConfig.Field(
        dtype=bool,
        doc="Save a copy of the pre-interpolated pixel values?",
        default=False,
    )

    # Initial masking options.
    doSetBadRegions = pexConfig.Field(
        dtype=bool,
        doc="Should we set the level of all BAD patches of the chip to the chip's average value?",
        default=True,
    )

    # Brighter-Fatter correction.
    doBrighterFatter = pexConfig.Field(
        dtype=bool,
        doc="Apply the brighter-fatter correction?",
        default=False,
    )
    brighterFatterLevel = pexConfig.ChoiceField(
        dtype=str,
        doc="The level at which to correct for brighter-fatter.",
        allowed={
            "AMP": "Every amplifier treated separately.",
            "DETECTOR": "One kernel per detector.",
        },
        default="DETECTOR",
    )
    brighterFatterMaxIter = pexConfig.Field(
        dtype=int,
        doc="Maximum number of iterations for the brighter-fatter correction.",
        default=10,
    )
    brighterFatterThreshold = pexConfig.Field(
        dtype=float,
        doc="Threshold used to stop iterating the brighter-fatter correction.  It is the "
        "absolute value of the difference between the current corrected image and the one "
        "from the previous iteration summed over all the pixels.",
        default=1000,
    )
    brighterFatterApplyGain = pexConfig.Field(
        dtype=bool,
        doc="Should the gain be applied when applying the brighter-fatter correction?",
        default=True,
    )
    brighterFatterMaskListToInterpolate = pexConfig.ListField(
        dtype=str,
        doc="List of mask planes that should be interpolated over when applying the brighter-fatter "
        "correction.",
        default=["SAT", "BAD", "NO_DATA", "UNMASKEDNAN"],
    )
    brighterFatterMaskGrowSize = pexConfig.Field(
        dtype=int,
        doc="Number of pixels to grow the masks listed in config.brighterFatterMaskListToInterpolate "
        "when brighter-fatter correction is applied.",
        default=0,
    )
    fwhm = pexConfig.Field(
        dtype=float,
        doc="FWHM of PSF in arcseconds.",
        default=1.0,
    )
    growSaturationFootprintSize = pexConfig.Field(
        dtype=int,
        doc="Number of pixels by which to grow the saturation footprints.",
        default=1,
    )
    brighterFatterMaskListToInterpolate = pexConfig.ListField(
        dtype=str,
        doc="List of mask planes that should be interpolated over when applying the brighter-fatter."
        "correction.",
        default=["SAT", "BAD", "NO_DATA", "UNMASKEDNAN"],
    )

    # Dark subtraction.
    doDark = pexConfig.Field(
        dtype=bool,
        doc="Apply dark frame correction.",
        default=True,
    )

    # Post-ISR operations to make fluence image.
    doFluence = pexConfig.Field(
        dtype=bool,
        doc="Apply post-ISR operations?",
        default=True,
    )

    # Write the outputs to disk. If ISR is run as a subtask, this may not
    # be needed.
    doWrite = pexConfig.Field(
        dtype=bool,
        doc="Persist postISRCCD?",
        default=True,
    )

    # Calculate image quality statistics?
    doStandardStatistics = pexConfig.Field(
        dtype=bool,
        doc="Should standard image quality statistics be calculated?",
        default=True,
    )
    # Calculate additional statistics?
    doCalculateStatistics = pexConfig.Field(
        dtype=bool,
        doc="Should additional ISR statistics be calculated?",
        default=False,
    )
    isrStats = pexConfig.ConfigurableField(
        target=IsrStatisticsTask,
        doc="Task to calculate additional statistics.",
    )


class IsrTaskLSST(pipeBase.PipelineTask):
    ConfigClass = IsrTaskLSSTConfig
    _DefaultName = "isr"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.makeSubtask("overscan")
        self.makeSubtask("assembleCcd")
        self.makeSubtask("deferredChargeCorrection")
        self.makeSubtask("crosstalk")
        self.makeSubtask("masking")
        self.makeSubtask("isrStats")

    def runQuantum(self, butlerQC, inputRefs, outputRefs):

        inputs = butlerQC.get(inputRefs)
        self.validateInput(inputs)
        super().runQuantum(butlerQC, inputRefs, outputRefs)

    def validateInput(self, inputs):
        """
        This is a check that all the inputs required by the config
        are available.
        """

        inputMap = {'dnlLUT': self.config.doDiffNonLinearCorrection,
                    'bias': self.config.doBias,
                    'deferredChargeCalib': self.config.doDeferredCharge,
                    'linearizer': self.config.doLinearize,
                    'ptc': self.config.doGainNormalize,
                    'crosstalk': self.config.doCrosstalk,
                    'defects': self.config.doDefect,
                    'bfKernel': self.config.doBrighterFatter,
                    'dark': self.config.doDark,
                    }

        for calibrationFile, configValue in inputMap.items():
            if configValue and inputs[calibrationFile] is None:
                raise RuntimeError("Must supply ", calibrationFile)

    def diffNonLinearCorrection(self, ccdExposure, dnlLUT, **kwargs):
        # TODO DM 36636
        # isrFunctions.diffNonLinearCorrection
        pass

    def overscanCorrection(self, ccd, ccdExposure):
        # TODO DM 36637 for per amp

        overscans = []
        for amp in ccd:

            # Overscan correction on amp-by-amp basis.
            if amp.getRawHorizontalOverscanBBox().isEmpty():
                self.log.info("ISR_OSCAN: No overscan region.  Not performing overscan correction.")
                overscans.append(None)
            else:

                # Perform overscan correction on subregions.
                overscanResults = self.overscan.run(ccdExposure, amp)

                self.log.debug("Corrected overscan for amplifier %s.", amp.getName())
                if len(overscans) == 0:
                    ccdExposure.getMetadata().set('OVERSCAN', "Overscan corrected")

                overscans.append(overscanResults if overscanResults is not None else None)

        return overscans

    def snapCombine(self,**kwargs):
        #TODO DM 36638
        pass

    def getLinearizer(self, detector):
        # Here we assume linearizer as dict or LUT are not supported
        # TODO DM 28741

        # TODO construct isrcalib input
        linearizer = linearize.Linearizer(detector=detector, log=self.log)
        self.log.warning("Constructing linearizer from cameraGeom information.")

        return linearizer

    def gainNormalize(self, **kwargs):
        # TODO DM 36639
        gains = []
        readNoise = []

        return gains, readNoise

    def updateVariance(self, ampExposure, amp, overscanImage=None, ptcDataset=None):
        """Set the variance plane using the gain and read noise

        The read noise is calculated from the ``overscanImage`` if the
        ``doEmpiricalReadNoise`` option is set in the configuration; otherwise
        the value from the amplifier data is used.

        Parameters
        ----------
        ampExposure : `lsst.afw.image.Exposure`
            Exposure to process.
        amp : `lsst.afw.cameraGeom.Amplifier` or `FakeAmp`
            Amplifier detector data.
        overscanImage : `lsst.afw.image.MaskedImage`, optional.
            Image of overscan, required only for empirical read noise.
        ptcDataset : `lsst.ip.isr.PhotonTransferCurveDataset`, optional
            PTC dataset containing the gains and read noise.

        Raises
        ------
        RuntimeError
            Raised if either ``usePtcGains`` or ``usePtcReadNoise``
            are ``True``, but ptcDataset is not provided.

            Raised if ```doEmpiricalReadNoise`` is ``True`` but
            ``overscanImage`` is ``None``.

        See also
        --------
        lsst.ip.isr.isrFunctions.updateVariance
        """
        maskPlanes = [self.config.saturatedMaskName, self.config.suspectMaskName]
        if self.config.usePtcGains:
            if ptcDataset is None:
                raise RuntimeError("No ptcDataset provided to use PTC gains.")
            else:
                gain = ptcDataset.gain[amp.getName()]
                self.log.info("Using gain from Photon Transfer Curve.")
        else:
            gain = amp.getGain()

        if math.isnan(gain):
            gain = 1.0
            self.log.warning("Gain set to NAN!  Updating to 1.0 to generate Poisson variance.")
        elif gain <= 0:
            patchedGain = 1.0
            self.log.warning("Gain for amp %s == %g <= 0; setting to %f.",
                             amp.getName(), gain, patchedGain)
            gain = patchedGain

        if self.config.doEmpiricalReadNoise and overscanImage is None:
            badPixels = isrFunctions.countMaskedPixels(ampExposure.getMaskedImage(),
                                                       [self.config.saturatedMaskName,
                                                        self.config.suspectMaskName,
                                                        "BAD", "NO_DATA"])
            allPixels = ampExposure.getWidth() * ampExposure.getHeight()
            if allPixels == badPixels:
                # If the image is bad, do not raise.
                self.log.info("Skipping empirical read noise for amp %s.  No good pixels.",
                              amp.getName())
            else:
                raise RuntimeError("Overscan is none for EmpiricalReadNoise.")

        if self.config.doEmpiricalReadNoise and overscanImage is not None:
            stats = afwMath.StatisticsControl()
            stats.setAndMask(overscanImage.mask.getPlaneBitMask(maskPlanes))
            readNoise = afwMath.makeStatistics(overscanImage.getImage(),
                                               afwMath.STDEVCLIP, stats).getValue()
            self.log.info("Calculated empirical read noise for amp %s: %f.",
                          amp.getName(), readNoise)
        elif self.config.usePtcReadNoise:
            if ptcDataset is None:
                raise RuntimeError("No ptcDataset provided to use PTC readnoise.")
            else:
                readNoise = ptcDataset.noise[amp.getName()]
            self.log.info("Using read noise from Photon Transfer Curve.")
        else:
            readNoise = amp.getReadNoise()

        metadata = ampExposure.getMetadata()
        metadata[f'LSST GAIN {amp.getName()}'] = gain
        metadata[f'LSST READNOISE {amp.getName()}'] = readNoise

        isrFunctions.updateVariance(
            maskedImage=ampExposure.getMaskedImage(),
            gain=gain,
            readNoise=readNoise,
        )

    def maskNegativeVariance(self, exposure):
        """Identify and mask pixels with negative variance values.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to process.

        See Also
        --------
        lsst.ip.isr.isrFunctions.updateVariance
        """
        maskPlane = exposure.getMask().getPlaneBitMask(self.config.negativeVarianceMaskName)
        bad = numpy.where(exposure.getVariance().getArray() <= 0.0)
        exposure.mask.array[bad] |= maskPlane

    # TODO check make stats is necessary or not
    def variancePlane(self, ccdExposure, ccd, overscans, ptc):
        for amp, overscanResults in zip(ccd, overscans):
            if ccdExposure.getBBox().contains(amp.getBBox()):
                self.log.debug("Constructing variance map for amplifer %s.", amp.getName())
                ampExposure = ccdExposure.Factory(ccdExposure, amp.getBBox())

                isrFunctions.updateVariance(
                        maskedImage=ampExposure.getMaskedImage(),
                        gain=gain,
                        readNoise=readNoise,
                        )

#                if self.config.maskNegativeVariance:
#                    self.maskNegativeVariance(ccdExposure)

    def run(self,**kwargs):

        ccd = ccdExposure.getDetector()
        filterLabel = ccdExposure.getFilter()
        physicalFilter = isrFunctions.getPhysicalFilter(filterLabel, self.log)

        self.validateInput(**kwargs)
        if self.config.doDiffNonLinearCorrection:
            self.diffNonLinearCorrection(ccdExposure,dnlLUT,**kwargs)
        if self.config.doOverscan:
            overscans = self.overscanCorrection(self, ccd, ccdExposure)

        if self.config.doAssembleCcd:
            self.log.info("Assembling CCD from amplifiers.")
            ccdExposure = self.assembleCcd.assembleCcd(ccdExposure)

            if self.config.expectWcs and not ccdExposure.getWcs():
                self.log.warning("No WCS found in input exposure.")
            self.debugView(ccdExposure, "doAssembleCcd")

        if self.config.doSnapCombine:
            self.snapCombine(**kwargs)

        if self.config.doBias:
            self.log.info("Applying bias correction.")
            isrFunctions.biasCorrection(ccdExposure.getMaskedImage(), bias.getMaskedImage(),
                                        trimToFit=self.config.doTrimToMatchCalib)
            self.debugView(ccdExposure, "doBias")

        if self.config.doDeferredCharge:
            self.log.info("Applying deferred charge/CTI correction.")
            self.deferredChargeCorrection.run(ccdExposure, deferredChargeCalib)
            self.debugView(ccdExposure, "doDeferredCharge")


        if self.config.doLinearize:
            self.log.info("Applying linearizer.")
            linearizer.applyLinearity(image=ccdExposure.getMaskedImage().getImage(),
                                      detector=ccd, log=self.log)

        if self.config.doGainNormalize:
            gains, readNoise = self.gainNormalize(**kwargs)

        if self.config.doVariance:
            self.variancePlane(gains, readNoise, **kwargs)

        if self.config.doCrosstalk:
            self.log.info("Applying crosstalk correction.")
            self.crosstalk.run(ccdExposure, crosstalk=crosstalk,
                               crosstalkSources=crosstalkSources, isTrimmed=True)
            self.debugView(ccdExposure, "doCrosstalk")











