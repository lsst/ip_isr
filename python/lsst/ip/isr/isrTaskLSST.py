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
import lsst.afw.image as afwImage
import lsst.pipe.base.connectionTypes as cT
from lsst.meas.algorithms.detection import SourceDetectionTask

from .overscan import SerialOverscanCorrectionTask, ParallelOverscanCorrectionTask
from .assembleCcdTask import AssembleCcdTask
from .deferredCharge import DeferredChargeTask
from .crosstalk import CrosstalkTask
from .masking import MaskingTask
from .isrStatistics import IsrStatisticsTask
from .isr import maskNans


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
        doc="Complete kernel + gain solutions.",
        storageClass="BrighterFatterKernel",
        dimensions=["instrument", "detector"],
        isCalibration=True,
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
    outputBin1Exposure = cT.Output(
        name="postIsrBin1",
        doc="First binned image.",
        storageClass="ExposureF",
        dimensions=["instrument", "exposure", "detector"],
    )
    outputBin2Exposure = cT.Output(
        name="postIsrBin2",
        doc="Second binned image.",
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
            del self.dnlLUT
        if config.doBias is not True:
            del self.bias
        if config.doDeferredCharge is not True:
            del self.deferredChargeCalib
        if config.doLinearize is not True:
            del self.linearizer
        if not config.doCrosstalk and not config.doParallelOverscanCrosstalk:
            del self.crosstalk
        if config.doDefect is not True:
            del self.defects
        if config.doBrighterFatter is not True:
            del self.bfKernel
        if config.doDark is not True:
            del self.dark

        if config.doBinnedExposures is not True:
            del self.outputBin1Exposure
            del self.outputBin2Exposure
        if config.doSaveInterpPixels is not True:
            del self.preInterpExposure

        if config.doCalculateStatistics is not True:
            del self.outputStatistics


class IsrTaskLSSTConfig(pipeBase.PipelineTaskConfig,
                        pipelineConnections=IsrTaskLSSTConnections):
    """Configuration parameters for IsrTaskLSST.

    Items are grouped in the order in which they are executed by the task.
    """
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
        default=False,
    )

    doSerialOverscan = pexConfig.Field(
        dtype=bool,
        doc="Do serial overscan subtraction?",
        default=True,
    )
    serialOverscan = pexConfig.ConfigurableField(
        target=SerialOverscanCorrectionTask,
        doc="Serial overscan subtraction task for image segments.",
    )
    doParallelOverscan = pexConfig.Field(
        dtype=bool,
        doc="Do parallel overscan subtraction?",
        default=True,
    )
    parallelOverscan = pexConfig.ConfigurableField(
        target=ParallelOverscanCorrectionTask,
        doc="Parallel overscan subtraction task for image segments.",
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
        default=True,
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
        default=True,
    )
    doParallelOverscanCrosstalk = pexConfig.Field(
        dtype=bool,
        doc="Apply crosstalk correction in parallel overscan region?",
        default=True,
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
        default=True,
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
    brighterFatterFwhmForInterpolation = pexConfig.Field(
        dtype=float,
        doc="FWHM of PSF in arcseconds used for interpolation in brighter-fatter correction "
        "(currently unused).",
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
        default=True,
    )
    isrStats = pexConfig.ConfigurableField(
        target=IsrStatisticsTask,
        doc="Task to calculate additional statistics.",
    )

    # Make binned images?
    doBinnedExposures = pexConfig.Field(
        dtype=bool,
        doc="Should binned exposures be calculated?",
        default=False,
    )
    binFactor1 = pexConfig.Field(
        dtype=int,
        doc="Binning factor for first binned exposure. This is intended for a finely binned output.",
        default=8,
        check=lambda x: x > 1,
    )
    binFactor2 = pexConfig.Field(
        dtype=int,
        doc="Binning factor for second binned exposure. This is intended for a coarsely binned output.",
        default=64,
        check=lambda x: x > 1,
    )

    def validate(self):
        super().validate()

        if self.doCalculateStatistics and self.isrStats.doCtiStatistics:
            # DM-41912: Implement doApplyGains in LSST IsrTask
            # if self.doApplyGains !=
            #      self.isrStats.doApplyGainsForCtiStatistics:
            raise ValueError("doApplyGains must match isrStats.applyGainForCtiStatistics.")

    def setDefaults(self):
        super().setDefaults()

        self.serialOverscan.fitType = "MEDIAN_PER_ROW"
        self.parallelOverscan.fitType = "MEDIAN_PER_ROW"


class IsrTaskLSST(pipeBase.PipelineTask):
    ConfigClass = IsrTaskLSSTConfig
    _DefaultName = "isr"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.makeSubtask("serialOverscan")
        self.makeSubtask("parallelOverscan")
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
                    'crosstalk': self.config.doCrosstalk or self.config.doParallelOverscanCrosstalk,
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

    def maskAmplifiers(self, ccdExposure, detector, defects):
        """
        Check for fully masked bad amplifiers, and generate masks
        for SUSPECT and SATURATED values.

        Parameters
        ----------
        ccdExposure : `lsst.afw.image.Exposure`
            Input exposure to be masked.
        detector : `lsst.afw.cameraGeom.Detector`
            Detector object.
        defects : `lsst.ip.isr.Defects`
            List of defects.  Used to determine if an entire
            amplifier is bad.

        Returns
        -------
        badAmpDict : `str`[`bool`]
            Dictionary of amplifiers, keyed by name.
        """
        badAmpDict = {}

        maskedImage = ccdExposure.getMaskedImage()

        for amp in detector:
            ampName = amp.getName()
            badAmpDict[ampName] = False

            # Check if entire amp region is defined as a defect
            # NB: need to use amp.getBBox() for correct comparison with current
            # defects definition.
            if defects is not None:
                badAmpDict[ampName] = bool(sum([v.getBBox().contains(amp.getBBox()) for v in defects]))

            # In the case of a bad amp, we will set mask to "BAD"
            # (here use amp.getRawBBox() for correct association with pixels in
            # current ccdExposure).
            if badAmpDict[ampName]:
                dataView = afwImage.MaskedImageF(maskedImage, amp.getRawBBox(),
                                                 afwImage.PARENT)
                maskView = dataView.getMask()
                maskView |= maskView.getPlaneBitMask("BAD")
                del maskView

                self.log.warning("Amplifier %s is bad (completely covered with defects)", ampName)

                continue

            # We only mask saturated and suspect pixels now. Full defect
            # masking must come after detector is assembled to allow for
            # defects that cross amplifier boundaries.
            limits = {}
            if self.config.doSaturation:
                limits.update({self.config.saturatedMaskName: amp.getSaturation()})
            if self.config.doSuspect:
                limits.update({self.config.suspectMaskName: amp.getSuspectLevel()})
            if math.isfinite(self.config.saturation):
                limits.update({self.config.saturatedMaskName: self.config.saturation})

            for maskName, maskThreshold in limits.items():
                if not math.isnan(maskThreshold):
                    dataView = maskedImage.Factory(maskedImage, amp.getRawBBox())
                    isrFunctions.makeThresholdMask(
                        maskedImage=dataView,
                        threshold=maskThreshold,
                        growFootprints=0,
                        maskName=maskName
                    )

            # Determine if we've fully masked this amplifier with SUSPECT and
            # SAT pixels.
            maskView = afwImage.Mask(maskedImage.getMask(), amp.getRawDataBBox(),
                                     afwImage.PARENT)
            maskVal = maskView.getPlaneBitMask([self.config.saturatedMaskName,
                                                self.config.suspectMaskName])
            if numpy.all(maskView.getArray() & maskVal > 0):
                self.log.warning("Amplifier %s is bad (completely SATURATED or SUSPECT)", ampName)
                badAmpDict[ampName] = True
                maskView |= maskView.getPlaneBitMask("BAD")

        return badAmpDict

    def overscanCorrection(self, mode, detector, badAmpDict, ccdExposure):
        """Apply serial overscan correction in place to all amps.

        The actual overscan subtraction is performed by the
        `lsst.ip.isr.overscan.OverscanTask`, which is called here.

        Parameters
        ----------
        mode : `str`
            Must be `SERIAL` or `PARALLEL`.
        detector : `lsst.afw.cameraGeom.Detector`
            Detector object.
        badAmpDict : `dict`
            Dictionary of amp name to whether it is a bad amp.
        ccdExposure : `lsst.afw.image.Exposure`
            Exposure to have overscan correction performed.

        Returns
        -------
        serialOverscans : `list` [`lsst.pipe.base.Struct` or None]
            Each result struct has components:

            ``imageFit``
                Value or fit subtracted from the amplifier image data.
                (scalar or `lsst.afw.image.Image`)
            ``overscanFit``
                Value or fit subtracted from the overscan image data.
                (scalar or `lsst.afw.image.Image`)
            ``overscanImage``
                Image of the overscan region with the overscan
                correction applied. This quantity is used to estimate
                the amplifier read noise empirically.
                (`lsst.afw.image.Image`)
            ``overscanMean``
                Mean overscan fit value. (`float`)
            ``overscanMedian``
                Median overscan fit value. (`float`)
            ``overscanSigma``
                Clipped standard deviation of the overscan fit. (`float`)
            ``residualMean``
                Mean of the overscan after fit subtraction. (`float`)
            ``residualMedian``
                Median of the overscan after fit subtraction. (`float`)
            ``residualSigma``
                Clipped standard deviation of the overscan after fit
                subtraction. (`float`)

        See Also
        --------
        lsst.ip.isr.overscan.OverscanTask
        """
        if mode not in ["SERIAL", "PARALLEL"]:
            raise ValueError("Mode must be SERIAL or PARALLEL")

        overscans = []

        for i, amp in enumerate(detector):
            ampName = amp.getName()
            if badAmpDict[ampName] or not ccdExposure.getBBox().contains(amp.getBBox()):
                results = None
            else:
                # Question: should this be something else when in PARALLEL
                # mode?
                if amp.getRawHorizontalOverscanBBox().isEmpty():
                    self.log.warning(
                        "ISR_OSCAN: No overscan region for amp %s. Not performing overscan correction.",
                        ampName,
                    )
                    results = None
                else:
                    if mode == "SERIAL":
                        results = self.serialOverscan.run(ccdExposure, amp)
                    else:
                        results = self.parallelOverscan.run(ccdExposure, amp)

                    metadata = ccdExposure.getMetadata()
                    keyBase = "LSST ISR OVERSCAN"
                    metadata[f"{keyBase} {mode} MEAN {ampName}"] = results.overscanMean
                    metadata[f"{keyBase} {mode} MEDIAN {ampName}"] = results.overscanMedian
                    metadata[f"{keyBase} {mode} STDEV {ampName}"] = results.overscanSigma

                    metadata[f"{keyBase} RESIDUAL {mode} MEAN {ampName}"] = results.residualMean
                    metadata[f"{keyBase} RESIDUAL {mode} MEDIAN {ampName}"] = results.residualMedian
                    metadata[f"{keyBase} RESIDUAL {mode} STDEV {ampName}"] = results.residualSigma

            overscans[i] = results

        # Question: should this be finer grained?
        ccdExposure.getMetadata().set("OVERSCAN", "Overscan corrected")

        return overscans

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

    def updateVariance(self, ampExposure, amp, ptcDataset=None):
        """Set the variance plane using the gain and read noise.

        Parameters
        ----------
        ampExposure : `lsst.afw.image.Exposure`
            Exposure to process.
        amp : `lsst.afw.cameraGeom.Amplifier` or `FakeAmp`
            Amplifier detector data.
        ptcDataset : `lsst.ip.isr.PhotonTransferCurveDataset`, optional
            PTC dataset containing the gains and read noise.

        Raises
        ------
        RuntimeError
            Raised if ptcDataset is not provided.

        See also
        --------
        lsst.ip.isr.isrFunctions.updateVariance
        """
        # Get gains from PTC
        if ptcDataset is None:
            raise RuntimeError("No ptcDataset provided to use PTC gains.")
        else:
            gain = ptcDataset.gain[amp.getName()]
            self.log.debug("Getting gain from Photon Transfer Curve.")

        if math.isnan(gain):
            gain = 1.0
            self.log.warning("Gain set to NAN!  Updating to 1.0 to generate Poisson variance.")
        elif gain <= 0:
            patchedGain = 1.0
            self.log.warning("Gain for amp %s == %g <= 0; setting to %f.",
                             amp.getName(), gain, patchedGain)
            gain = patchedGain

        # Get read noise from PTC
        if ptcDataset is None:
            raise RuntimeError("No ptcDataset provided to use PTC readnoise.")
        else:
            readNoise = ptcDataset.noise[amp.getName()]
            self.log.debug("Getting read noise from Photon Transfer Curve.")

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

                self.updateVariance(ampExposure, amp, ptcDataset=ptc)

                if self.config.qa is not None and self.config.qa.saveStats is True:
                    qaStats = afwMath.makeStatistics(ampExposure.getVariance(),
                                                     afwMath.MEDIAN | afwMath.STDEVCLIP)
                    self.log.debug("  Variance stats for amplifer %s: %f +/- %f.",
                                   amp.getName(), qaStats.getValue(afwMath.MEDIAN),
                                   qaStats.getValue(afwMath.STDEVCLIP))
        if self.config.maskNegativeVariance:
            self.maskNegativeVariance(ccdExposure)

    def maskDefect(self, exposure, defectBaseList):
        """Mask defects using mask plane "BAD", in place.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to process.

        defectBaseList : defect-type
            List of defects to mask. Can be of type  `lsst.ip.isr.Defects`
            or `list` of `lsst.afw.image.DefectBase`.
        """
        maskedImage = exposure.getMaskedImage()
        if not isinstance(defectBaseList, Defects):
            # Promotes DefectBase to Defect
            defectList = Defects(defectBaseList)
        else:
            defectList = defectBaseList
        defectList.maskPixels(maskedImage, maskName="BAD")

    def maskEdges(self, exposure, numEdgePixels=0, maskPlane="SUSPECT", level='DETECTOR'):
        """Mask edge pixels with applicable mask plane.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to process.
        numEdgePixels : `int`, optional
            Number of edge pixels to mask.
        maskPlane : `str`, optional
            Mask plane name to use.
        level : `str`, optional
            Level at which to mask edges.
        """
        maskedImage = exposure.getMaskedImage()
        maskBitMask = maskedImage.getMask().getPlaneBitMask(maskPlane)

        if numEdgePixels > 0:
            if level == 'DETECTOR':
                boxes = [maskedImage.getBBox()]
            elif level == 'AMP':
                boxes = [amp.getBBox() for amp in exposure.getDetector()]

            for box in boxes:
                # This makes a bbox numEdgeSuspect pixels smaller than the
                # image on each side
                subImage = maskedImage[box]
                box.grow(-numEdgePixels)
                # Mask pixels outside box
                SourceDetectionTask.setEdgeBits(
                    subImage,
                    box,
                    maskBitMask)

    def maskNan(self, exposure):
        """Mask NaNs using mask plane "UNMASKEDNAN", in place.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to process.

        Notes
        -----
        We mask over all non-finite values (NaN, inf), including those
        that are masked with other bits (because those may or may not be
        interpolated over later, and we want to remove all NaN/infs).
        Despite this behaviour, the "UNMASKEDNAN" mask plane is used to
        preserve the historical name.
        """
        maskedImage = exposure.getMaskedImage()

        # Find and mask NaNs
        maskedImage.getMask().addMaskPlane("UNMASKEDNAN")
        maskVal = maskedImage.getMask().getPlaneBitMask("UNMASKEDNAN")
        numNans = maskNans(maskedImage, maskVal)
        self.metadata["NUMNANS"] = numNans
        if numNans > 0:
            self.log.warning("There were %d unmasked NaNs.", numNans)

    def countBadPixels(self, exposure):
        """
        Notes
        -----
        Reset and interpolate bad pixels.

        Large contiguous bad regions (which should have the BAD mask
        bit set) should have their values set to the image median.
        This group should include defects and bad amplifiers. As the
        area covered by these defects are large, there's little
        reason to expect that interpolation would provide a more
        useful value.

        Smaller defects can be safely interpolated after the larger
        regions have had their pixel values reset.  This ensures
        that the remaining defects adjacent to bad amplifiers (as an
        example) do not attempt to interpolate extreme values.
        """
        badPixelCount, badPixelValue = isrFunctions.setBadRegions(exposure)
        if badPixelCount > 0:
            self.log.info("Set %d BAD pixels to %f.", badPixelCount, badPixelValue)

    @contextmanager
    def flatContext(self, exp, flat, dark=None):
        """Context manager that applies and removes flats and darks,
        if the task is configured to apply them.

        Parameters
        ----------
        exp : `lsst.afw.image.Exposure`
            Exposure to process.
        flat : `lsst.afw.image.Exposure`
            Flat exposure the same size as ``exp``.
        dark : `lsst.afw.image.Exposure`, optional
            Dark exposure the same size as ``exp``.

        Yields
        ------
        exp : `lsst.afw.image.Exposure`
            The flat and dark corrected exposure.
        """
        if self.config.doDark and dark is not None:
            self.darkCorrection(exp, dark)
        if self.config.doFlat:
            self.flatCorrection(exp, flat)
        try:
            yield exp
        finally:
            if self.config.doFlat:
                self.flatCorrection(exp, flat, invert=True)
            if self.config.doDark and dark is not None:
                self.darkCorrection(exp, dark, invert=True)

    def getBrighterFatterKernel(self, detector, bfKernel):
        detName = detector.getName()

        # This is expected to be a dictionary of amp-wise gains.
        bfGains = bfKernel.gain
        if bfKernel.level == 'DETECTOR':
            if detName in bfKernel.detKernels:
                bfKernelOut = bfKernel.detKernels[detName]
                return bfKernelOut, bfGains
            else:
                raise RuntimeError("Failed to extract kernel from new-style BF kernel.")
        elif bfKernel.level == 'AMP':
            self.log.warning("Making DETECTOR level kernel from AMP based brighter "
                             "fatter kernels.")
            bfKernel.makeDetectorKernelFromAmpwiseKernels(detName)
            bfKernelOut = bfKernel.detKernels[detName]
            return bfKernelOut, bfGains

    def applyBrighterFatterCorrection(self, ccdExposure, flat, dark, bfKernel, bfGains):
        # We need to apply flats and darks before we can interpolate, and
        # we need to interpolate before we do B-F, but we do B-F without
        # the flats and darks applied so we can work in units of electrons
        # or holes. This context manager applies and then removes the darks
        # and flats.
        #
        # We also do not want to interpolate values here, so operate on
        # temporary images so we can apply only the BF-correction and roll
        # back the interpolation.
        # This won't be necessary once the gain normalization
        # is done appropriately.
        interpExp = ccdExposure.clone()
        with self.flatContext(interpExp, flat, dark):
            isrFunctions.interpolateFromMask(
                maskedImage=interpExp.getMaskedImage(),
                fwhm=self.config.brighterFatterFwhmForInterpolation,
                growSaturatedFootprints=self.config.growSaturationFootprintSize,
                maskNameList=list(self.config.brighterFatterMaskListToInterpolate)
            )
        bfExp = interpExp.clone()
        self.log.info("Applying brighter-fatter correction using kernel type %s / gains %s.",
                      type(bfKernel), type(bfGains))
        bfResults = isrFunctions.brighterFatterCorrection(bfExp, bfKernel,
                                                          self.config.brighterFatterMaxIter,
                                                          self.config.brighterFatterThreshold,
                                                          self.config.brighterFatterApplyGain,
                                                          bfGains)
        if bfResults[1] == self.config.brighterFatterMaxIter:
            self.log.warning("Brighter-fatter correction did not converge, final difference %f.",
                             bfResults[0])
        else:
            self.log.info("Finished brighter-fatter correction in %d iterations.",
                          bfResults[1])

        image = ccdExposure.getMaskedImage().getImage()
        bfCorr = bfExp.getMaskedImage().getImage()
        bfCorr -= interpExp.getMaskedImage().getImage()
        image += bfCorr

        # Applying the brighter-fatter correction applies a
        # convolution to the science image. At the edges this
        # convolution may not have sufficient valid pixels to
        # produce a valid correction. Mark pixels within the size
        # of the brighter-fatter kernel as EDGE to warn of this
        # fact.
        self.log.info("Ensuring image edges are masked as EDGE to the brighter-fatter kernel size.")
        self.maskEdges(ccdExposure, numEdgePixels=numpy.max(bfKernel.shape) // 2,
                       maskPlane="EDGE")

        if self.config.brighterFatterMaskGrowSize > 0:
            self.log.info("Growing masks to account for brighter-fatter kernel convolution.")
            for maskPlane in self.config.brighterFatterMaskListToInterpolate:
                isrFunctions.growMasks(ccdExposure.getMask(),
                                       radius=self.config.brighterFatterMaskGrowSize,
                                       maskNameList=maskPlane,
                                       maskValue=maskPlane)

        return ccdExposure

    def darkCorrection(self, exposure, darkExposure, invert=False):
        """Apply dark correction in place.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to process.
        darkExposure : `lsst.afw.image.Exposure`
            Dark exposure of the same size as ``exposure``.
        invert : `Bool`, optional
            If True, re-add the dark to an already corrected image.

        Raises
        ------
        RuntimeError
            Raised if either ``exposure`` or ``darkExposure`` do not
            have their dark time defined.

        See Also
        --------
        lsst.ip.isr.isrFunctions.darkCorrection
        """
        expScale = exposure.getInfo().getVisitInfo().getDarkTime()
        if math.isnan(expScale):
            raise RuntimeError("Exposure darktime is NAN.")
        if darkExposure.getInfo().getVisitInfo() is not None \
                and not math.isnan(darkExposure.getInfo().getVisitInfo().getDarkTime()):
            darkScale = darkExposure.getInfo().getVisitInfo().getDarkTime()
        else:
            # DM-17444: darkExposure.getInfo.getVisitInfo() is None
            #           so getDarkTime() does not exist.
            self.log.warning("darkExposure.getInfo().getVisitInfo() does not exist. Using darkScale = 1.0.")
            darkScale = 1.0

        isrFunctions.darkCorrection(
            maskedImage=exposure.getMaskedImage(),
            darkMaskedImage=darkExposure.getMaskedImage(),
            expScale=expScale,
            darkScale=darkScale,
            invert=invert,
        )

    @staticmethod
    def extractCalibDate(calib):
        """Extract common calibration metadata values that will be written to
        output header.

        Parameters
        ----------
        calib : `lsst.afw.image.Exposure` or `lsst.ip.isr.IsrCalib`
            Calibration to pull date information from.

        Returns
        -------
        dateString : `str`
            Calibration creation date string to add to header.
        """
        if hasattr(calib, "getMetadata"):
            if 'CALIB_CREATION_DATE' in calib.getMetadata():
                return " ".join((calib.getMetadata().get("CALIB_CREATION_DATE", "Unknown"),
                                 calib.getMetadata().get("CALIB_CREATION_TIME", "Unknown")))
            else:
                return " ".join((calib.getMetadata().get("CALIB_CREATE_DATE", "Unknown"),
                                 calib.getMetadata().get("CALIB_CREATE_TIME", "Unknown")))
        else:
            return "Unknown Unknown"

    def doLinearize(self, detector):
        """Check if linearization is needed for the detector cameraGeom.

        Checks config.doLinearize and the linearity type of the first
        amplifier.

        Parameters
        ----------
        detector : `lsst.afw.cameraGeom.Detector`
            Detector to get linearity type from.

        Returns
        -------
        doLinearize : `Bool`
            If True, linearization should be performed.
        """
        return self.config.doLinearize and \
            detector.getAmplifiers()[0].getLinearityType() != NullLinearityType

    def makeBinnedImages(self, exposure):
        """Make visualizeVisit style binned exposures.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to bin.

        Returns
        -------
        bin1 : `lsst.afw.image.Exposure`
            Binned exposure using binFactor1.
        bin2 : `lsst.afw.image.Exposure`
            Binned exposure using binFactor2.
        """
        mi = exposure.getMaskedImage()

        bin1 = afwMath.binImage(mi, self.config.binFactor1)
        bin2 = afwMath.binImage(mi, self.config.binFactor2)

        return bin1, bin2

    def run(self, *, ccdExposure, dnlLUT=None, bias=None, deferredChargeCalib=None, linearizer=None,
            ptc=None, crosstalk=None, defects=None, bfKernel=None, bfGains=None, dark=None,
            flat=None, camera=None, **kwargs
            ):

        detector = ccdExposure.getDetector()

        if self.config.doHeaderProvenance:
            # Inputs have been validated, so we can add their date
            # information to the output header.
            exposureMetadata = ccdExposure.getMetadata()
            exposureMetadata["LSST CALIB DATE PTC"] = self.extractCalibDate(ptc)
            if self.config.doDiffNonLinearCorrection:
                exposureMetadata["LSST CALIB DATE DNL"] = self.extractCalibDate(dnlLUT)
            if self.config.doBias:
                exposureMetadata["LSST CALIB DATE BIAS"] = self.extractCalibDate(bias)
            if self.config.doDeferredCharge:
                exposureMetadata["LSST CALIB DATE CTI"] = self.extractCalibDate(deferredChargeCalib)
            if self.doLinearize(detector):
                exposureMetadata["LSST CALIB DATE LINEARIZER"] = self.extractCalibDate(linearizer)
            if self.config.doCrosstalk or self.config.doParallelOverscanCrosstalk:
                exposureMetadata["LSST CALIB DATE CROSSTALK"] = self.extractCalibDate(crosstalk)
            if self.config.doDefect:
                exposureMetadata["LSST CALIB DATE DEFECTS"] = self.extractCalibDate(defects)
            if self.config.doBrighterFatter:
                exposureMetadata["LSST CALIB DATE BFK"] = self.extractCalibDate(bfKernel)
            if self.config.doDark:
                exposureMetadata["LSST CALIB DATE DARK"] = self.extractCalibDate(dark)

        # First we must mark which amplifiers are bad.
        # TODO: Should this be here or after overscan to ensure
        # correct units of the saturation?
        badAmpDict = self.maskAmplifiers(ccdExposure, detector, defects)

        if self.config.doDiffNonLinearCorrection:
            self.diffNonLinearCorrection(ccdExposure, dnlLUT)

        if self.config.doSerialOverscan:
            # Input units: ADU
            serialOverscans = self.overscanCorrection("SERIAL", detector, badAmpDict, ccdExposure)

        if self.config.doParallelOverscanCrosstalk:
            # Input units: ADU
            # Make sure that the units here are consistent with later
            # application.
            self.crosstalk.run(
                ccdExposure,
                crosstalk=crosstalk,
                camera=camera,
                parallelOverscanRegion=True,
            )

        if self.config.doParallelOverscans:
            # Input units: ADU
            # At the moment we do not use the parallelOverscans return.
            _ = self.overscanCorrection("PARALLEL", detector, badAmpDict, ccdExposure)

        if self.config.doAssembleCcd:
            # Input units: ADU
            self.log.info("Assembling CCD from amplifiers.")
            ccdExposure = self.assembleCcd.assembleCcd(ccdExposure)

            if self.config.expectWcs and not ccdExposure.getWcs():
                self.log.warning("No WCS found in input exposure.")

        if self.config.doBias:
            # Input units: ADU
            self.log.info("Applying bias correction.")
            isrFunctions.biasCorrection(ccdExposure.getMaskedImage(), bias.getMaskedImage())

        if self.config.doDeferredCharge:
            # Input units: ADU
            self.log.info("Applying deferred charge/CTI correction.")
            self.deferredChargeCorrection.run(ccdExposure, deferredChargeCalib)

        if self.config.doLinearize:
            # Input units: ADU
            self.log.info("Applying linearizer.")
            linearizer = self.getLinearizer(detector=detector)
            linearizer.applyLinearity(image=ccdExposure.getMaskedImage().getImage(),
                                      detector=detector, log=self.log)

        if self.config.doGainNormalize:
            # Input units: ADU
            # Output units: electrons
            # TODO DM 36639
            gains, readNoise = self.gainNormalize(**kwargs)

        if self.config.doVariance:
            # Input units: electrons
            self.variancePlane(ccdExposure, detector, serialOverscans, ptc)

        if self.config.doCrosstalk:
            # Input units: electrons
            self.log.info("Applying crosstalk correction.")
            self.crosstalk.run(ccdExposure, crosstalk=crosstalk)

        # Masking block (defects, NAN pixels and trails).
        # Saturated and suspect pixels have already been masked.
        if self.config.doDefect:
            # Input units: electrons
            self.log.info("Applying defects masking.")
            self.maskDefect(ccdExposure, defects)

        if self.config.doNanMasking:
            self.log.info("Masking non-finite (NAN, inf) value pixels.")
            self.maskNan(ccdExposure)

        if self.config.doWidenSaturationTrails:
            self.log.info("Widening saturation trails.")
            isrFunctions.widenSaturationTrails(ccdExposure.getMaskedImage().getMask())

        preInterpExp = None
        if self.config.doSaveInterpPixels:
            preInterpExp = ccdExposure.clone()

        if self.config.doSetBadRegions:
            self.log.info('Counting pixels in BAD regions.')
            self.countBadPixels(ccdExposure)

        if self.config.doInterpolate:
            self.log.info("Interpolating masked pixels.")
            isrFunctions.interpolateFromMask(
                maskedImage=ccdExposure.getMaskedImage(),
                fwhm=self.config.brighterFatterFwhmForInterpolation,
                growSaturatedFootprints=self.config.growSaturationFootprintSize,
                maskNameList=list(self.config.maskListToInterpolate)
            )

        if self.config.doBrighterFatter:
            # Input units: electrons
            self.log.info("Applying Bright-Fatter kernels.")
            bfKernelOut, bfGains = self.getBrighterFatterKernel(detector, bfKernel)
            ccdExposure = self.applyBrighterFatterCorrection(ccdExposure, flat, dark, bfKernelOut, bfGains)

        if self.config.doDark:
            # Input units: electrons
            self.log.info("Applying dark subtraction.")
            self.darkCorrection(ccdExposure, dark)

        # Calculate standard image quality statistics
        if self.config.doStandardStatistics:
            metadata = ccdExposure.getMetadata()
            for amp in detector:
                ampExposure = ccdExposure.Factory(ccdExposure, amp.getBBox())
                ampName = amp.getName()
                metadata[f"LSST ISR MASK SAT {ampName}"] = isrFunctions.countMaskedPixels(
                    ampExposure.getMaskedImage(),
                    [self.config.saturatedMaskName]
                )
                metadata[f"LSST ISR MASK BAD {ampName}"] = isrFunctions.countMaskedPixels(
                    ampExposure.getMaskedImage(),
                    ["BAD"]
                )
                qaStats = afwMath.makeStatistics(ampExposure.getImage(),
                                                 afwMath.MEAN | afwMath.MEDIAN | afwMath.STDEVCLIP)

                metadata[f"LSST ISR FINAL MEAN {ampName}"] = qaStats.getValue(afwMath.MEAN)
                metadata[f"LSST ISR FINAL MEDIAN {ampName}"] = qaStats.getValue(afwMath.MEDIAN)
                metadata[f"LSST ISR FINAL STDEV {ampName}"] = qaStats.getValue(afwMath.STDEVCLIP)

                k1 = f"LSST ISR FINAL MEDIAN {ampName}"
                k2 = f"LSST ISR OVERSCAN SERIAL MEDIAN {ampName}"
                if self.config.doOverscan and k1 in metadata and k2 in metadata:
                    metadata[f"LSST ISR LEVEL {ampName}"] = metadata[k1] - metadata[k2]
                else:
                    metadata[f"LSST ISR LEVEL {ampName}"] = numpy.nan

        # calculate additional statistics.
        outputStatistics = None
        if self.config.doCalculateStatistics:
            outputStatistics = self.isrStats.run(ccdExposure, overscanResults=serialOverscans,
                                                 ptc=ptc).results

        # do image binning.
        outputBin1Exposure = None
        outputBin2Exposure = None
        if self.config.doBinnedExposures:
            outputBin1Exposure, outputBin2Exposure = self.makeBinnedImages(ccdExposure)

        return pipeBase.Struct(
            exposure=ccdExposure,

            outputBin1Exposure=outputBin1Exposure,
            outputBin2Exposure=outputBin2Exposure,

            preInterpExposure=preInterpExp,
            outputExposure=ccdExposure,
            outputStatistics=outputStatistics,
        )
