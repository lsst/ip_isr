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

from .ampOffset import AmpOffsetTask
from .overscan import SerialOverscanCorrectionTask, ParallelOverscanCorrectionTask
from .overscanAmpConfig import OverscanCameraConfig
from .assembleCcdTask import AssembleCcdTask
from .deferredCharge import DeferredChargeTask
from .crosstalk import CrosstalkTask
from .masking import MaskingTask
from .isrStatistics import IsrStatisticsTask
from .isr import maskNans
from .ptcDataset import PhotonTransferCurveDataset


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
        if not config.doCrosstalk and not config.overscanCamera.doAnyParallelOverscanCrosstalk:
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

    # Calib checking configuration:
    doRaiseOnCalibMismatch = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Should IsrTaskLSST halt if exposure and calibration header values do not match?",
    )
    cameraKeywordsToCompare = pexConfig.ListField(
        dtype=str,
        doc="List of header keywords to compare between exposure and calibrations.",
        default=[],
    )

    # Differential non-linearity correction.
    doDiffNonLinearCorrection = pexConfig.Field(
        dtype=bool,
        doc="Do differential non-linearity correction?",
        default=False,
    )

    integerDitherMode = pexConfig.ChoiceField(
        dtype=str,
        doc="Dithering mode to cancel integerization of counts.",
        default="SYMMETRIC",
        allowed={
            "POSITIVE": "Dithering is done with a uniform random in the range [0, 1).",
            "NEGATIVE": "Dithering is done with a uniform random in the range [-1, 0).",
            "SYMMETRIC": "Dithering is done with a uniform random in the range [-0.5, 0.5).",
            "NONE": "No dithering is performed.",
        },
    )

    nominalGain = pexConfig.Field(
        dtype=float,
        default=1.7,
        doc="Nominal gain to use if no PTC is supplied.",
    )

    overscanCamera = pexConfig.ConfigField(
        dtype=OverscanCameraConfig,
        doc="Per-detector and per-amplifier overscan configurations.",
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

    # Gains.
    doCorrectGains = pexConfig.Field(
        dtype=bool,
        doc="Apply temperature correction to the gains?",
        default=False,
    )
    doApplyGains = pexConfig.Field(
        dtype=bool,
        doc="Apply gains to the image?",
        default=True,
    )

    # Variance construction.
    doVariance = pexConfig.Field(
        dtype=bool,
        doc="Calculate variance?",
        default=True
    )
    # gain = pexConfig.Field(
    #     dtype=float,
    #     doc="The gain to use if no Detector is present in the Exposure
    #  (ignored if NaN).",
    #     default=float("NaN"),
    # )
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
    doSaturation = pexConfig.Field(
        dtype=bool,
        doc="Mask saturated pixels? NB: this is totally independent of the"
        " interpolation option - this is ONLY setting the bits in the mask."
        " To have them interpolated make sure doSaturationInterpolation=True",
        default=True,
    )
    saturation = pexConfig.Field(
        dtype=float,
        doc="The saturation level to use if no Detector is present in the Exposure (ignored if NaN)",
        default=float("NaN"),
    )
    saturatedMaskName = pexConfig.Field(
        dtype=str,
        doc="Name of mask plane to use in saturation detection and interpolation.",
        default="SAT",
    )
    doSuspect = pexConfig.Field(
        dtype=bool,
        doc="Mask suspect pixels?",
        default=False,
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

    # Amp offset correction.
    doAmpOffset = pexConfig.Field(
        doc="Calculate amp offset corrections?",
        dtype=bool,
        default=False,
    )
    ampOffset = pexConfig.ConfigurableField(
        doc="Amp offset correction task.",
        target=AmpOffsetTask,
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

    # Flat correction.
    doFlat = pexConfig.Field(
        dtype=bool,
        doc="Apply flat field correction.",
        default=True,
    )
    flatScalingType = pexConfig.ChoiceField(
        dtype=str,
        doc="The method for scaling the flat on the fly.",
        default='USER',
        allowed={
            "USER": "Scale by flatUserScale",
            "MEAN": "Scale by the inverse of the mean",
            "MEDIAN": "Scale by the inverse of the median",
        },
    )
    flatUserScale = pexConfig.Field(
        dtype=float,
        doc="If flatScalingType is 'USER' then scale flat by this amount; ignored otherwise.",
        default=1.0,
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
        if self.ampOffset.doApplyAmpOffset and not self.doAmpOffset:
            raise ValueError("ampOffset.doApplyAmpOffset requires doAmpOffset to be True.")

    def setDefaults(self):
        super().setDefaults()


class IsrTaskLSST(pipeBase.PipelineTask):
    ConfigClass = IsrTaskLSSTConfig
    _DefaultName = "isr"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.makeSubtask("assembleCcd")
        self.makeSubtask("deferredChargeCorrection")
        self.makeSubtask("crosstalk")
        self.makeSubtask("masking")
        self.makeSubtask("isrStats")

    def runQuantum(self, butlerQC, inputRefs, outputRefs):

        inputs = butlerQC.get(inputRefs)
        self.validateInput(inputs)

        if self.config.doHeaderProvenance:
            # Add calibration provenanace info to header.
            exposureMetadata = inputs['ccdExposure'].getMetadata()
            for inputName in sorted(list(inputs.keys())):
                reference = getattr(inputRefs, inputName, None)
                if reference is not None and hasattr(reference, "run"):
                    runKey = f"LSST CALIB RUN {inputName.upper()}"
                    runValue = reference.run
                    idKey = f"LSST CALIB UUID {inputName.upper()}"
                    idValue = str(reference.id)
                    dateKey = f"LSST CALIB DATE {inputName.upper()}"
                    dateValue = self.extractCalibDate(inputs[inputName])

                    exposureMetadata[runKey] = runValue
                    exposureMetadata[idKey] = idValue
                    exposureMetadata[dateKey] = dateValue

        super().runQuantum(butlerQC, inputRefs, outputRefs)

    def validateInput(self, inputs):
        """
        This is a check that all the inputs required by the config
        are available.
        """

        doCrosstalk = self.config.doCrosstalk or self.config.overscanCamera.doAnyParallelOverscanCrosstalk

        inputMap = {'dnlLUT': self.config.doDiffNonLinearCorrection,
                    'bias': self.config.doBias,
                    'deferredChargeCalib': self.config.doDeferredCharge,
                    'linearizer': self.config.doLinearize,
                    'ptc': self.config.doApplyGains,
                    'crosstalk': doCrosstalk,
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

    def maskFullDefectAmplifiers(self, ccdExposure, detector, defects):
        """
        Check for fully masked bad amplifiers and mask them.

        Full defect masking happens later to allow for defects which
        cross amplifier boundaries.

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
            Dictionary of amplifiers, keyed by name, value is True if
            amplifier is fully masked.
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

        return badAmpDict

    def maskSaturatedPixels(self, badAmpDict, ccdExposure, detector, ptcGains=None):
        """
        Mask SATURATED and SUSPECT pixels and check if any amplifiers
        are fully masked.

        Parameters
        ----------
        badAmpDict : `str` [`bool`]
            Dictionary of amplifiers, keyed by name, value is True if
            amplifier is fully masked.
        ccdExposure : `lsst.afw.image.Exposure`
            Input exposure to be masked.
        detector : `lsst.afw.cameraGeom.Detector`
            Detector object.
        defects : `lsst.ip.isr.Defects`
            List of defects.  Used to determine if an entire
            amplifier is bad.
        ptcGains : `dict`[`str`], optional
            Dictionary keyed by amp name containing the PTC gains.

        Returns
        -------
        badAmpDict : `str`[`bool`]
            Dictionary of amplifiers, keyed by name.
        """
        maskedImage = ccdExposure.getMaskedImage()

        for amp in detector:
            ampName = amp.getName()

            if badAmpDict[ampName]:
                # No need to check fully bad amplifiers.
                continue

            if ptcGains is not None and ampName in ptcGains:
                gain = ptcGains[ampName]
            else:
                gain = 1.0

            # Mask saturated and suspect pixels.
            limits = {}
            if self.config.doSaturation:
                # Set to the default from the camera model.
                limits.update({self.config.saturatedMaskName: amp.getSaturation() * gain})
                # And update if it is set in the config.
                if math.isfinite(self.config.saturation):
                    limits.update({self.config.saturatedMaskName: self.config.saturation})
            if self.config.doSuspect:
                limits.update({self.config.suspectMaskName: amp.getSuspectLevel()})

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

    def overscanCorrection(self, mode, detectorConfig, detector, badAmpDict, ccdExposure):
        """Apply serial overscan correction in place to all amps.

        The actual overscan subtraction is performed by the
        `lsst.ip.isr.overscan.OverscanTask`, which is called here.

        Parameters
        ----------
        mode : `str`
            Must be `SERIAL` or `PARALLEL`.
        detectorConfig : `lsst.ip.isr.OverscanDetectorConfig`
            Per-amplifier configurations.
        detector : `lsst.afw.cameraGeom.Detector`
            Detector object.
        badAmpDict : `dict`
            Dictionary of amp name to whether it is a bad amp.
        ccdExposure : `lsst.afw.image.Exposure`
            Exposure to have overscan correction performed.

        Returns
        -------
        overscans : `list` [`lsst.pipe.base.Struct` or None]
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

        # This returns a list in amp order, with None for uncorrected amps.
        overscans = []

        for i, amp in enumerate(detector):
            ampName = amp.getName()

            ampConfig = detectorConfig.getOverscanAmpConfig(amp)

            if mode == "SERIAL" and not ampConfig.doSerialOverscan:
                self.log.debug(
                    "ISR_OSCAN: Amplifier %s/%s configured to skip serial overscan.",
                    detector.getName(),
                    ampName,
                )
                results = None
            elif mode == "PARALLEL" and not ampConfig.doParallelOverscan:
                self.log.debug(
                    "ISR_OSCAN: Amplifier %s configured to skip parallel overscan.",
                    detector.getName(),
                    ampName,
                )
                results = None
            elif badAmpDict[ampName] or not ccdExposure.getBBox().contains(amp.getBBox()):
                results = None
            else:
                # This check is to confirm that we are not trying to run
                # overscan on an already trimmed image. Therefore, always
                # checking just the horizontal overscan bounding box is
                # sufficient.
                if amp.getRawHorizontalOverscanBBox().isEmpty():
                    self.log.warning(
                        "ISR_OSCAN: No overscan region for amp %s. Not performing overscan correction.",
                        ampName,
                    )
                    results = None
                else:
                    if mode == "SERIAL":
                        # We need to set up the subtask here with a custom
                        # configuration.
                        serialOverscan = SerialOverscanCorrectionTask(config=ampConfig.serialOverscanConfig)
                        results = serialOverscan.run(ccdExposure, amp)
                    else:
                        parallelOverscan = ParallelOverscanCorrectionTask(
                            config=ampConfig.parallelOverscanConfig,
                        )
                        results = parallelOverscan.run(ccdExposure, amp)

                    metadata = ccdExposure.getMetadata()
                    keyBase = "LSST ISR OVERSCAN"
                    metadata[f"{keyBase} {mode} MEAN {ampName}"] = results.overscanMean
                    metadata[f"{keyBase} {mode} MEDIAN {ampName}"] = results.overscanMedian
                    metadata[f"{keyBase} {mode} STDEV {ampName}"] = results.overscanSigma

                    metadata[f"{keyBase} RESIDUAL {mode} MEAN {ampName}"] = results.residualMean
                    metadata[f"{keyBase} RESIDUAL {mode} MEDIAN {ampName}"] = results.residualMedian
                    metadata[f"{keyBase} RESIDUAL {mode} STDEV {ampName}"] = results.residualSigma

            overscans.append(results)

        # Question: should this be finer grained?
        ccdExposure.getMetadata().set("OVERSCAN", "Overscan corrected")

        return overscans

    def correctGains(self, exposure, ptc, gains):
        # TODO DM 36639
        gains = []
        readNoise = []

        return gains, readNoise

    def updateVariance(self, ampExposure, amp, ptcDataset):
        """Set the variance plane using the gain and read noise.

        Parameters
        ----------
        ampExposure : `lsst.afw.image.Exposure`
            Exposure to process.
        amp : `lsst.afw.cameraGeom.Amplifier` or `FakeAmp`
            Amplifier detector data.
        ptcDataset : `lsst.ip.isr.PhotonTransferCurveDataset`
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
        gain = ptcDataset.gain[amp.getName()]

        if math.isnan(gain):
            gain = 1.0
            self.log.warning("Gain set to NAN!  Updating to 1.0 to generate Poisson variance.")
        elif gain <= 0:
            patchedGain = 1.0
            self.log.warning("Gain for amp %s == %g <= 0; setting to %f.",
                             amp.getName(), gain, patchedGain)
            gain = patchedGain

        # Get read noise from PTC
        readNoise = ptcDataset.noise[amp.getName()]

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

    def variancePlane(self, ccdExposure, ccd, ptc):
        for amp in ccd:
            if ccdExposure.getBBox().contains(amp.getBBox()):
                self.log.debug("Constructing variance map for amplifer %s.", amp.getName())
                ampExposure = ccdExposure.Factory(ccdExposure, amp.getBBox())

                self.updateVariance(ampExposure, amp, ptc)

                if self.config.qa is not None and self.config.qa.saveStats is True:
                    qaStats = afwMath.makeStatistics(ampExposure.getVariance(),
                                                     afwMath.MEDIAN | afwMath.STDEVCLIP)
                    self.log.debug("  Variance stats for amplifer %s: %f +/- %f.",
                                   amp.getName(), qaStats.getValue(afwMath.MEDIAN),
                                   qaStats.getValue(afwMath.STDEVCLIP))
        if self.config.maskNegativeVariance:
            self.maskNegativeVariance(ccdExposure)

    def maskDefects(self, exposure, defectBaseList):
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
        if self.config.doFlat and flat is not None:
            self.flatCorrection(exp, flat)
        try:
            yield exp
        finally:
            if self.config.doFlat and flat is not None:
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
            maskedImage=exposure.maskedImage,
            darkMaskedImage=darkExposure.maskedImage,
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

    def compareCameraKeywords(self, exposureMetadata, calib, calibName):
        """Compare header keywords to confirm camera states match.

        Parameters
        ----------
        exposureMetadata : `lsst.daf.base.PropertySet`
            Header for the exposure being processed.
        calib : `lsst.afw.image.Exposure` or `lsst.ip.isr.IsrCalib`
            Calibration to be applied.
        calibName : `str`
            Calib type for log message.
        """
        try:
            calibMetadata = calib.getMetadata()
        except AttributeError:
            return
        for keyword in self.config.cameraKeywordsToCompare:
            if keyword in exposureMetadata and keyword in calibMetadata:
                if exposureMetadata[keyword] != calibMetadata[keyword]:
                    if self.config.doRaiseOnCalibMismatch:
                        raise RuntimeError("Sequencer mismatch for %s [%s]: exposure: %s calib: %s",
                                           calibName, keyword,
                                           exposureMetadata[keyword], calibMetadata[keyword])
                    else:
                        self.log.warning("Sequencer mismatch for %s [%s]: exposure: %s calib: %s",
                                         calibName, keyword,
                                         exposureMetadata[keyword], calibMetadata[keyword])
            else:
                self.log.debug("Sequencer keyword %s not found.", keyword)

    def convertIntToFloat(self, exposure):
        """Convert exposure image from uint16 to float.

        If the exposure does not need to be converted, the input is
        immediately returned.  For exposures that are converted to use
        floating point pixels, the variance is set to unity and the
        mask to zero.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
           The raw exposure to be converted.

        Returns
        -------
        newexposure : `lsst.afw.image.Exposure`
           The input ``exposure``, converted to floating point pixels.

        Raises
        ------
        RuntimeError
            Raised if the exposure type cannot be converted to float.

        """
        if isinstance(exposure, afwImage.ExposureF):
            # Nothing to be done
            self.log.debug("Exposure already of type float.")
            return exposure
        if not hasattr(exposure, "convertF"):
            raise RuntimeError("Unable to convert exposure (%s) to float." % type(exposure))

        newexposure = exposure.convertF()
        newexposure.variance[:] = 1
        newexposure.mask[:] = 0x0

        return newexposure

    def ditherCounts(self, exposure, seed=None):
        """Dither the counts in the exposure.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
           The raw exposure to be dithered.
        """
        if self.config.integerDitherMode == "NONE":
            # Nothing to do here.
            return

        # FIXME: we need a seed that combines some exposure
        # metadata with detector.
        # Need to ask.
        seed = exposure.getInfo().getId()
        if seed is None:
            seed = 12345
        rng = numpy.random.RandomState(seed=seed)

        if self.config.integerDitherMode == "POSITIVE":
            low = 0.0
            high = 1.0
        elif self.config.integerDitherMode == "NEGATIVE":
            low = -1.0
            high = 0.0
        elif self.config.integerDitherMode == "SYMMETRIC":
            low = -0.5
            high = 0.5
        else:
            raise RuntimeError("Invalid config")

        exposure.image.array[:, :] += rng.uniform(low=low, high=high, size=exposure.image.array.shape)

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

    def flatCorrection(self, exposure, flatExposure, invert=False):
        """Apply flat correction in place.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to process.
        flatExposure : `lsst.afw.image.Exposure`
            Flat exposure of the same size as ``exposure``.
        invert : `Bool`, optional
            If True, unflatten an already flattened image.

        See Also
        --------
        lsst.ip.isr.isrFunctions.flatCorrection
        """
        isrFunctions.flatCorrection(
            maskedImage=exposure.getMaskedImage(),
            flatMaskedImage=flatExposure.getMaskedImage(),
            scalingType=self.config.flatScalingType,
            userScale=self.config.flatUserScale,
            invert=invert
        )

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

        bin1 = afwImage.makeExposure(bin1)
        bin2 = afwImage.makeExposure(bin2)

        bin1.setInfo(exposure.getInfo())
        bin2.setInfo(exposure.getInfo())

        return bin1, bin2

    def run(self, ccdExposure, *, dnlLUT=None, bias=None, deferredChargeCalib=None, linearizer=None,
            ptc=None, crosstalk=None, defects=None, bfKernel=None, bfGains=None, dark=None,
            flat=None, camera=None, **kwargs
            ):

        detector = ccdExposure.getDetector()

        overscanDetectorConfig = self.config.overscanCamera.getOverscanDetectorConfig(detector)

        # Validation step: check inputs match exposure configuration.
        exposureMetadata = ccdExposure.getMetadata()
        if ptc is not None:
            self.compareCameraKeywords(exposureMetadata, ptc, "PTC")

        if self.config.doDiffNonLinearCorrection:
            if dnlLUT is None:
                raise RuntimeError("doDiffNonLinearCorrection is True but no dnlLUT provided.")
            self.compareCameraKeywords(exposureMetadata, dnlLUT, "dnlLUT")
        if self.doLinearize(detector):
            if linearizer is None:
                raise RuntimeError("doLinearize is True but no linearizer provided.")
            self.compareCameraKeywords(exposureMetadata, linearizer, "linearizer")
        if self.config.doBias:
            if bias is None:
                raise RuntimeError("doBias is True but no bias provided.")
            self.compareCameraKeywords(exposureMetadata, bias, "bias")
        if self.config.doCrosstalk or overscanDetectorConfig.doAnyParallelOverscanCrosstalk:
            if crosstalk is None:
                raise RuntimeError("doCrosstalk is True but no crosstalk provided.")
            self.compareCameraKeywords(exposureMetadata, crosstalk, "crosstalk")
        if self.config.doDeferredCharge:
            if deferredChargeCalib is None:
                raise RuntimeError("doDeferredCharge is True but no deferredChargeCalib provided.")
            self.compareCameraKeywords(exposureMetadata, deferredChargeCalib, "CTI")
        if self.config.doDefect:
            if defects is None:
                raise RuntimeError("doDefect is True but no defects provided.")
            self.compareCameraKeywords(exposureMetadata, defects, "defects")
        if self.config.doDark:
            if dark is None:
                raise RuntimeError("doDark is True but no dark frame provided.")
            self.compareCameraKeywords(exposureMetadata, dark, "dark")
        if self.config.doBrighterFatter:
            if bfKernel is None:
                raise RuntimeError("doBrighterFatter is True not no bfKernel provided.")
            self.compareCameraKeywords(exposureMetadata, bfKernel, "brighter-fatter")
        if self.config.doFlat:
            if flat is None:
                raise RuntimeError("doFlat is True but not flat provided.")
            self.compareCameraKeywords(exposureMetadata, flat, "flat")

        # FIXME: Make sure that if linearity is done then it is matched
        # with the right PTC.

        # We keep track of units: start in ADU.
        exposureMetadata["LSST ISR UNITS"] = "adu"

        # First we convert the exposure to floating point values
        # (if necessary).
        self.log.debug("Converting exposure to floating point values.")
        ccdExposure = self.convertIntToFloat(ccdExposure)

        # Then we mark which amplifiers are completely bad from defects.
        badAmpDict = self.maskFullDefectAmplifiers(ccdExposure, detector, defects)

        # Now we go through ISR steps.

        # Differential non-linearity correction.
        # Units: ADU
        if self.config.doDiffNonLinearCorrection:
            self.diffNonLinearCorrection(ccdExposure, dnlLUT)

        # Dither the integer counts.
        # Input units: integerized ADU
        # Output units: floating-point ADU
        self.ditherCounts(ccdExposure)

        nominalPtcUsed = False
        if ptc is None:
            self.log.warning("No PTC provided; using nominal gain of %.3f.", self.config.nominalGain)
            nominalPtcUsed = True
            ptc = PhotonTransferCurveDataset([amp.getName() for amp in detector], "NOMINAL_PTC", 1)
            for amp in detector:
                ptc.gain[amp.getName()] = self.config.nominalGain
                ptc.noise[amp.getName()] = 0.0

        gains = ptc.gain

        # Serial overscan correction.
        # Units: ADU
        if overscanDetectorConfig.doAnySerialOverscan:
            serialOverscans = self.overscanCorrection(
                "SERIAL",
                overscanDetectorConfig,
                detector,
                badAmpDict,
                ccdExposure,
            )

            if nominalPtcUsed:
                # Get the empirical read noise
                # Log this; also put in metadata.
                for amp, serialOverscan in zip(detector, serialOverscans):
                    if serialOverscan is None:
                        ptc.noise[amp.getName()] = 0.0
                    else:
                        ptc.noise[amp.getName()] = serialOverscan.residualSigma
        else:
            serialOverscans = [None]*len(detector)

        # After serial overscan correction, we can mask SATURATED and
        # SUSPECT pixels. This updates badAmpDict if any amplifier
        # is fully saturated after serial overscan correction.

        # The saturation is currently assumed to be recorded in
        # overscan-corrected ADU.
        badAmpDict = self.maskSaturatedPixels(badAmpDict, ccdExposure, detector)

        # Parallel overscan crosstalk correction.
        # Units: ADU
        if overscanDetectorConfig.doAnyParallelOverscanCrosstalk:
            # We supply the gains to the crosstalk task to allow for adu to be
            # corrected with an electron matrix.
            self.crosstalk.run(
                ccdExposure,
                crosstalk=crosstalk,
                camera=camera,
                parallelOverscanRegion=True,
                detectorConfig=overscanDetectorConfig,
                doSqrCrosstalk=False,
                gains=gains,
            )

        # Parallel overscan correction.
        # Units: ADU
        if overscanDetectorConfig.doAnyParallelOverscan:
            # At the moment we do not use the return values from this task.
            _ = self.overscanCorrection(
                "PARALLEL",
                overscanDetectorConfig,
                detector,
                badAmpDict,
                ccdExposure,
            )

        if self.config.doCorrectGains:
            # TODO DM 36639
            # This requires the PTC (?) with the temperature dependence.
            self.log.info("Apply temperature dependence to the gains.")
            gains, readNoise = self.correctGains(ccdExposure, ptc, gains)

        # Do gain normalization; this may be the nominal gains.
        if self.config.doApplyGains:
            # Input units: ADU
            # Output units: electrons
            self.log.info("Using gain values to convert from ADU to electrons.")
            isrFunctions.applyGains(ccdExposure, normalizeGains=False, ptcGains=gains, isTrimmed=False)
            # The units are now electrons.
            exposureMetadata["LSST ISR UNITS"] = "electron"

        # Do crosstalk correction in the imaging region.
        # Units: electrons
        if self.config.doCrosstalk:
            self.log.info("Applying crosstalk corrections.")
            self.crosstalk.run(ccdExposure, crosstalk=crosstalk, isTrimmed=False, gains=gains)

        # Linearity correction
        # FIXME: watch those units here; linearity code may need update.
        # Units: electrons
        if self.config.doLinearize:
            self.log.info("Applying linearizer.")
            linearizer.applyLinearity(image=ccdExposure.image, detector=detector, log=self.log, gains=gains)

        # Serial CTI (deferred charge)
        # FIXME: watch out for gain units; cti code may need update
        # (to make it simpler!)
        # Units: electrons
        if self.config.doDeferredCharge:
            self.log.info("Applying deferred charge/CTI correction.")
            self.deferredChargeCorrection.run(ccdExposure, deferredChargeCalib)

        # Assemble/trim
        # Units: electrons
        if self.config.doAssembleCcd:
            self.log.info("Assembling CCD from amplifiers.")
            ccdExposure = self.assembleCcd.assembleCcd(ccdExposure)

            if self.config.expectWcs and not ccdExposure.getWcs():
                self.log.warning("No WCS found in input exposure.")

        # Bias subtraction
        # Units: electrons
        if self.config.doBias:
            self.log.info("Applying bias correction.")
            isrFunctions.biasCorrection(ccdExposure.maskedImage, bias.maskedImage)

        # Dark subtraction
        # Units: electrons
        if self.config.doDark:
            self.log.info("Applying dark subtraction.")
            self.darkCorrection(ccdExposure, dark)

        # Defect masking
        # Masking block (defects, NAN pixels and trails).
        # Saturated and suspect pixels have already been masked.
        # Units: electrons
        if self.config.doDefect:
            self.log.info("Applying defect masking.")
            self.maskDefects(ccdExposure, defects)

        if self.config.doNanMasking:
            self.log.info("Masking non-finite (NAN, inf) value pixels.")
            self.maskNan(ccdExposure)

        if self.config.doWidenSaturationTrails:
            self.log.info("Widening saturation trails.")
            isrFunctions.widenSaturationTrails(ccdExposure.getMaskedImage().getMask())

        # Brighter/Fatter
        # FIXME: watch out (again) for gain units.
        # Units: electrons
        if self.config.doBrighterFatter:
            self.log.info("Applying Bright-Fatter kernels.")
            bfKernelOut, bfGains = self.getBrighterFatterKernel(detector, bfKernel)
            ccdExposure = self.applyBrighterFatterCorrection(ccdExposure, flat, dark, bfKernelOut, bfGains)

        # Variance plane creation
        # Units: electrons
        if self.config.doVariance:
            self.variancePlane(ccdExposure, detector, ptc)

        # Flat-fielding
        # This may move elsewhere.
        # Placeholder while the LSST flat procedure is done.
        # Units: electrons
        if self.config.doFlat:
            self.log.info("Applying flat correction.")
            self.flatCorrection(ccdExposure, flat)

        # Pixel values for masked regions are set here
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

        # Calculate amp offset corrections within the CCD.
        if self.config.doAmpOffset:
            if self.config.ampOffset.doApplyAmpOffset:
                self.log.info("Calculating and applying amp offset corrections.")
            else:
                self.log.info("Calculating amp offset corrections without applying them.")
            self.ampOffset.run(ccdExposure)

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
                if overscanDetectorConfig.doAnySerialOverscan and k1 in metadata and k2 in metadata:
                    metadata[f"LSST ISR LEVEL {ampName}"] = metadata[k1] - metadata[k2]
                else:
                    metadata[f"LSST ISR LEVEL {ampName}"] = numpy.nan

        # calculate additional statistics.
        outputStatistics = None
        if self.config.doCalculateStatistics:
            outputStatistics = self.isrStats.run(ccdExposure, overscanResults=serialOverscans,
                                                 bias=bias, dark=dark, flat=flat, ptc=ptc,
                                                 defects=defects).results

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
