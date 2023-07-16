# This file is part of ip_isr.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

__all__ = ["IsrTask", "IsrTaskConfig"]

import math
import numpy

import lsst.geom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as cT

from contextlib import contextmanager
from lsstDebug import getDebugFrame

from lsst.afw.cameraGeom import NullLinearityType
from lsst.afw.display import getDisplay
from lsst.meas.algorithms.detection import SourceDetectionTask
from lsst.utils.timer import timeMethod

from . import isrFunctions
from . import isrQa
from . import linearize
from .defects import Defects

from .assembleCcdTask import AssembleCcdTask
from .crosstalk import CrosstalkTask, CrosstalkCalib
from .fringe import FringeTask
from .isr import maskNans
from .masking import MaskingTask
from .overscan import OverscanCorrectionTask
from .straylight import StrayLightTask
from .vignette import VignetteTask
from .ampOffset import AmpOffsetTask
from .deferredCharge import DeferredChargeTask
from .isrStatistics import IsrStatisticsTask
from lsst.daf.butler import DimensionGraph


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


class IsrTaskConnections(pipeBase.PipelineTaskConnections,
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
    bias = cT.PrerequisiteInput(
        name="bias",
        doc="Input bias calibration.",
        storageClass="ExposureF",
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
    flat = cT.PrerequisiteInput(
        name="flat",
        doc="Input flat calibration.",
        storageClass="ExposureF",
        dimensions=["instrument", "physical_filter", "detector"],
        isCalibration=True,
    )
    ptc = cT.PrerequisiteInput(
        name="ptc",
        doc="Input Photon Transfer Curve dataset",
        storageClass="PhotonTransferCurveDataset",
        dimensions=["instrument", "detector"],
        isCalibration=True,
    )
    fringes = cT.PrerequisiteInput(
        name="fringe",
        doc="Input fringe calibration.",
        storageClass="ExposureF",
        dimensions=["instrument", "physical_filter", "detector"],
        isCalibration=True,
        minimum=0,  # only needed for some bands, even when enabled
    )
    strayLightData = cT.PrerequisiteInput(
        name='yBackground',
        doc="Input stray light calibration.",
        storageClass="StrayLightData",
        dimensions=["instrument", "physical_filter", "detector"],
        deferLoad=True,
        isCalibration=True,
        minimum=0,  # only needed for some bands, even when enabled
    )
    bfKernel = cT.PrerequisiteInput(
        name='bfKernel',
        doc="Input brighter-fatter kernel.",
        storageClass="NumpyArray",
        dimensions=["instrument"],
        isCalibration=True,
        minimum=0,  # can use either bfKernel or newBFKernel
    )
    newBFKernel = cT.PrerequisiteInput(
        name='brighterFatterKernel',
        doc="Newer complete kernel + gain solutions.",
        storageClass="BrighterFatterKernel",
        dimensions=["instrument", "detector"],
        isCalibration=True,
        minimum=0,  # can use either bfKernel or newBFKernel
    )
    defects = cT.PrerequisiteInput(
        name='defects',
        doc="Input defect tables.",
        storageClass="Defects",
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
    opticsTransmission = cT.PrerequisiteInput(
        name="transmission_optics",
        storageClass="TransmissionCurve",
        doc="Transmission curve due to the optics.",
        dimensions=["instrument"],
        isCalibration=True,
    )
    filterTransmission = cT.PrerequisiteInput(
        name="transmission_filter",
        storageClass="TransmissionCurve",
        doc="Transmission curve due to the filter.",
        dimensions=["instrument", "physical_filter"],
        isCalibration=True,
    )
    sensorTransmission = cT.PrerequisiteInput(
        name="transmission_sensor",
        storageClass="TransmissionCurve",
        doc="Transmission curve due to the sensor.",
        dimensions=["instrument", "detector"],
        isCalibration=True,
    )
    atmosphereTransmission = cT.PrerequisiteInput(
        name="transmission_atmosphere",
        storageClass="TransmissionCurve",
        doc="Transmission curve due to the atmosphere.",
        dimensions=["instrument"],
        isCalibration=True,
    )
    illumMaskedImage = cT.PrerequisiteInput(
        name="illum",
        doc="Input illumination correction.",
        storageClass="MaskedImageF",
        dimensions=["instrument", "physical_filter", "detector"],
        isCalibration=True,
    )
    deferredChargeCalib = cT.PrerequisiteInput(
        name="cpCtiCalib",
        doc="Deferred charge/CTI correction dataset.",
        storageClass="IsrCalib",
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
    outputOssThumbnail = cT.Output(
        name="OssThumb",
        doc="Output Overscan-subtracted thumbnail image.",
        storageClass="Thumbnail",
        dimensions=["instrument", "exposure", "detector"],
    )
    outputFlattenedThumbnail = cT.Output(
        name="FlattenedThumb",
        doc="Output flat-corrected thumbnail image.",
        storageClass="Thumbnail",
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

        if config.doBias is not True:
            self.prerequisiteInputs.remove("bias")
        if config.doLinearize is not True:
            self.prerequisiteInputs.remove("linearizer")
        if config.doCrosstalk is not True:
            self.prerequisiteInputs.remove("crosstalkSources")
            self.prerequisiteInputs.remove("crosstalk")
        if config.doBrighterFatter is not True:
            self.prerequisiteInputs.remove("bfKernel")
            self.prerequisiteInputs.remove("newBFKernel")
        if config.doDefect is not True:
            self.prerequisiteInputs.remove("defects")
        if config.doDark is not True:
            self.prerequisiteInputs.remove("dark")
        if config.doFlat is not True:
            self.prerequisiteInputs.remove("flat")
        if config.doFringe is not True:
            self.prerequisiteInputs.remove("fringes")
        if config.doStrayLight is not True:
            self.prerequisiteInputs.remove("strayLightData")
        if config.usePtcGains is not True and config.usePtcReadNoise is not True:
            self.prerequisiteInputs.remove("ptc")
        if config.doAttachTransmissionCurve is not True:
            self.prerequisiteInputs.remove("opticsTransmission")
            self.prerequisiteInputs.remove("filterTransmission")
            self.prerequisiteInputs.remove("sensorTransmission")
            self.prerequisiteInputs.remove("atmosphereTransmission")
        else:
            if config.doUseOpticsTransmission is not True:
                self.prerequisiteInputs.remove("opticsTransmission")
            if config.doUseFilterTransmission is not True:
                self.prerequisiteInputs.remove("filterTransmission")
            if config.doUseSensorTransmission is not True:
                self.prerequisiteInputs.remove("sensorTransmission")
            if config.doUseAtmosphereTransmission is not True:
                self.prerequisiteInputs.remove("atmosphereTransmission")
        if config.doIlluminationCorrection is not True:
            self.prerequisiteInputs.remove("illumMaskedImage")
        if config.doDeferredCharge is not True:
            self.prerequisiteInputs.remove("deferredChargeCalib")

        if config.doWrite is not True:
            self.outputs.remove("outputExposure")
            self.outputs.remove("preInterpExposure")
            self.outputs.remove("outputFlattenedThumbnail")
            self.outputs.remove("outputOssThumbnail")
            self.outputs.remove("outputStatistics")

        if config.doSaveInterpPixels is not True:
            self.outputs.remove("preInterpExposure")
        if config.qa.doThumbnailOss is not True:
            self.outputs.remove("outputOssThumbnail")
        if config.qa.doThumbnailFlattened is not True:
            self.outputs.remove("outputFlattenedThumbnail")
        if config.doCalculateStatistics is not True:
            self.outputs.remove("outputStatistics")


class IsrTaskConfig(pipeBase.PipelineTaskConfig,
                    pipelineConnections=IsrTaskConnections):
    """Configuration parameters for IsrTask.

    Items are grouped in the order in which they are executed by the task.
    """
    datasetType = pexConfig.Field(
        dtype=str,
        doc="Dataset type for input data; users will typically leave this alone, "
        "but camera-specific ISR tasks will override it",
        default="raw",
    )

    fallbackFilterName = pexConfig.Field(
        dtype=str,
        doc="Fallback default filter name for calibrations.",
        optional=True
    )
    useFallbackDate = pexConfig.Field(
        dtype=bool,
        doc="Pass observation date when using fallback filter.",
        default=False,
    )
    expectWcs = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Expect input science images to have a WCS (set False for e.g. spectrographs)."
    )
    fwhm = pexConfig.Field(
        dtype=float,
        doc="FWHM of PSF in arcseconds (currently unused).",
        default=1.0,
    )
    qa = pexConfig.ConfigField(
        dtype=isrQa.IsrQaConfig,
        doc="QA related configuration options.",
    )
    doHeaderProvenance = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Write calibration identifiers into output exposure header?",
    )

    # Calib checking configuration:
    doRaiseOnCalibMismatch = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Should IsrTask halt if exposure and calibration header values do not match?",
    )
    cameraKeywordsToCompare = pexConfig.ListField(
        dtype=str,
        doc="List of header keywords to compare between exposure and calibrations.",
        default=[],
    )

    # Image conversion configuration
    doConvertIntToFloat = pexConfig.Field(
        dtype=bool,
        doc="Convert integer raw images to floating point values?",
        default=True,
    )

    # Saturated pixel handling.
    doSaturation = pexConfig.Field(
        dtype=bool,
        doc="Mask saturated pixels? NB: this is totally independent of the"
        " interpolation option - this is ONLY setting the bits in the mask."
        " To have them interpolated make sure doSaturationInterpolation=True",
        default=True,
    )
    saturatedMaskName = pexConfig.Field(
        dtype=str,
        doc="Name of mask plane to use in saturation detection and interpolation",
        default="SAT",
    )
    saturation = pexConfig.Field(
        dtype=float,
        doc="The saturation level to use if no Detector is present in the Exposure (ignored if NaN)",
        default=float("NaN"),
    )
    growSaturationFootprintSize = pexConfig.Field(
        dtype=int,
        doc="Number of pixels by which to grow the saturation footprints",
        default=1,
    )

    # Suspect pixel handling.
    doSuspect = pexConfig.Field(
        dtype=bool,
        doc="Mask suspect pixels?",
        default=False,
    )
    suspectMaskName = pexConfig.Field(
        dtype=str,
        doc="Name of mask plane to use for suspect pixels",
        default="SUSPECT",
    )
    numEdgeSuspect = pexConfig.Field(
        dtype=int,
        doc="Number of edge pixels to be flagged as untrustworthy.",
        default=0,
    )
    edgeMaskLevel = pexConfig.ChoiceField(
        dtype=str,
        doc="Mask edge pixels in which coordinate frame: DETECTOR or AMP?",
        default="DETECTOR",
        allowed={
            'DETECTOR': 'Mask only the edges of the full detector.',
            'AMP': 'Mask edges of each amplifier.',
        },
    )

    # Initial masking options.
    doSetBadRegions = pexConfig.Field(
        dtype=bool,
        doc="Should we set the level of all BAD patches of the chip to the chip's average value?",
        default=True,
    )
    badStatistic = pexConfig.ChoiceField(
        dtype=str,
        doc="How to estimate the average value for BAD regions.",
        default='MEANCLIP',
        allowed={
            "MEANCLIP": "Correct using the (clipped) mean of good data",
            "MEDIAN": "Correct using the median of the good data",
        },
    )

    # Overscan subtraction configuration.
    doOverscan = pexConfig.Field(
        dtype=bool,
        doc="Do overscan subtraction?",
        default=True,
    )
    overscan = pexConfig.ConfigurableField(
        target=OverscanCorrectionTask,
        doc="Overscan subtraction task for image segments.",
    )

    # Amplifier to CCD assembly configuration
    doAssembleCcd = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Assemble amp-level exposures into a ccd-level exposure?"
    )
    assembleCcd = pexConfig.ConfigurableField(
        target=AssembleCcdTask,
        doc="CCD assembly task",
    )

    # General calibration configuration.
    doAssembleIsrExposures = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Assemble amp-level calibration exposures into ccd-level exposure?"
    )
    doTrimToMatchCalib = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Trim raw data to match calibration bounding boxes?"
    )

    # Bias subtraction.
    doBias = pexConfig.Field(
        dtype=bool,
        doc="Apply bias frame correction?",
        default=True,
    )
    biasDataProductName = pexConfig.Field(
        dtype=str,
        doc="Name of the bias data product",
        default="bias",
    )
    doBiasBeforeOverscan = pexConfig.Field(
        dtype=bool,
        doc="Reverse order of overscan and bias correction.",
        default=False
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

    # Variance construction
    doVariance = pexConfig.Field(
        dtype=bool,
        doc="Calculate variance?",
        default=True
    )
    gain = pexConfig.Field(
        dtype=float,
        doc="The gain to use if no Detector is present in the Exposure (ignored if NaN)",
        default=float("NaN"),
    )
    readNoise = pexConfig.Field(
        dtype=float,
        doc="The read noise to use if no Detector is present in the Exposure",
        default=0.0,
    )
    doEmpiricalReadNoise = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Calculate empirical read noise instead of value from AmpInfo data?"
    )
    usePtcReadNoise = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Use readnoise values from the Photon Transfer Curve?"
    )
    maskNegativeVariance = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Mask pixels that claim a negative variance?  This likely indicates a failure "
        "in the measurement of the overscan at an edge due to the data falling off faster "
        "than the overscan model can account for it."
    )
    negativeVarianceMaskName = pexConfig.Field(
        dtype=str,
        default="BAD",
        doc="Mask plane to use to mark pixels with negative variance, if `maskNegativeVariance` is True.",
    )
    # Linearization.
    doLinearize = pexConfig.Field(
        dtype=bool,
        doc="Correct for nonlinearity of the detector's response?",
        default=True,
    )

    # Crosstalk.
    doCrosstalk = pexConfig.Field(
        dtype=bool,
        doc="Apply intra-CCD crosstalk correction?",
        default=False,
    )
    doCrosstalkBeforeAssemble = pexConfig.Field(
        dtype=bool,
        doc="Apply crosstalk correction before CCD assembly, and before trimming?",
        default=False,
    )
    crosstalk = pexConfig.ConfigurableField(
        target=CrosstalkTask,
        doc="Intra-CCD crosstalk correction",
    )

    # Masking options.
    doDefect = pexConfig.Field(
        dtype=bool,
        doc="Apply correction for CCD defects, e.g. hot pixels?",
        default=True,
    )
    doNanMasking = pexConfig.Field(
        dtype=bool,
        doc="Mask non-finite (NAN, inf) pixels?",
        default=True,
    )
    doWidenSaturationTrails = pexConfig.Field(
        dtype=bool,
        doc="Widen bleed trails based on their width?",
        default=True
    )

    # Brighter-Fatter correction.
    doBrighterFatter = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Apply the brighter-fatter correction?"
    )
    doFluxConservingBrighterFatterCorrection = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Apply the flux-conserving BFE correction by Miller et al.?"
    )
    brighterFatterLevel = pexConfig.ChoiceField(
        dtype=str,
        default="DETECTOR",
        doc="The level at which to correct for brighter-fatter.",
        allowed={
            "AMP": "Every amplifier treated separately.",
            "DETECTOR": "One kernel per detector",
        }
    )
    brighterFatterMaxIter = pexConfig.Field(
        dtype=int,
        default=10,
        doc="Maximum number of iterations for the brighter-fatter correction"
    )
    brighterFatterThreshold = pexConfig.Field(
        dtype=float,
        default=1000,
        doc="Threshold used to stop iterating the brighter-fatter correction.  It is the "
        "absolute value of the difference between the current corrected image and the one "
        "from the previous iteration summed over all the pixels."
    )
    brighterFatterApplyGain = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Should the gain be applied when applying the brighter-fatter correction?"
    )
    brighterFatterMaskListToInterpolate = pexConfig.ListField(
        dtype=str,
        doc="List of mask planes that should be interpolated over when applying the brighter-fatter "
        "correction.",
        default=["SAT", "BAD", "NO_DATA", "UNMASKEDNAN"],
    )
    brighterFatterMaskGrowSize = pexConfig.Field(
        dtype=int,
        default=0,
        doc="Number of pixels to grow the masks listed in config.brighterFatterMaskListToInterpolate "
        "when brighter-fatter correction is applied."
    )

    # Dark subtraction.
    doDark = pexConfig.Field(
        dtype=bool,
        doc="Apply dark frame correction?",
        default=True,
    )
    darkDataProductName = pexConfig.Field(
        dtype=str,
        doc="Name of the dark data product",
        default="dark",
    )

    # Camera-specific stray light removal.
    doStrayLight = pexConfig.Field(
        dtype=bool,
        doc="Subtract stray light in the y-band (due to encoder LEDs)?",
        default=False,
    )
    strayLight = pexConfig.ConfigurableField(
        target=StrayLightTask,
        doc="y-band stray light correction"
    )

    # Flat correction.
    doFlat = pexConfig.Field(
        dtype=bool,
        doc="Apply flat field correction?",
        default=True,
    )
    flatDataProductName = pexConfig.Field(
        dtype=str,
        doc="Name of the flat data product",
        default="flat",
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
        doc="If flatScalingType is 'USER' then scale flat by this amount; ignored otherwise",
        default=1.0,
    )
    doTweakFlat = pexConfig.Field(
        dtype=bool,
        doc="Tweak flats to match observed amplifier ratios?",
        default=False
    )

    # Amplifier normalization based on gains instead of using flats
    # configuration.
    doApplyGains = pexConfig.Field(
        dtype=bool,
        doc="Correct the amplifiers for their gains instead of applying flat correction",
        default=False,
    )
    usePtcGains = pexConfig.Field(
        dtype=bool,
        doc="Use the gain values from the Photon Transfer Curve?",
        default=False,
    )
    normalizeGains = pexConfig.Field(
        dtype=bool,
        doc="Normalize all the amplifiers in each CCD to have the same median value.",
        default=False,
    )

    # Fringe correction.
    doFringe = pexConfig.Field(
        dtype=bool,
        doc="Apply fringe correction?",
        default=True,
    )
    fringe = pexConfig.ConfigurableField(
        target=FringeTask,
        doc="Fringe subtraction task",
    )
    fringeAfterFlat = pexConfig.Field(
        dtype=bool,
        doc="Do fringe subtraction after flat-fielding?",
        default=True,
    )

    # Amp offset correction.
    doAmpOffset = pexConfig.Field(
        doc="Calculate and apply amp offset corrections?",
        dtype=bool,
        default=False,
    )
    ampOffset = pexConfig.ConfigurableField(
        doc="Amp offset correction task.",
        target=AmpOffsetTask,
    )

    # Initial CCD-level background statistics options.
    doMeasureBackground = pexConfig.Field(
        dtype=bool,
        doc="Measure the background level on the reduced image?",
        default=False,
    )

    # Camera-specific masking configuration.
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
    doSaturationInterpolation = pexConfig.Field(
        dtype=bool,
        doc="Perform interpolation over pixels masked as saturated?"
        " NB: This is independent of doSaturation; if that is False this plane"
        " will likely be blank, resulting in a no-op here.",
        default=True,
    )
    doNanInterpolation = pexConfig.Field(
        dtype=bool,
        doc="Perform interpolation over pixels masked as NaN?"
        " NB: This is independent of doNanMasking; if that is False this plane"
        " will likely be blank, resulting in a no-op here.",
        default=True,
    )
    doNanInterpAfterFlat = pexConfig.Field(
        dtype=bool,
        doc=("If True, ensure we interpolate NaNs after flat-fielding, even if we "
             "also have to interpolate them before flat-fielding."),
        default=False,
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

    # Default photometric calibration options.
    fluxMag0T1 = pexConfig.DictField(
        keytype=str,
        itemtype=float,
        doc="The approximate flux of a zero-magnitude object in a one-second exposure, per filter.",
        default=dict((f, pow(10.0, 0.4*m)) for f, m in (("Unknown", 28.0),
                                                        ))
    )
    defaultFluxMag0T1 = pexConfig.Field(
        dtype=float,
        doc="Default value for fluxMag0T1 (for an unrecognized filter).",
        default=pow(10.0, 0.4*28.0)
    )

    # Vignette correction configuration.
    doVignette = pexConfig.Field(
        dtype=bool,
        doc=("Compute and attach the validPolygon defining the unvignetted region to the exposure "
             "according to vignetting parameters?"),
        default=False,
    )
    doMaskVignettePolygon = pexConfig.Field(
        dtype=bool,
        doc=("Add a mask bit for pixels within the vignetted region.  Ignored if doVignette "
             "is False"),
        default=True,
    )
    vignetteValue = pexConfig.Field(
        dtype=float,
        doc="Value to replace image array pixels with in the vignetted region?  Ignored if None.",
        optional=True,
        default=None,
    )
    vignette = pexConfig.ConfigurableField(
        target=VignetteTask,
        doc="Vignetting task.",
    )

    # Transmission curve configuration.
    doAttachTransmissionCurve = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Construct and attach a wavelength-dependent throughput curve for this CCD image?"
    )
    doUseOpticsTransmission = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Load and use transmission_optics (if doAttachTransmissionCurve is True)?"
    )
    doUseFilterTransmission = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Load and use transmission_filter (if doAttachTransmissionCurve is True)?"
    )
    doUseSensorTransmission = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Load and use transmission_sensor (if doAttachTransmissionCurve is True)?"
    )
    doUseAtmosphereTransmission = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Load and use transmission_atmosphere (if doAttachTransmissionCurve is True)?"
    )

    # Illumination correction.
    doIlluminationCorrection = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Perform illumination correction?"
    )
    illuminationCorrectionDataProductName = pexConfig.Field(
        dtype=str,
        doc="Name of the illumination correction data product.",
        default="illumcor",
    )
    illumScale = pexConfig.Field(
        dtype=float,
        doc="Scale factor for the illumination correction.",
        default=1.0,
    )
    illumFilters = pexConfig.ListField(
        dtype=str,
        default=[],
        doc="Only perform illumination correction for these filters."
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

    # Write the outputs to disk. If ISR is run as a subtask, this may not
    # be needed.
    doWrite = pexConfig.Field(
        dtype=bool,
        doc="Persist postISRCCD?",
        default=True,
    )

    def validate(self):
        super().validate()
        if self.doFlat and self.doApplyGains:
            raise ValueError("You may not specify both doFlat and doApplyGains")
        if self.doBiasBeforeOverscan and self.doTrimToMatchCalib:
            raise ValueError("You may not specify both doBiasBeforeOverscan and doTrimToMatchCalib")
        if self.doSaturationInterpolation and self.saturatedMaskName not in self.maskListToInterpolate:
            self.maskListToInterpolate.append(self.saturatedMaskName)
        if not self.doSaturationInterpolation and self.saturatedMaskName in self.maskListToInterpolate:
            self.maskListToInterpolate.remove(self.saturatedMaskName)
        if self.doNanInterpolation and "UNMASKEDNAN" not in self.maskListToInterpolate:
            self.maskListToInterpolate.append("UNMASKEDNAN")


class IsrTask(pipeBase.PipelineTask):
    """Apply common instrument signature correction algorithms to a raw frame.

    The process for correcting imaging data is very similar from
    camera to camera.  This task provides a vanilla implementation of
    doing these corrections, including the ability to turn certain
    corrections off if they are not needed.  The inputs to the primary
    method, `run()`, are a raw exposure to be corrected and the
    calibration data products. The raw input is a single chip sized
    mosaic of all amps including overscans and other non-science
    pixels.

    The __init__ method sets up the subtasks for ISR processing, using
    the defaults from `lsst.ip.isr`.

    Parameters
    ----------
    args : `list`
        Positional arguments passed to the Task constructor.
        None used at this time.
    kwargs : `dict`, optional
        Keyword arguments passed on to the Task constructor.
        None used at this time.
    """
    ConfigClass = IsrTaskConfig
    _DefaultName = "isr"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.makeSubtask("assembleCcd")
        self.makeSubtask("crosstalk")
        self.makeSubtask("strayLight")
        self.makeSubtask("fringe")
        self.makeSubtask("masking")
        self.makeSubtask("overscan")
        self.makeSubtask("vignette")
        self.makeSubtask("ampOffset")
        self.makeSubtask("deferredChargeCorrection")
        self.makeSubtask("isrStats")

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        try:
            inputs['detectorNum'] = inputRefs.ccdExposure.dataId['detector']
        except Exception as e:
            raise ValueError("Failure to find valid detectorNum value for Dataset %s: %s." %
                             (inputRefs, e))

        detector = inputs['ccdExposure'].getDetector()

        if self.config.doCrosstalk is True:
            # Crosstalk sources need to be defined by the pipeline
            # yaml if they exist.
            if 'crosstalk' in inputs and inputs['crosstalk'] is not None:
                if not isinstance(inputs['crosstalk'], CrosstalkCalib):
                    inputs['crosstalk'] = CrosstalkCalib.fromTable(inputs['crosstalk'])
            else:
                coeffVector = (self.config.crosstalk.crosstalkValues
                               if self.config.crosstalk.useConfigCoefficients else None)
                crosstalkCalib = CrosstalkCalib().fromDetector(detector, coeffVector=coeffVector)
                inputs['crosstalk'] = crosstalkCalib
            if inputs['crosstalk'].interChip and len(inputs['crosstalk'].interChip) > 0:
                if 'crosstalkSources' not in inputs:
                    self.log.warning("No crosstalkSources found for chip with interChip terms!")

        if self.doLinearize(detector) is True:
            if 'linearizer' in inputs:
                if isinstance(inputs['linearizer'], dict):
                    linearizer = linearize.Linearizer(detector=detector, log=self.log)
                    linearizer.fromYaml(inputs['linearizer'])
                    self.log.warning("Dictionary linearizers will be deprecated in DM-28741.")
                elif isinstance(inputs['linearizer'], numpy.ndarray):
                    linearizer = linearize.Linearizer(table=inputs.get('linearizer', None),
                                                      detector=detector,
                                                      log=self.log)
                    self.log.warning("Bare lookup table linearizers will be deprecated in DM-28741.")
                else:
                    linearizer = inputs['linearizer']
                    linearizer.log = self.log
                inputs['linearizer'] = linearizer
            else:
                inputs['linearizer'] = linearize.Linearizer(detector=detector, log=self.log)
                self.log.warning("Constructing linearizer from cameraGeom information.")

        if self.config.doDefect is True:
            if "defects" in inputs and inputs['defects'] is not None:
                # defects is loaded as a BaseCatalog with columns
                # x0, y0, width, height. Masking expects a list of defects
                # defined by their bounding box
                if not isinstance(inputs["defects"], Defects):
                    inputs["defects"] = Defects.fromTable(inputs["defects"])

        # Load the correct style of brighter-fatter kernel, and repack
        # the information as a numpy array.
        if self.config.doBrighterFatter:
            brighterFatterKernel = inputs.pop('newBFKernel', None)
            if brighterFatterKernel is None:
                # This type of kernel must be in (y, x) index
                # ordering, as it used directly as the .array
                # component of the afwImage kernel.
                brighterFatterKernel = inputs.get('bfKernel', None)

            if brighterFatterKernel is not None and not isinstance(brighterFatterKernel, numpy.ndarray):
                # This is a ISR calib kernel.  These kernels are
                # generated in (x, y) index ordering, and need to be
                # transposed to be used directly as the .array
                # component of the afwImage kernel.  This is done
                # explicitly below when setting the ``bfKernel``
                # input.
                detName = detector.getName()
                level = brighterFatterKernel.level

                # This is expected to be a dictionary of amp-wise gains.
                inputs['bfGains'] = brighterFatterKernel.gain
                if self.config.brighterFatterLevel == 'DETECTOR':
                    kernel = None
                    if level == 'DETECTOR':
                        if detName in brighterFatterKernel.detKernels:
                            kernel = brighterFatterKernel.detKernels[detName]
                        else:
                            raise RuntimeError("Failed to extract kernel from new-style BF kernel.")
                    elif level == 'AMP':
                        self.log.warning("Making DETECTOR level kernel from AMP based brighter "
                                         "fatter kernels.")
                        brighterFatterKernel.makeDetectorKernelFromAmpwiseKernels(detName)
                        kernel = brighterFatterKernel.detKernels[detName]
                    if kernel is None:
                        raise RuntimeError("Could not identify brighter-fatter kernel!")
                    # Do the one single transpose here so the kernel
                    # can be directly loaded into the afwImage .array
                    # component.
                    inputs['bfKernel'] = numpy.transpose(kernel)
                elif self.config.brighterFatterLevel == 'AMP':
                    raise NotImplementedError("Per-amplifier brighter-fatter correction not implemented")

        if self.config.doFringe is True and self.fringe.checkFilter(inputs['ccdExposure']):
            expId = inputs['ccdExposure'].info.id
            inputs['fringes'] = self.fringe.loadFringes(inputs['fringes'],
                                                        expId=expId,
                                                        assembler=self.assembleCcd
                                                        if self.config.doAssembleIsrExposures else None)
        else:
            inputs['fringes'] = pipeBase.Struct(fringes=None)

        if self.config.doStrayLight is True and self.strayLight.checkFilter(inputs['ccdExposure']):
            if 'strayLightData' not in inputs:
                inputs['strayLightData'] = None

        if self.config.doHeaderProvenance:
            # Add calibration provenanace info to header.
            exposureMetadata = inputs['ccdExposure'].getMetadata()
            for inputName in sorted(inputs.keys()):
                reference = getattr(inputRefs, inputName, None)
                if reference is not None and hasattr(reference, "run"):
                    runKey = f"LSST CALIB RUN {inputName.upper()}"
                    runValue = reference.run
                    idKey = f"LSST CALIB UUID {inputName.upper()}"
                    idValue = str(reference.id)

                    exposureMetadata[runKey] = runValue
                    exposureMetadata[idKey] = idValue

        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)

    @timeMethod
    def run(self, ccdExposure, *, camera=None, bias=None, linearizer=None,
            crosstalk=None, crosstalkSources=None,
            dark=None, flat=None, ptc=None, bfKernel=None, bfGains=None, defects=None,
            fringes=pipeBase.Struct(fringes=None), opticsTransmission=None, filterTransmission=None,
            sensorTransmission=None, atmosphereTransmission=None,
            detectorNum=None, strayLightData=None, illumMaskedImage=None,
            deferredChargeCalib=None,
            ):
        """Perform instrument signature removal on an exposure.

        Steps included in the ISR processing, in order performed, are:

        - saturation and suspect pixel masking
        - overscan subtraction
        - CCD assembly of individual amplifiers
        - bias subtraction
        - variance image construction
        - linearization of non-linear response
        - crosstalk masking
        - brighter-fatter correction
        - dark subtraction
        - fringe correction
        - stray light subtraction
        - flat correction
        - masking of known defects and camera specific features
        - vignette calculation
        - appending transmission curve and distortion model

        Parameters
        ----------
        ccdExposure : `lsst.afw.image.Exposure`
            The raw exposure that is to be run through ISR.  The
            exposure is modified by this method.
        camera : `lsst.afw.cameraGeom.Camera`, optional
            The camera geometry for this exposure. Required if
            one or more of ``ccdExposure``, ``bias``, ``dark``, or
            ``flat`` does not have an associated detector.
        bias : `lsst.afw.image.Exposure`, optional
            Bias calibration frame.
        linearizer : `lsst.ip.isr.linearize.LinearizeBase`, optional
            Functor for linearization.
        crosstalk : `lsst.ip.isr.crosstalk.CrosstalkCalib`, optional
            Calibration for crosstalk.
        crosstalkSources : `list`, optional
            List of possible crosstalk sources.
        dark : `lsst.afw.image.Exposure`, optional
            Dark calibration frame.
        flat : `lsst.afw.image.Exposure`, optional
            Flat calibration frame.
        ptc : `lsst.ip.isr.PhotonTransferCurveDataset`, optional
            Photon transfer curve dataset, with, e.g., gains
            and read noise.
        bfKernel : `numpy.ndarray`, optional
            Brighter-fatter kernel.
        bfGains : `dict` of `float`, optional
            Gains used to override the detector's nominal gains for the
            brighter-fatter correction. A dict keyed by amplifier name for
            the detector in question.
        defects : `lsst.ip.isr.Defects`, optional
            List of defects.
        fringes : `lsst.pipe.base.Struct`, optional
            Struct containing the fringe correction data, with
            elements:

            ``fringes``
                fringe calibration frame (`lsst.afw.image.Exposure`)
            ``seed``
                random seed derived from the ``ccdExposureId`` for random
                number generator (`numpy.uint32`)
        opticsTransmission: `lsst.afw.image.TransmissionCurve`, optional
            A ``TransmissionCurve`` that represents the throughput of the,
            optics, to be evaluated in focal-plane coordinates.
        filterTransmission : `lsst.afw.image.TransmissionCurve`
            A ``TransmissionCurve`` that represents the throughput of the
            filter itself, to be evaluated in focal-plane coordinates.
        sensorTransmission : `lsst.afw.image.TransmissionCurve`
            A ``TransmissionCurve`` that represents the throughput of the
            sensor itself, to be evaluated in post-assembly trimmed detector
            coordinates.
        atmosphereTransmission : `lsst.afw.image.TransmissionCurve`
            A ``TransmissionCurve`` that represents the throughput of the
            atmosphere, assumed to be spatially constant.
        detectorNum : `int`, optional
            The integer number for the detector to process.
        strayLightData : `object`, optional
            Opaque object containing calibration information for stray-light
            correction.  If `None`, no correction will be performed.
        illumMaskedImage : `lsst.afw.image.MaskedImage`, optional
            Illumination correction image.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with component:

            ``exposure``
                The fully ISR corrected exposure.
                (`lsst.afw.image.Exposure`)
            ``outputExposure``
                An alias for ``exposure``. (`lsst.afw.image.Exposure`)
            ``ossThumb``
                Thumbnail image of the exposure after overscan subtraction.
                (`numpy.ndarray`)
            ``flattenedThumb``
                Thumbnail image of the exposure after flat-field correction.
                (`numpy.ndarray`)
            ``outputStatistics``
                Values of the additional statistics calculated.

        Raises
        ------
        RuntimeError
            Raised if a configuration option is set to `True`, but the
            required calibration data has not been specified.

        Notes
        -----
        The current processed exposure can be viewed by setting the
        appropriate `lsstDebug` entries in the ``debug.display``
        dictionary.  The names of these entries correspond to some of
        the `IsrTaskConfig` Boolean options, with the value denoting the
        frame to use.  The exposure is shown inside the matching
        option check and after the processing of that step has
        finished.  The steps with debug points are:

        * doAssembleCcd
        * doBias
        * doCrosstalk
        * doBrighterFatter
        * doDark
        * doFringe
        * doStrayLight
        * doFlat

        In addition, setting the ``postISRCCD`` entry displays the
        exposure after all ISR processing has finished.
        """

        ccdExposure = self.ensureExposure(ccdExposure, camera, detectorNum)
        bias = self.ensureExposure(bias, camera, detectorNum)
        dark = self.ensureExposure(dark, camera, detectorNum)
        flat = self.ensureExposure(flat, camera, detectorNum)

        ccd = ccdExposure.getDetector()
        filterLabel = ccdExposure.getFilter()
        physicalFilter = isrFunctions.getPhysicalFilter(filterLabel, self.log)

        if not ccd:
            assert not self.config.doAssembleCcd, "You need a Detector to run assembleCcd."
            ccd = [FakeAmp(ccdExposure, self.config)]

        # Validate Input
        if self.config.doBias and bias is None:
            raise RuntimeError("Must supply a bias exposure if config.doBias=True.")
        if self.doLinearize(ccd) and linearizer is None:
            raise RuntimeError("Must supply a linearizer if config.doLinearize=True for this detector.")
        if self.config.doBrighterFatter and bfKernel is None:
            raise RuntimeError("Must supply a kernel if config.doBrighterFatter=True.")
        if self.config.doDark and dark is None:
            raise RuntimeError("Must supply a dark exposure if config.doDark=True.")
        if self.config.doFlat and flat is None:
            raise RuntimeError("Must supply a flat exposure if config.doFlat=True.")
        if self.config.doDefect and defects is None:
            raise RuntimeError("Must supply defects if config.doDefect=True.")
        if (self.config.doFringe and physicalFilter in self.fringe.config.filters
                and fringes.fringes is None):
            # The `fringes` object needs to be a pipeBase.Struct, as
            # we use it as a `dict` for the parameters of
            # `FringeTask.run()`.  The `fringes.fringes` `list` may
            # not be `None` if `doFringe=True`.  Otherwise, raise.
            raise RuntimeError("Must supply fringe exposure as a pipeBase.Struct.")
        if (self.config.doIlluminationCorrection and physicalFilter in self.config.illumFilters
                and illumMaskedImage is None):
            raise RuntimeError("Must supply an illumcor if config.doIlluminationCorrection=True.")
        if (self.config.doDeferredCharge and deferredChargeCalib is None):
            raise RuntimeError("Must supply a deferred charge calibration if config.doDeferredCharge=True.")

        if self.config.doHeaderProvenance:
            # Inputs have been validated, so we can add their date
            # information to the output header.
            exposureMetadata = ccdExposure.getMetadata()
            if self.config.doBias:
                exposureMetadata["LSST CALIB DATE BIAS"] = self.extractCalibDate(bias)
                self.compareCameraKeywords(exposureMetadata, bias, "bias")
            if self.config.doBrighterFatter:
                exposureMetadata["LSST CALIB DATE BFK"] = self.extractCalibDate(bfKernel)
                self.compareCameraKeywords(exposureMetadata, bfKernel, "brighter-fatter")
            if self.config.doCrosstalk:
                exposureMetadata["LSST CALIB DATE CROSSTALK"] = self.extractCalibDate(crosstalk)
                self.compareCameraKeywords(exposureMetadata, crosstalk, "crosstalk")
            if self.config.doDark:
                exposureMetadata["LSST CALIB DATE DARK"] = self.extractCalibDate(dark)
                self.compareCameraKeywords(exposureMetadata, dark, "dark")
            if self.config.doDefect:
                exposureMetadata["LSST CALIB DATE DEFECTS"] = self.extractCalibDate(defects)
                self.compareCameraKeywords(exposureMetadata, defects, "defects")
            if self.config.doDeferredCharge:
                exposureMetadata["LSST CALIB DATE CTI"] = self.extractCalibDate(deferredChargeCalib)
                self.compareCameraKeywords(exposureMetadata, deferredChargeCalib, "CTI")
            if self.config.doFlat:
                exposureMetadata["LSST CALIB DATE FLAT"] = self.extractCalibDate(flat)
                self.compareCameraKeywords(exposureMetadata, flat, "flat")
            if (self.config.doFringe and physicalFilter in self.fringe.config.filters):
                exposureMetadata["LSST CALIB DATE FRINGE"] = self.extractCalibDate(fringes.fringes)
                self.compareCameraKeywords(exposureMetadata, fringes.fringes, "fringe")
            if (self.config.doIlluminationCorrection and physicalFilter in self.config.illumFilters):
                exposureMetadata["LSST CALIB DATE ILLUMINATION"] = self.extractCalibDate(illumMaskedImage)
                self.compareCameraKeywords(exposureMetadata, illumMaskedImage, "illumination")
            if self.doLinearize(ccd):
                exposureMetadata["LSST CALIB DATE LINEARIZER"] = self.extractCalibDate(linearizer)
                self.compareCameraKeywords(exposureMetadata, linearizer, "linearizer")
            if self.config.usePtcGains or self.config.usePtcReadNoise:
                exposureMetadata["LSST CALIB DATE PTC"] = self.extractCalibDate(ptc)
                self.compareCameraKeywords(exposureMetadata, ptc, "PTC")
            if self.config.doStrayLight:
                exposureMetadata["LSST CALIB DATE STRAYLIGHT"] = self.extractCalibDate(strayLightData)
                self.compareCameraKeywords(exposureMetadata, strayLightData, "straylight")
            if self.config.doAttachTransmissionCurve:
                exposureMetadata["LSST CALIB DATE OPTICS_TR"] = self.extractCalibDate(opticsTransmission)
                exposureMetadata["LSST CALIB DATE FILTER_TR"] = self.extractCalibDate(filterTransmission)
                exposureMetadata["LSST CALIB DATE SENSOR_TR"] = self.extractCalibDate(sensorTransmission)
                exposureMetadata["LSST CALIB DATE ATMOSP_TR"] = self.extractCalibDate(atmosphereTransmission)

        # Begin ISR processing.
        if self.config.doConvertIntToFloat:
            self.log.info("Converting exposure to floating point values.")
            ccdExposure = self.convertIntToFloat(ccdExposure)

        if self.config.doBias and self.config.doBiasBeforeOverscan:
            self.log.info("Applying bias correction.")
            isrFunctions.biasCorrection(ccdExposure.getMaskedImage(), bias.getMaskedImage(),
                                        trimToFit=self.config.doTrimToMatchCalib)
            self.debugView(ccdExposure, "doBias")

        # Amplifier level processing.
        overscans = []

        if self.config.doOverscan and self.config.overscan.doParallelOverscan:
            # This will attempt to mask bleed pixels across all amplifiers.
            self.overscan.maskParallelOverscan(ccdExposure, ccd)

        for amp in ccd:
            # if ccdExposure is one amp,
            # check for coverage to prevent performing ops multiple times
            if ccdExposure.getBBox().contains(amp.getBBox()):
                # Check for fully masked bad amplifiers,
                # and generate masks for SUSPECT and SATURATED values.
                badAmp = self.maskAmplifier(ccdExposure, amp, defects)

                if self.config.doOverscan and not badAmp:
                    # Overscan correction on amp-by-amp basis.
                    overscanResults = self.overscanCorrection(ccdExposure, amp)
                    self.log.debug("Corrected overscan for amplifier %s.", amp.getName())
                    if overscanResults is not None and \
                       self.config.qa is not None and self.config.qa.saveStats is True:
                        if isinstance(overscanResults.overscanMean, float):
                            # Only serial overscan was run
                            mean = overscanResults.overscanMean
                            sigma = overscanResults.overscanSigma
                            residMean = overscanResults.residualMean
                            residSigma = overscanResults.residualSigma
                        else:
                            # Both serial and parallel overscan were
                            # run.  Only report serial here.
                            mean = overscanResults.overscanMean[0]
                            sigma = overscanResults.overscanSigma[0]
                            residMean = overscanResults.residualMean[0]
                            residSigma = overscanResults.residualSigma[0]

                        self.metadata[f"FIT MEDIAN {amp.getName()}"] = mean
                        self.metadata[f"FIT STDEV {amp.getName()}"] = sigma
                        self.log.debug("  Overscan stats for amplifer %s: %f +/- %f",
                                       amp.getName(), mean, sigma)

                        self.metadata[f"RESIDUAL MEDIAN {amp.getName()}"] = residMean
                        self.metadata[f"RESIDUAL STDEV {amp.getName()}"] = residSigma
                        self.log.debug("  Overscan stats for amplifer %s after correction: %f +/- %f",
                                       amp.getName(), residMean, residSigma)

                        ccdExposure.getMetadata().set('OVERSCAN', "Overscan corrected")
                else:
                    if badAmp:
                        self.log.warning("Amplifier %s is bad.", amp.getName())
                    overscanResults = None

                overscans.append(overscanResults if overscanResults is not None else None)
            else:
                self.log.info("Skipped OSCAN for %s.", amp.getName())

        if self.config.doDeferredCharge:
            self.log.info("Applying deferred charge/CTI correction.")
            self.deferredChargeCorrection.run(ccdExposure, deferredChargeCalib)
            self.debugView(ccdExposure, "doDeferredCharge")

        if self.config.doCrosstalk and self.config.doCrosstalkBeforeAssemble:
            self.log.info("Applying crosstalk correction.")
            self.crosstalk.run(ccdExposure, crosstalk=crosstalk,
                               crosstalkSources=crosstalkSources, camera=camera)
            self.debugView(ccdExposure, "doCrosstalk")

        if self.config.doAssembleCcd:
            self.log.info("Assembling CCD from amplifiers.")
            ccdExposure = self.assembleCcd.assembleCcd(ccdExposure)

            if self.config.expectWcs and not ccdExposure.getWcs():
                self.log.warning("No WCS found in input exposure.")
            self.debugView(ccdExposure, "doAssembleCcd")

        ossThumb = None
        if self.config.qa.doThumbnailOss:
            ossThumb = isrQa.makeThumbnail(ccdExposure, isrQaConfig=self.config.qa)

        if self.config.doBias and not self.config.doBiasBeforeOverscan:
            self.log.info("Applying bias correction.")
            isrFunctions.biasCorrection(ccdExposure.getMaskedImage(), bias.getMaskedImage(),
                                        trimToFit=self.config.doTrimToMatchCalib)
            self.debugView(ccdExposure, "doBias")

        if self.config.doVariance:
            for amp, overscanResults in zip(ccd, overscans):
                if ccdExposure.getBBox().contains(amp.getBBox()):
                    self.log.debug("Constructing variance map for amplifer %s.", amp.getName())
                    ampExposure = ccdExposure.Factory(ccdExposure, amp.getBBox())
                    if overscanResults is not None:
                        self.updateVariance(ampExposure, amp,
                                            overscanImage=overscanResults.overscanImage,
                                            ptcDataset=ptc)
                    else:
                        self.updateVariance(ampExposure, amp,
                                            overscanImage=None,
                                            ptcDataset=ptc)
                    if self.config.qa is not None and self.config.qa.saveStats is True:
                        qaStats = afwMath.makeStatistics(ampExposure.getVariance(),
                                                         afwMath.MEDIAN | afwMath.STDEVCLIP)
                        self.metadata[f"ISR VARIANCE {amp.getName()} MEDIAN"] = \
                            qaStats.getValue(afwMath.MEDIAN)
                        self.metadata[f"ISR VARIANCE {amp.getName()} STDEV"] = \
                            qaStats.getValue(afwMath.STDEVCLIP)
                        self.log.debug("  Variance stats for amplifer %s: %f +/- %f.",
                                       amp.getName(), qaStats.getValue(afwMath.MEDIAN),
                                       qaStats.getValue(afwMath.STDEVCLIP))
            if self.config.maskNegativeVariance:
                self.maskNegativeVariance(ccdExposure)

        if self.doLinearize(ccd):
            self.log.info("Applying linearizer.")
            linearizer.applyLinearity(image=ccdExposure.getMaskedImage().getImage(),
                                      detector=ccd, log=self.log)

        if self.config.doCrosstalk and not self.config.doCrosstalkBeforeAssemble:
            self.log.info("Applying crosstalk correction.")
            self.crosstalk.run(ccdExposure, crosstalk=crosstalk,
                               crosstalkSources=crosstalkSources, isTrimmed=True)
            self.debugView(ccdExposure, "doCrosstalk")

        # Masking block. Optionally mask known defects, NAN/inf pixels,
        # widen trails, and do anything else the camera needs. Saturated and
        # suspect pixels have already been masked.
        if self.config.doDefect:
            self.log.info("Masking defects.")
            self.maskDefect(ccdExposure, defects)

            if self.config.numEdgeSuspect > 0:
                self.log.info("Masking edges as SUSPECT.")
                self.maskEdges(ccdExposure, numEdgePixels=self.config.numEdgeSuspect,
                               maskPlane="SUSPECT", level=self.config.edgeMaskLevel)

        if self.config.doNanMasking:
            self.log.info("Masking non-finite (NAN, inf) value pixels.")
            self.maskNan(ccdExposure)

        if self.config.doWidenSaturationTrails:
            self.log.info("Widening saturation trails.")
            isrFunctions.widenSaturationTrails(ccdExposure.getMaskedImage().getMask())

        if self.config.doCameraSpecificMasking:
            self.log.info("Masking regions for camera specific reasons.")
            self.masking.run(ccdExposure)

        if self.config.doBrighterFatter:
            # We need to apply flats and darks before we can interpolate, and
            # we need to interpolate before we do B-F, but we do B-F without
            # the flats and darks applied so we can work in units of electrons
            # or holes. This context manager applies and then removes the darks
            # and flats.
            #
            # We also do not want to interpolate values here, so operate on
            # temporary images so we can apply only the BF-correction and roll
            # back the interpolation.
            interpExp = ccdExposure.clone()
            with self.flatContext(interpExp, flat, dark):
                isrFunctions.interpolateFromMask(
                    maskedImage=interpExp.getMaskedImage(),
                    fwhm=self.config.fwhm,
                    growSaturatedFootprints=self.config.growSaturationFootprintSize,
                    maskNameList=list(self.config.brighterFatterMaskListToInterpolate)
                )
            bfExp = interpExp.clone()

            self.log.info("Applying brighter-fatter correction using kernel type %s / gains %s.",
                          type(bfKernel), type(bfGains))
            if self.config.doFluxConservingBrighterFatterCorrection:
                bfResults = isrFunctions.fluxConservingBrighterFatterCorrection(
                    bfExp,
                    bfKernel,
                    self.config.brighterFatterMaxIter,
                    self.config.brighterFatterThreshold,
                    self.config.brighterFatterApplyGain,
                    bfGains
                )
            else:
                bfResults = isrFunctions.brighterFatterCorrection(
                    bfExp,
                    bfKernel,
                    self.config.brighterFatterMaxIter,
                    self.config.brighterFatterThreshold,
                    self.config.brighterFatterApplyGain,
                    bfGains
                )
            if bfResults[1] == self.config.brighterFatterMaxIter-1:
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

            self.debugView(ccdExposure, "doBrighterFatter")

        if self.config.doDark:
            self.log.info("Applying dark correction.")
            self.darkCorrection(ccdExposure, dark)
            self.debugView(ccdExposure, "doDark")

        if self.config.doFringe and not self.config.fringeAfterFlat:
            self.log.info("Applying fringe correction before flat.")
            self.fringe.run(ccdExposure, **fringes.getDict())
            self.debugView(ccdExposure, "doFringe")

        if self.config.doStrayLight and self.strayLight.check(ccdExposure):
            self.log.info("Checking strayLight correction.")
            self.strayLight.run(ccdExposure, strayLightData)
            self.debugView(ccdExposure, "doStrayLight")

        if self.config.doFlat:
            self.log.info("Applying flat correction.")
            self.flatCorrection(ccdExposure, flat)
            self.debugView(ccdExposure, "doFlat")

        if self.config.doApplyGains:
            self.log.info("Applying gain correction instead of flat.")
            if self.config.usePtcGains:
                self.log.info("Using gains from the Photon Transfer Curve.")
                isrFunctions.applyGains(ccdExposure, self.config.normalizeGains,
                                        ptcGains=ptc.gain)
            else:
                isrFunctions.applyGains(ccdExposure, self.config.normalizeGains)

        if self.config.doFringe and self.config.fringeAfterFlat:
            self.log.info("Applying fringe correction after flat.")
            self.fringe.run(ccdExposure, **fringes.getDict())

        if self.config.doVignette:
            if self.config.doMaskVignettePolygon:
                self.log.info("Constructing, attaching, and masking vignette polygon.")
            else:
                self.log.info("Constructing and attaching vignette polygon.")
            self.vignettePolygon = self.vignette.run(
                exposure=ccdExposure, doUpdateMask=self.config.doMaskVignettePolygon,
                vignetteValue=self.config.vignetteValue, log=self.log)

        if self.config.doAttachTransmissionCurve:
            self.log.info("Adding transmission curves.")
            isrFunctions.attachTransmissionCurve(ccdExposure, opticsTransmission=opticsTransmission,
                                                 filterTransmission=filterTransmission,
                                                 sensorTransmission=sensorTransmission,
                                                 atmosphereTransmission=atmosphereTransmission)

        flattenedThumb = None
        if self.config.qa.doThumbnailFlattened:
            flattenedThumb = isrQa.makeThumbnail(ccdExposure, isrQaConfig=self.config.qa)

        if self.config.doIlluminationCorrection and physicalFilter in self.config.illumFilters:
            self.log.info("Performing illumination correction.")
            isrFunctions.illuminationCorrection(ccdExposure.getMaskedImage(),
                                                illumMaskedImage, illumScale=self.config.illumScale,
                                                trimToFit=self.config.doTrimToMatchCalib)

        preInterpExp = None
        if self.config.doSaveInterpPixels:
            preInterpExp = ccdExposure.clone()

        # Reset and interpolate bad pixels.
        #
        # Large contiguous bad regions (which should have the BAD mask
        # bit set) should have their values set to the image median.
        # This group should include defects and bad amplifiers. As the
        # area covered by these defects are large, there's little
        # reason to expect that interpolation would provide a more
        # useful value.
        #
        # Smaller defects can be safely interpolated after the larger
        # regions have had their pixel values reset.  This ensures
        # that the remaining defects adjacent to bad amplifiers (as an
        # example) do not attempt to interpolate extreme values.
        if self.config.doSetBadRegions:
            badPixelCount, badPixelValue = isrFunctions.setBadRegions(ccdExposure)
            if badPixelCount > 0:
                self.log.info("Set %d BAD pixels to %f.", badPixelCount, badPixelValue)

        if self.config.doInterpolate:
            self.log.info("Interpolating masked pixels.")
            isrFunctions.interpolateFromMask(
                maskedImage=ccdExposure.getMaskedImage(),
                fwhm=self.config.fwhm,
                growSaturatedFootprints=self.config.growSaturationFootprintSize,
                maskNameList=list(self.config.maskListToInterpolate)
            )

        self.roughZeroPoint(ccdExposure)

        # correct for amp offsets within the CCD
        if self.config.doAmpOffset:
            self.log.info("Correcting amp offsets.")
            self.ampOffset.run(ccdExposure)

        if self.config.doMeasureBackground:
            self.log.info("Measuring background level.")
            self.measureBackground(ccdExposure, self.config.qa)

            if self.config.qa is not None and self.config.qa.saveStats is True:
                for amp in ccd:
                    ampExposure = ccdExposure.Factory(ccdExposure, amp.getBBox())
                    qaStats = afwMath.makeStatistics(ampExposure.getImage(),
                                                     afwMath.MEDIAN | afwMath.STDEVCLIP)
                    self.metadata[f"ISR BACKGROUND {amp.getName()} MEDIAN"] = qaStats.getValue(afwMath.MEDIAN)
                    self.metadata[f"ISR BACKGROUND {amp.getName()} STDEV"] = \
                        qaStats.getValue(afwMath.STDEVCLIP)
                    self.log.debug("  Background stats for amplifer %s: %f +/- %f",
                                   amp.getName(), qaStats.getValue(afwMath.MEDIAN),
                                   qaStats.getValue(afwMath.STDEVCLIP))

        # Calculate standard image quality statistics
        if self.config.doStandardStatistics:
            metadata = ccdExposure.getMetadata()
            for amp in ccd:
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
            outputStatistics = self.isrStats.run(ccdExposure, overscanResults=overscans,
                                                 ptc=ptc).results

        self.debugView(ccdExposure, "postISRCCD")

        return pipeBase.Struct(
            exposure=ccdExposure,
            ossThumb=ossThumb,
            flattenedThumb=flattenedThumb,

            preInterpExposure=preInterpExp,
            outputExposure=ccdExposure,
            outputOssThumbnail=ossThumb,
            outputFlattenedThumbnail=flattenedThumb,
            outputStatistics=outputStatistics,
        )

    def ensureExposure(self, inputExp, camera=None, detectorNum=None):
        """Ensure that the data returned by Butler is a fully constructed exp.

        ISR requires exposure-level image data for historical reasons, so if we
        did not recieve that from Butler, construct it from what we have,
        modifying the input in place.

        Parameters
        ----------
        inputExp : `lsst.afw.image` image-type.
            The input data structure obtained from Butler.
            Can be  `lsst.afw.image.Exposure`,
            `lsst.afw.image.DecoratedImageU`,
            or `lsst.afw.image.ImageF`
        camera : `lsst.afw.cameraGeom.camera`, optional
            The camera associated with the image.  Used to find the appropriate
            detector if detector is not already set.
        detectorNum : `int`, optional
            The detector in the camera to attach, if the detector is not
            already set.

        Returns
        -------
        inputExp : `lsst.afw.image.Exposure`
            The re-constructed exposure, with appropriate detector parameters.

        Raises
        ------
        TypeError
            Raised if the input data cannot be used to construct an exposure.
        """
        if isinstance(inputExp, afwImage.DecoratedImageU):
            inputExp = afwImage.makeExposure(afwImage.makeMaskedImage(inputExp))
        elif isinstance(inputExp, afwImage.ImageF):
            inputExp = afwImage.makeExposure(afwImage.makeMaskedImage(inputExp))
        elif isinstance(inputExp, afwImage.MaskedImageF):
            inputExp = afwImage.makeExposure(inputExp)
        elif isinstance(inputExp, afwImage.Exposure):
            pass
        elif inputExp is None:
            # Assume this will be caught by the setup if it is a problem.
            return inputExp
        else:
            raise TypeError("Input Exposure is not known type in isrTask.ensureExposure: %s." %
                            (type(inputExp), ))

        if inputExp.getDetector() is None:
            if camera is None or detectorNum is None:
                raise RuntimeError('Must supply both a camera and detector number when using exposures '
                                   'without a detector set.')
            inputExp.setDetector(camera[detectorNum])

        return inputExp

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

    def maskAmplifier(self, ccdExposure, amp, defects):
        """Identify bad amplifiers, saturated and suspect pixels.

        Parameters
        ----------
        ccdExposure : `lsst.afw.image.Exposure`
            Input exposure to be masked.
        amp : `lsst.afw.cameraGeom.Amplifier`
            Catalog of parameters defining the amplifier on this
            exposure to mask.
        defects : `lsst.ip.isr.Defects`
            List of defects.  Used to determine if the entire
            amplifier is bad.

        Returns
        -------
        badAmp : `Bool`
            If this is true, the entire amplifier area is covered by
            defects and unusable.

        """
        maskedImage = ccdExposure.getMaskedImage()

        badAmp = False

        # Check if entire amp region is defined as a defect
        # NB: need to use amp.getBBox() for correct comparison with current
        # defects definition.
        if defects is not None:
            badAmp = bool(sum([v.getBBox().contains(amp.getBBox()) for v in defects]))

        # In the case of a bad amp, we will set mask to "BAD"
        # (here use amp.getRawBBox() for correct association with pixels in
        # current ccdExposure).
        if badAmp:
            dataView = afwImage.MaskedImageF(maskedImage, amp.getRawBBox(),
                                             afwImage.PARENT)
            maskView = dataView.getMask()
            maskView |= maskView.getPlaneBitMask("BAD")
            del maskView
            return badAmp

        # Mask remaining defects after assembleCcd() to allow for defects that
        # cross amplifier boundaries. Saturation and suspect pixels can be
        # masked now, though.
        limits = dict()
        if self.config.doSaturation and not badAmp:
            limits.update({self.config.saturatedMaskName: amp.getSaturation()})
        if self.config.doSuspect and not badAmp:
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
            badAmp = True
            maskView |= maskView.getPlaneBitMask("BAD")

        return badAmp

    def overscanCorrection(self, ccdExposure, amp):
        """Apply overscan correction in place.

        This method does initial pixel rejection of the overscan
        region.  The overscan can also be optionally segmented to
        allow for discontinuous overscan responses to be fit
        separately.  The actual overscan subtraction is performed by
        the `lsst.ip.isr.overscan.OverscanTask`, which is called here
        after the amplifier is preprocessed.

        Parameters
        ----------
        ccdExposure : `lsst.afw.image.Exposure`
            Exposure to have overscan correction performed.
        amp : `lsst.afw.cameraGeom.Amplifer`
            The amplifier to consider while correcting the overscan.

        Returns
        -------
        overscanResults : `lsst.pipe.base.Struct`
            Result struct with components:

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
            ``edgeMask``
                Mask of the suspect pixels. (`lsst.afw.image.Mask`)
            ``overscanMean``
                Median overscan fit value. (`float`)
            ``overscanSigma``
                Clipped standard deviation of the overscan after
                correction. (`float`)

        Raises
        ------
        RuntimeError
            Raised if the ``amp`` does not contain raw pixel information.

        See Also
        --------
        lsst.ip.isr.overscan.OverscanTask
        """
        if amp.getRawHorizontalOverscanBBox().isEmpty():
            self.log.info("ISR_OSCAN: No overscan region.  Not performing overscan correction.")
            return None

        # Perform overscan correction on subregions.
        overscanResults = self.overscan.run(ccdExposure, amp)

        metadata = ccdExposure.getMetadata()
        ampName = amp.getName()

        keyBase = "LSST ISR OVERSCAN"
        # Updated quantities
        if isinstance(overscanResults.overscanMean, float):
            # Serial overscan correction only:
            metadata[f"{keyBase} SERIAL MEAN {ampName}"] = overscanResults.overscanMean
            metadata[f"{keyBase} SERIAL MEDIAN {ampName}"] = overscanResults.overscanMedian
            metadata[f"{keyBase} SERIAL STDEV {ampName}"] = overscanResults.overscanSigma

            metadata[f"{keyBase} RESIDUAL SERIAL MEAN {ampName}"] = overscanResults.residualMean
            metadata[f"{keyBase} RESIDUAL SERIAL MEDIAN {ampName}"] = overscanResults.residualMedian
            metadata[f"{keyBase} RESIDUAL SERIAL STDEV {ampName}"] = overscanResults.residualSigma
        elif isinstance(overscanResults.overscanMean, tuple):
            # Both serial and parallel overscan have run:
            metadata[f"{keyBase} SERIAL MEAN {ampName}"] = overscanResults.overscanMean[0]
            metadata[f"{keyBase} SERIAL MEDIAN {ampName}"] = overscanResults.overscanMedian[0]
            metadata[f"{keyBase} SERIAL STDEV {ampName}"] = overscanResults.overscanSigma[0]

            metadata[f"{keyBase} PARALLEL MEAN {ampName}"] = overscanResults.overscanMean[1]
            metadata[f"{keyBase} PARALLEL MEDIAN {ampName}"] = overscanResults.overscanMedian[1]
            metadata[f"{keyBase} PARALLEL STDEV {ampName}"] = overscanResults.overscanSigma[1]

            metadata[f"{keyBase} RESIDUAL SERIAL MEAN {ampName}"] = overscanResults.residualMean[0]
            metadata[f"{keyBase} RESIDUAL SERIAL MEDIAN {ampName}"] = overscanResults.residualMedian[0]
            metadata[f"{keyBase} RESIDUAL SERIAL STDEV {ampName}"] = overscanResults.residualSigma[0]

            metadata[f"{keyBase} RESIDUAL PARALLEL MEAN {ampName}"] = overscanResults.residualMean[1]
            metadata[f"{keyBase} RESIDUAL PARALLEL MEDIAN {ampName}"] = overscanResults.residualMedian[1]
            metadata[f"{keyBase} RESIDUAL PARALLEL STDEV {ampName}"] = overscanResults.residualSigma[1]
        else:
            self.log.warning("Unexpected type for overscan values; none added to header.")

        return overscanResults

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
            Raised if either ``usePtcGains`` of ``usePtcReadNoise``
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
            trimToFit=self.config.doTrimToMatchCalib
        )

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
            invert=invert,
            trimToFit=self.config.doTrimToMatchCalib
        )

    def saturationDetection(self, exposure, amp):
        """Detect and mask saturated pixels in config.saturatedMaskName.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to process.  Only the amplifier DataSec is processed.
        amp : `lsst.afw.cameraGeom.Amplifier`
            Amplifier detector data.

        See Also
        --------
        lsst.ip.isr.isrFunctions.makeThresholdMask
        """
        if not math.isnan(amp.getSaturation()):
            maskedImage = exposure.getMaskedImage()
            dataView = maskedImage.Factory(maskedImage, amp.getRawBBox())
            isrFunctions.makeThresholdMask(
                maskedImage=dataView,
                threshold=amp.getSaturation(),
                growFootprints=0,
                maskName=self.config.saturatedMaskName,
            )

    def saturationInterpolation(self, exposure):
        """Interpolate over saturated pixels, in place.

        This method should be called after `saturationDetection`, to
        ensure that the saturated pixels have been identified in the
        SAT mask.  It should also be called after `assembleCcd`, since
        saturated regions may cross amplifier boundaries.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to process.

        See Also
        --------
        lsst.ip.isr.isrTask.saturationDetection
        lsst.ip.isr.isrFunctions.interpolateFromMask
        """
        isrFunctions.interpolateFromMask(
            maskedImage=exposure.getMaskedImage(),
            fwhm=self.config.fwhm,
            growSaturatedFootprints=self.config.growSaturationFootprintSize,
            maskNameList=list(self.config.saturatedMaskName),
        )

    def suspectDetection(self, exposure, amp):
        """Detect and mask suspect pixels in config.suspectMaskName.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to process.  Only the amplifier DataSec is processed.
        amp : `lsst.afw.cameraGeom.Amplifier`
            Amplifier detector data.

        See Also
        --------
        lsst.ip.isr.isrFunctions.makeThresholdMask

        Notes
        -----
        Suspect pixels are pixels whose value is greater than
        amp.getSuspectLevel(). This is intended to indicate pixels that may be
        affected by unknown systematics; for example if non-linearity
        corrections above a certain level are unstable then that would be a
        useful value for suspectLevel. A value of `nan` indicates that no such
        level exists and no pixels are to be masked as suspicious.
        """
        suspectLevel = amp.getSuspectLevel()
        if math.isnan(suspectLevel):
            return

        maskedImage = exposure.getMaskedImage()
        dataView = maskedImage.Factory(maskedImage, amp.getRawBBox())
        isrFunctions.makeThresholdMask(
            maskedImage=dataView,
            threshold=suspectLevel,
            growFootprints=0,
            maskName=self.config.suspectMaskName,
        )

    def maskDefect(self, exposure, defectBaseList):
        """Mask defects using mask plane "BAD", in place.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to process.
        defectBaseList : defect-type
            List of defects to mask. Can be of type  `lsst.ip.isr.Defects`
            or `list` of `lsst.afw.image.DefectBase`.

        Notes
        -----
        Call this after CCD assembly, since defects may cross amplifier
        boundaries.
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

    def maskAndInterpolateDefects(self, exposure, defectBaseList):
        """Mask and interpolate defects using mask plane "BAD", in place.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to process.
        defectBaseList : defects-like
            List of defects to mask and interpolate. Can be
            `lsst.ip.isr.Defects` or `list` of `lsst.afw.image.DefectBase`.

        See Also
        --------
        lsst.ip.isr.isrTask.maskDefect
        """
        self.maskDefect(exposure, defectBaseList)
        self.maskEdges(exposure, numEdgePixels=self.config.numEdgeSuspect,
                       maskPlane="SUSPECT", level=self.config.edgeMaskLevel)
        isrFunctions.interpolateFromMask(
            maskedImage=exposure.getMaskedImage(),
            fwhm=self.config.fwhm,
            growSaturatedFootprints=0,
            maskNameList=["BAD"],
        )

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

    def maskAndInterpolateNan(self, exposure):
        """"Mask and interpolate NaN/infs using mask plane "UNMASKEDNAN",
        in place.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to process.

        See Also
        --------
        lsst.ip.isr.isrTask.maskNan
        """
        self.maskNan(exposure)
        isrFunctions.interpolateFromMask(
            maskedImage=exposure.getMaskedImage(),
            fwhm=self.config.fwhm,
            growSaturatedFootprints=0,
            maskNameList=["UNMASKEDNAN"],
        )

    def measureBackground(self, exposure, IsrQaConfig=None):
        """Measure the image background in subgrids, for quality control.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to process.
        IsrQaConfig : `lsst.ip.isr.isrQa.IsrQaConfig`
            Configuration object containing parameters on which background
            statistics and subgrids to use.
        """
        if IsrQaConfig is not None:
            statsControl = afwMath.StatisticsControl(IsrQaConfig.flatness.clipSigma,
                                                     IsrQaConfig.flatness.nIter)
            maskVal = exposure.getMaskedImage().getMask().getPlaneBitMask(["BAD", "SAT", "DETECTED"])
            statsControl.setAndMask(maskVal)
            maskedImage = exposure.getMaskedImage()
            stats = afwMath.makeStatistics(maskedImage, afwMath.MEDIAN | afwMath.STDEVCLIP, statsControl)
            skyLevel = stats.getValue(afwMath.MEDIAN)
            skySigma = stats.getValue(afwMath.STDEVCLIP)
            self.log.info("Flattened sky level: %f +/- %f.", skyLevel, skySigma)
            metadata = exposure.getMetadata()
            metadata["SKYLEVEL"] = skyLevel
            metadata["SKYSIGMA"] = skySigma

            # calcluating flatlevel over the subgrids
            stat = afwMath.MEANCLIP if IsrQaConfig.flatness.doClip else afwMath.MEAN
            meshXHalf = int(IsrQaConfig.flatness.meshX/2.)
            meshYHalf = int(IsrQaConfig.flatness.meshY/2.)
            nX = int((exposure.getWidth() + meshXHalf) / IsrQaConfig.flatness.meshX)
            nY = int((exposure.getHeight() + meshYHalf) / IsrQaConfig.flatness.meshY)
            skyLevels = numpy.zeros((nX, nY))

            for j in range(nY):
                yc = meshYHalf + j * IsrQaConfig.flatness.meshY
                for i in range(nX):
                    xc = meshXHalf + i * IsrQaConfig.flatness.meshX

                    xLLC = xc - meshXHalf
                    yLLC = yc - meshYHalf
                    xURC = xc + meshXHalf - 1
                    yURC = yc + meshYHalf - 1

                    bbox = lsst.geom.Box2I(lsst.geom.Point2I(xLLC, yLLC), lsst.geom.Point2I(xURC, yURC))
                    miMesh = maskedImage.Factory(exposure.getMaskedImage(), bbox, afwImage.LOCAL)

                    skyLevels[i, j] = afwMath.makeStatistics(miMesh, stat, statsControl).getValue()

            good = numpy.where(numpy.isfinite(skyLevels))
            skyMedian = numpy.median(skyLevels[good])
            flatness = (skyLevels[good] - skyMedian) / skyMedian
            flatness_rms = numpy.std(flatness)
            flatness_pp = flatness.max() - flatness.min() if len(flatness) > 0 else numpy.nan

            self.log.info("Measuring sky levels in %dx%d grids: %f.", nX, nY, skyMedian)
            self.log.info("Sky flatness in %dx%d grids - pp: %f rms: %f.",
                          nX, nY, flatness_pp, flatness_rms)

            metadata["FLATNESS_PP"] = float(flatness_pp)
            metadata["FLATNESS_RMS"] = float(flatness_rms)
            metadata["FLATNESS_NGRIDS"] = '%dx%d' % (nX, nY)
            metadata["FLATNESS_MESHX"] = IsrQaConfig.flatness.meshX
            metadata["FLATNESS_MESHY"] = IsrQaConfig.flatness.meshY

    def roughZeroPoint(self, exposure):
        """Set an approximate magnitude zero point for the exposure.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to process.
        """
        filterLabel = exposure.getFilter()
        physicalFilter = isrFunctions.getPhysicalFilter(filterLabel, self.log)

        if physicalFilter in self.config.fluxMag0T1:
            fluxMag0 = self.config.fluxMag0T1[physicalFilter]
        else:
            self.log.warning("No rough magnitude zero point defined for filter %s.", physicalFilter)
            fluxMag0 = self.config.defaultFluxMag0T1

        expTime = exposure.getInfo().getVisitInfo().getExposureTime()
        if not expTime > 0:  # handle NaN as well as <= 0
            self.log.warning("Non-positive exposure time; skipping rough zero point.")
            return

        self.log.info("Setting rough magnitude zero point for filter %s: %f",
                      physicalFilter, 2.5*math.log10(fluxMag0*expTime))
        exposure.setPhotoCalib(afwImage.makePhotoCalibFromCalibZeroPoint(fluxMag0*expTime, 0.0))

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

    def debugView(self, exposure, stepname):
        """Utility function to examine ISR exposure at different stages.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to view.
        stepname : `str`
            State of processing to view.
        """
        frame = getDebugFrame(self._display, stepname)
        if frame:
            display = getDisplay(frame)
            display.scale('asinh', 'zscale')
            display.mtv(exposure)
            prompt = "Press Enter to continue [c]... "
            while True:
                ans = input(prompt).lower()
                if ans in ("", "c",):
                    break


class FakeAmp(object):
    """A Detector-like object that supports returning gain and saturation level

    This is used when the input exposure does not have a detector.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure to generate a fake amplifier for.
    config : `lsst.ip.isr.isrTaskConfig`
        Configuration to apply to the fake amplifier.
    """

    def __init__(self, exposure, config):
        self._bbox = exposure.getBBox(afwImage.LOCAL)
        self._RawHorizontalOverscanBBox = lsst.geom.Box2I()
        self._gain = config.gain
        self._readNoise = config.readNoise
        self._saturation = config.saturation

    def getBBox(self):
        return self._bbox

    def getRawBBox(self):
        return self._bbox

    def getRawHorizontalOverscanBBox(self):
        return self._RawHorizontalOverscanBBox

    def getGain(self):
        return self._gain

    def getReadNoise(self):
        return self._readNoise

    def getSaturation(self):
        return self._saturation

    def getSuspectLevel(self):
        return float("NaN")
