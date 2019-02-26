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

import math
import numpy

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.meas.algorithms as measAlg
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from contextlib import contextmanager
from lsstDebug import getDebugFrame

from lsst.afw.cameraGeom import PIXELS, FOCAL_PLANE, NullLinearityType
from lsst.afw.display import getDisplay
from lsst.afw.geom import Polygon
from lsst.daf.persistence import ButlerDataRef
from lsst.meas.algorithms.detection import SourceDetectionTask
from lsst.meas.algorithms import Defect

from . import isrFunctions
from . import isrQa
from . import linearize

from .assembleCcdTask import AssembleCcdTask
from .crosstalk import CrosstalkTask
from .fringe import FringeTask
from .isr import maskNans
from .masking import MaskingTask
from .straylight import StrayLightTask
from .vignette import VignetteTask

__all__ = ["IsrTask", "RunIsrTask"]


class IsrTaskConfig(pexConfig.Config):
    """Configuration parameters for IsrTask.

    Items are grouped in the order in which they are executed by the task.
    """
    # General ISR configuration

    # gen3 options
    isrName = pexConfig.Field(
        dtype=str,
        doc="Name of ISR",
        default="ISR",
    )

    # input datasets
    ccdExposure = pipeBase.InputDatasetField(
        doc="Input exposure to process",
        name="raw",
        scalar=True,
        storageClass="ExposureU",
        dimensions=["Instrument", "Exposure", "Detector"],
    )
    camera = pipeBase.InputDatasetField(
        doc="Input camera to construct complete exposures.",
        name="camera",
        scalar=True,
        storageClass="TablePersistableCamera",
        dimensions=["Instrument", "ExposureRange"],
    )
    bias = pipeBase.InputDatasetField(
        doc="Input bias calibration.",
        name="bias",
        scalar=True,
        storageClass="ImageF",
        dimensions=["Instrument", "ExposureRange", "Detector"],
    )
    dark = pipeBase.InputDatasetField(
        doc="Input dark calibration.",
        name="dark",
        scalar=True,
        storageClass="ImageF",
        dimensions=["Instrument", "ExposureRange", "Detector"],
    )
    flat = pipeBase.InputDatasetField(
        doc="Input flat calibration.",
        name="flat",
        scalar=True,
        storageClass="MaskedImageF",
        dimensions=["Instrument", "PhysicalFilter", "ExposureRange", "Detector"],
    )
    bfKernel = pipeBase.InputDatasetField(
        doc="Input brighter-fatter kernel.",
        name="bfKernel",
        scalar=True,
        storageClass="NumpyArray",
        dimensions=["Instrument", "ExposureRange"],
    )
    defects = pipeBase.InputDatasetField(
        doc="Input defect tables.",
        name="defects",
        scalar=True,
        storageClass="Catalog",
        dimensions=["Instrument", "ExposureRange", "Detector"],
    )
    opticsTransmission = pipeBase.InputDatasetField(
        doc="Transmission curve due to the optics.",
        name="transmission_optics",
        scalar=True,
        storageClass="TablePersistableTransmissionCurve",
        dimensions=["Instrument", "ExposureRange"],
    )
    filterTransmission = pipeBase.InputDatasetField(
        doc="Transmission curve due to the filter.",
        name="transmission_filter",
        scalar=True,
        storageClass="TablePersistableTransmissionCurve",
        dimensions=["Instrument", "PhysicalFilter", "ExposureRange"],
    )
    sensorTransmission = pipeBase.InputDatasetField(
        doc="Transmission curve due to the sensor.",
        name="transmission_sensor",
        scalar=True,
        storageClass="TablePersistableTransmissionCurve",
        dimensions=["Instrument", "ExposureRange", "Detector"],
    )
    atmosphereTransmission = pipeBase.InputDatasetField(
        doc="Transmission curve due to the atmosphere.",
        name="transmission_atmosphere",
        scalar=True,
        storageClass="TablePersistableTransmissionCurve",
        dimensions=["Instrument"],
    )

    # output datasets
    outputExposure = pipeBase.OutputDatasetField(
        doc="Output ISR processed exposure.",
        name="postISRCCD",
        scalar=True,
        storageClass="ExposureF",
        dimensions=["Instrument", "Visit", "Detector"],
    )
    outputOssThumbnail = pipeBase.OutputDatasetField(
        doc="Output Overscan-subtracted thumbnail image.",
        name="OssThumb",
        scalar=True,
        storageClass="Thumbnail",
        dimensions=["Instrument", "Visit", "Detector"],
    )
    outputFlattenedThumbnail = pipeBase.OutputDatasetField(
        doc="Output flat-corrected thumbnail image.",
        name="FlattenedThumb",
        scalar=True,
        storageClass="TextStorage",
        dimensions=["Instrument", "Visit", "Detector"],
    )

    quantum = pipeBase.QuantumConfig(
        dimensions=["Visit", "Detector", "Instrument"],
    )

    ## Original gen2 config
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
    expectWcs = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Expect input science images to have a WCS (set False for e.g. spectrographs)."
    )
    fwhm = pexConfig.Field(
        dtype=float,
        doc="FWHM of PSF in arcseconds.",
        default=1.0,
    )
    qa = pexConfig.ConfigField(
        dtype=isrQa.IsrQaConfig,
        doc="QA related configuration options.",
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
        default=True,
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
    overscanFitType = pexConfig.ChoiceField(
        dtype=str,
        doc="The method for fitting the overscan bias level.",
        default='MEDIAN',
        allowed={
            "POLY": "Fit ordinary polynomial to the longest axis of the overscan region",
            "CHEB": "Fit Chebyshev polynomial to the longest axis of the overscan region",
            "LEG": "Fit Legendre polynomial to the longest axis of the overscan region",
            "NATURAL_SPLINE": "Fit natural spline to the longest axis of the overscan region",
            "CUBIC_SPLINE": "Fit cubic spline to the longest axis of the overscan region",
            "AKIMA_SPLINE": "Fit Akima spline to the longest axis of the overscan region",
            "MEAN": "Correct using the mean of the overscan region",
            "MEANCLIP": "Correct using a clipped mean of the overscan region",
            "MEDIAN": "Correct using the median of the overscan region",
        },
    )
    overscanOrder = pexConfig.Field(
        dtype=int,
        doc=("Order of polynomial or to fit if overscan fit type is a polynomial, " +
             "or number of spline knots if overscan fit type is a spline."),
        default=1,
    )
    overscanNumSigmaClip = pexConfig.Field(
        dtype=float,
        doc="Rejection threshold (sigma) for collapsing overscan before fit",
        default=3.0,
    )
    overscanIsInt = pexConfig.Field(
        dtype=bool,
        doc="Treat overscan as an integer image for purposes of overscan.FitType=MEDIAN",
        default=True,
    )
    overscanNumLeadingColumnsToSkip = pexConfig.Field(
        dtype=int,
        doc="Number of columns to skip in overscan, i.e. those closest to amplifier",
        default=0,
    )
    overscanNumTrailingColumnsToSkip = pexConfig.Field(
        dtype=int,
        doc="Number of columns to skip in overscan, i.e. those farthest from amplifier",
        default=0,
    )
    overscanMaxDev = pexConfig.Field(
        dtype=float,
        doc="Maximum deviation from the median for overscan",
        default=1000.0, check=lambda x: x > 0
    )
    overscanBiasJump = pexConfig.Field(
        dtype=bool,
        doc="Fit the overscan in a piecewise-fashion to correct for bias jumps?",
        default=False,
    )
    overscanBiasJumpKeyword = pexConfig.Field(
        dtype=str,
        doc="Header keyword containing information about devices.",
        default="NO_SUCH_KEY",
    )
    overscanBiasJumpDevices = pexConfig.ListField(
        dtype=str,
        doc="List of devices that need piecewise overscan correction.",
        default=(),
    )
    overscanBiasJumpLocation = pexConfig.Field(
        dtype=int,
        doc="Location of bias jump along y-axis.",
        default=0,
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
        default=True,
    )
    crosstalk = pexConfig.ConfigurableField(
        target=CrosstalkTask,
        doc="Intra-CCD crosstalk correction",
    )

    # Masking option prior to brighter-fatter.
    doWidenSaturationTrails = pexConfig.Field(
        dtype=bool,
        doc="Widen bleed trails based on their width?",
        default=True
    )

    # Brighter-Fatter correction.
    doBrighterFatter = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Apply the brighter fatter correction"
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
    brighterFatterKernelFile = pexConfig.Field(
        dtype=str,
        default='',
        doc="Kernel file used for the brighter fatter correction"
    )
    brighterFatterMaxIter = pexConfig.Field(
        dtype=int,
        default=10,
        doc="Maximum number of iterations for the brighter fatter correction"
    )
    brighterFatterThreshold = pexConfig.Field(
        dtype=float,
        default=1000,
        doc="Threshold used to stop iterating the brighter fatter correction.  It is the "
        " absolute value of the difference between the current corrected image and the one"
        " from the previous iteration summed over all the pixels."
    )
    brighterFatterApplyGain = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Should the gain be applied when applying the brighter fatter correction?"
    )

    # Defect and bad pixel correction options.
    doDefect = pexConfig.Field(
        dtype=bool,
        doc="Apply correction for CCD defects, e.g. hot pixels?",
        default=True,
    )
    doSaturationInterpolation = pexConfig.Field(
        dtype=bool,
        doc="Perform interpolation over pixels masked as saturated?"
        " NB: This is independent of doSaturation; if that is False this plane"
        " will likely be blank, resulting in a no-op here.",
        default=True,
    )
    numEdgeSuspect = pexConfig.Field(
        dtype=int,
        doc="Number of edge pixels to be flagged as untrustworthy.",
        default=0,
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

    # Amplifier normalization based on gains instead of using flats configuration.
    doApplyGains = pexConfig.Field(
        dtype=bool,
        doc="Correct the amplifiers for their gains instead of applying flat correction",
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

    # NAN pixel interpolation option.
    doNanInterpAfterFlat = pexConfig.Field(
        dtype=bool,
        doc=("If True, ensure we interpolate NaNs after flat-fielding, even if we "
             "also have to interpolate them before flat-fielding."),
        default=False,
    )

    # Distortion model application.
    doAddDistortionModel = pexConfig.Field(
        dtype=bool,
        doc="Apply a distortion model based on camera geometry to the WCS?",
        default=True,
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
        doc="Apply vignetting parameters?",
        default=False,
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

    # Write the outputs to disk.  If ISR is run as a subtask, this may not be needed.
    doWrite = pexConfig.Field(
        dtype=bool,
        doc="Persist postISRCCD?",
        default=True,
    )

    def validate(self):
        super().validate()
        if self.doFlat and self.doApplyGains:
            raise ValueError("You may not specify both doFlat and doApplyGains")


class IsrTask(pipeBase.PipelineTask, pipeBase.CmdLineTask):
    r"""Apply common instrument signature correction algorithms to a raw frame.

    The process for correcting imaging data is very similar from
    camera to camera.  This task provides a vanilla implementation of
    doing these corrections, including the ability to turn certain
    corrections off if they are not needed.  The inputs to the primary
    method, `run()`, are a raw exposure to be corrected and the
    calibration data products. The raw input is a single chip sized
    mosaic of all amps including overscans and other non-science
    pixels.  The method `runDataRef()` identifies and defines the
    calibration data products, and is intended for use by a
    `lsst.pipe.base.cmdLineTask.CmdLineTask` and takes as input only a
    `daf.persistence.butlerSubset.ButlerDataRef`.  This task may be
    subclassed for different camera, although the most camera specific
    methods have been split into subtasks that can be redirected
    appropriately.

    The __init__ method sets up the subtasks for ISR processing, using
    the defaults from `lsst.ip.isr`.

    Parameters
    ----------
    args : `list`
        Positional arguments passed to the Task constructor. None used at this time.
    kwargs : `dict`, optional
        Keyword arguments passed on to the Task constructor. None used at this time.
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
        self.makeSubtask("vignette")

    @classmethod
    def getInputDatasetTypes(cls, config):
        inputTypeDict = super().getInputDatasetTypes(config)

        # Delete entries from the dictionary of InputDatasetTypes that we know we don't
        # need because the configuration tells us we will not be bothering with the
        # correction that uses that IDT.
        if config.doBias is not True:
            inputTypeDict.pop("bias", None)
        if config.doLinearize is not True:
            inputTypeDict.pop("linearizer", None)
        if config.doCrosstalk is not True:
            inputTypeDict.pop("crosstalkSources", None)
        if config.doBrighterFatter is not True:
            inputTypeDict.pop("bfKernel", None)
        if config.doDefect is not True:
            inputTypeDict.pop("defects", None)
        if config.doDark is not True:
            inputTypeDict.pop("dark", None)
        if config.doFlat is not True:
            inputTypeDict.pop("flat", None)
        if config.doAttachTransmissionCurve is not True:
            inputTypeDict.pop("opticsTransmission", None)
            inputTypeDict.pop("filterTransmission", None)
            inputTypeDict.pop("sensorTransmission", None)
            inputTypeDict.pop("atmosphereTransmission", None)
        if config.doUseOpticsTransmission is not True:
            inputTypeDict.pop("opticsTransmission", None)
        if config.doUseFilterTransmission is not True:
            inputTypeDict.pop("filterTransmission", None)
        if config.doUseSensorTransmission is not True:
            inputTypeDict.pop("sensorTransmission", None)
        if config.doUseAtmosphereTransmission is not True:
            inputTypeDict.pop("atmosphereTransmission", None)

        return inputTypeDict

    @classmethod
    def getOutputDatasetTypes(cls, config):
        outputTypeDict = super().getOutputDatasetTypes(config)

        if config.qa.doThumbnailOss is not True:
            outputTypeDict.pop("outputOssThumbnail", None)
        if config.qa.doThumbnailFlattened is not True:
            outputTypeDict.pop("outputFlattenedThumbnail", None)
        if config.doWrite is not True:
            outputTypeDict.pop("outputExposure", None)

        return outputTypeDict

    def adaptArgsAndRun(self, inputData, inputDataIds, outputDataIds, butler):
        try:
            inputData['detectorNum'] = int(inputDataIds['ccdExposure']['detector'])
        except Exception as e:
            raise ValueError(f"Failure to find valid detectorNum value for Dataset {inputDataIds}: {e}")

        inputData['isGen3'] = True

        if self.config.doLinearize is True:
            if 'linearizer' not in inputData.keys():
                detector = inputData['camera'][inputData['detectorNum']]
                linearityName = detector.getAmpInfoCatalog()[0].getLinearityType()
                inputData['linearizer'] = linearize.getLinearityTypeByName(linearityName)()

        if inputData['defects'] is not None:
            # defects is loaded as a BaseCatalog with columns x0, y0, width, height.
            # masking expects a list of defects defined by their bounding box
            defectList = []

            for r in inputData['defects']:
                bbox = afwGeom.BoxI(afwGeom.PointI(r.get("x0"), r.get("y0")),
                                    afwGeom.ExtentI(r.get("width"), r.get("height")))
                defectList.append(Defect(bbox))

            inputData['defects'] = defectList

        # Broken: DM-17169
        # ci_hsc does not use crosstalkSources, as it's intra-CCD CT only.  This needs to be
        # fixed for non-HSC cameras in the future.
        # inputData['crosstalkSources'] = (self.crosstalk.prepCrosstalk(inputDataIds['ccdExposure'])
        #                        if self.config.doCrosstalk else None)

        # Broken: DM-17152
        # Fringes are not tested to be handled correctly by Gen3 butler.
        # inputData['fringes'] = (self.fringe.readFringes(inputDataIds['ccdExposure'],
        #                                                 assembler=self.assembleCcd
        #                                                 if self.config.doAssembleIsrExposures else None)
        #                         if self.config.doFringe and
        #                            self.fringe.checkFilter(inputData['ccdExposure'])
        #                         else pipeBase.Struct(fringes=None))

        return super().adaptArgsAndRun(inputData, inputDataIds, outputDataIds, butler)

    def makeDatasetType(self, dsConfig):
        return super().makeDatasetType(dsConfig)

    def readIsrData(self, dataRef, rawExposure):
        """!Retrieve necessary frames for instrument signature removal.

        Pre-fetching all required ISR data products limits the IO
        required by the ISR. Any conflict between the calibration data
        available and that needed for ISR is also detected prior to
        doing processing, allowing it to fail quickly.

        Parameters
        ----------
        dataRef : `daf.persistence.butlerSubset.ButlerDataRef`
            Butler reference of the detector data to be processed
        rawExposure : `afw.image.Exposure`
            The raw exposure that will later be corrected with the
            retrieved calibration data; should not be modified in this
            method.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with components (which may be `None`):
            - ``bias``: bias calibration frame (`afw.image.Exposure`)
            - ``linearizer``: functor for linearization (`ip.isr.linearize.LinearizeBase`)
            - ``crosstalkSources``: list of possible crosstalk sources (`list`)
            - ``dark``: dark calibration frame (`afw.image.Exposure`)
            - ``flat``: flat calibration frame (`afw.image.Exposure`)
            - ``bfKernel``: Brighter-Fatter kernel (`numpy.ndarray`)
            - ``defects``: list of defects (`list`)
            - ``fringes``: `lsst.pipe.base.Struct` with components:
                - ``fringes``: fringe calibration frame (`afw.image.Exposure`)
                - ``seed``: random seed derived from the ccdExposureId for random
                    number generator (`uint32`)
            - ``opticsTransmission``: `lsst.afw.image.TransmissionCurve`
                A ``TransmissionCurve`` that represents the throughput of the optics,
                to be evaluated in focal-plane coordinates.
            - ``filterTransmission`` : `lsst.afw.image.TransmissionCurve`
                A ``TransmissionCurve`` that represents the throughput of the filter
                itself, to be evaluated in focal-plane coordinates.
            - ``sensorTransmission`` : `lsst.afw.image.TransmissionCurve`
                A ``TransmissionCurve`` that represents the throughput of the sensor
                itself, to be evaluated in post-assembly trimmed detector coordinates.
            - ``atmosphereTransmission`` : `lsst.afw.image.TransmissionCurve`
                A ``TransmissionCurve`` that represents the throughput of the
                atmosphere, assumed to be spatially constant.

        """
        ccd = rawExposure.getDetector()
        rawExposure.mask.addMaskPlane("UNMASKEDNAN")  # needed to match pre DM-15862 processing.
        biasExposure = (self.getIsrExposure(dataRef, self.config.biasDataProductName)
                        if self.config.doBias else None)
        # immediate=True required for functors and linearizers are functors; see ticket DM-6515
        linearizer = (dataRef.get("linearizer", immediate=True)
                      if self.doLinearize(ccd) else None)
        crosstalkSources = (self.crosstalk.prepCrosstalk(dataRef)
                            if self.config.doCrosstalk else None)
        darkExposure = (self.getIsrExposure(dataRef, self.config.darkDataProductName)
                        if self.config.doDark else None)
        flatExposure = (self.getIsrExposure(dataRef, self.config.flatDataProductName)
                        if self.config.doFlat else None)
        brighterFatterKernel = (dataRef.get("bfKernel")
                                if self.config.doBrighterFatter else None)
        defectList = (dataRef.get("defects")
                      if self.config.doDefect else None)
        fringeStruct = (self.fringe.readFringes(dataRef, assembler=self.assembleCcd
                                                if self.config.doAssembleIsrExposures else None)
                        if self.config.doFringe and self.fringe.checkFilter(rawExposure)
                        else pipeBase.Struct(fringes=None))

        if self.config.doAttachTransmissionCurve:
            opticsTransmission = (dataRef.get("transmission_optics")
                                  if self.config.doUseOpticsTransmission else None)
            filterTransmission = (dataRef.get("transmission_filter")
                                  if self.config.doUseFilterTransmission else None)
            sensorTransmission = (dataRef.get("transmission_sensor")
                                  if self.config.doUseSensorTransmission else None)
            atmosphereTransmission = (dataRef.get("transmission_atmosphere")
                                      if self.config.doUseAtmosphereTransmission else None)
        else:
            opticsTransmission = None
            filterTransmission = None
            sensorTransmission = None
            atmosphereTransmission = None

        # Struct should include only kwargs to run()
        return pipeBase.Struct(bias=biasExposure,
                               linearizer=linearizer,
                               crosstalkSources=crosstalkSources,
                               dark=darkExposure,
                               flat=flatExposure,
                               bfKernel=brighterFatterKernel,
                               defects=defectList,
                               fringes=fringeStruct,
                               opticsTransmission=opticsTransmission,
                               filterTransmission=filterTransmission,
                               sensorTransmission=sensorTransmission,
                               atmosphereTransmission=atmosphereTransmission,
                               )

    @pipeBase.timeMethod
    def run(self, ccdExposure, camera=None, bias=None, linearizer=None, crosstalkSources=None,
            dark=None, flat=None, bfKernel=None, defects=None, fringes=None,
            opticsTransmission=None, filterTransmission=None,
            sensorTransmission=None, atmosphereTransmission=None,
            detectorNum=None, isGen3=False
            ):
        """!Perform instrument signature removal on an exposure.

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
            The camera geometry for this exposure.  Used to select the
            distortion model appropriate for this data.
        bias : `lsst.afw.image.Exposure`, optional
            Bias calibration frame.
        linearizer : `lsst.ip.isr.linearize.LinearizeBase`, optional
            Functor for linearization.
        crosstalkSources : `list`, optional
            List of possible crosstalk sources.
        dark : `lsst.afw.image.Exposure`, optional
            Dark calibration frame.
        flat : `lsst.afw.image.Exposure`, optional
            Flat calibration frame.
        bfKernel : `numpy.ndarray`, optional
            Brighter-fatter kernel.
        defects : `list`, optional
            List of defects.
        fringes : `lsst.pipe.base.Struct`, optional
            Struct containing the fringe correction data, with
            elements:
            - ``fringes``: fringe calibration frame (`afw.image.Exposure`)
            - ``seed``: random seed derived from the ccdExposureId for random
                number generator (`uint32`)
        opticsTransmission: `lsst.afw.image.TransmissionCurve`, optional
            A ``TransmissionCurve`` that represents the throughput of the optics,
            to be evaluated in focal-plane coordinates.
        filterTransmission : `lsst.afw.image.TransmissionCurve`
            A ``TransmissionCurve`` that represents the throughput of the filter
            itself, to be evaluated in focal-plane coordinates.
        sensorTransmission : `lsst.afw.image.TransmissionCurve`
            A ``TransmissionCurve`` that represents the throughput of the sensor
            itself, to be evaluated in post-assembly trimmed detector coordinates.
        atmosphereTransmission : `lsst.afw.image.TransmissionCurve`
            A ``TransmissionCurve`` that represents the throughput of the
            atmosphere, assumed to be spatially constant.
        detectorNum : `int`, optional
            The integer number for the detector to process.
        isGen3 : bool, optional
            Flag this call to run() as using the Gen3 butler environment.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with component:
            - ``exposure`` : `afw.image.Exposure`
                The fully ISR corrected exposure.
            - ``outputExposure`` : `afw.image.Exposure`
                An alias for `exposure`
            - ``ossThumb`` : `numpy.ndarray`
                Thumbnail image of the exposure after overscan subtraction.
            - ``flattenedThumb`` : `numpy.ndarray`
                Thumbnail image of the exposure after flat-field correction.

        Raises
        ------
        RuntimeError
            Raised if a configuration option is set to True, but the
            required calibration data has not been specified.

        Notes
        -----
        The current processed exposure can be viewed by setting the
        appropriate lsstDebug entries in the `debug.display`
        dictionary.  The names of these entries correspond to some of
        the IsrTaskConfig Boolean options, with the value denoting the
        frame to use.  The exposure is shown inside the matching
        option check and after the processing of that step has
        finished.  The steps with debug points are:

        doAssembleCcd
        doBias
        doCrosstalk
        doBrighterFatter
        doDark
        doFringe
        doStrayLight
        doFlat

        In addition, setting the "postISRCCD" entry displays the
        exposure after all ISR processing has finished.

        """

        if isGen3 is True:
            # Gen3 currently cannot automatically do configuration overrides.
            # DM-15257 looks to discuss this issue.

            self.config.doFringe = False

            # Configure input exposures;
            if detectorNum is None:
                raise RuntimeError("Must supply the detectorNum if running as Gen3")

            ccdExposure = self.ensureExposure(ccdExposure, camera, detectorNum)
            bias = self.ensureExposure(bias, camera, detectorNum)
            dark = self.ensureExposure(dark, camera, detectorNum)
            flat = self.ensureExposure(flat, camera, detectorNum)
        else:
            if isinstance(ccdExposure, ButlerDataRef):
                return self.runDataRef(ccdExposure)

        ccd = ccdExposure.getDetector()

        if not ccd:
            assert not self.config.doAssembleCcd, "You need a Detector to run assembleCcd"
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
        if fringes is None:
            fringes = pipeBase.Struct(fringes=None)
        if self.config.doFringe and not isinstance(fringes, pipeBase.Struct):
            raise RuntimeError("Must supply fringe exposure as a pipeBase.Struct.")
        if self.config.doFlat and flat is None:
            raise RuntimeError("Must supply a flat exposure if config.doFlat=True.")
        if self.config.doDefect and defects is None:
            raise RuntimeError("Must supply defects if config.doDefect=True.")
        if self.config.doAddDistortionModel and camera is None:
            raise RuntimeError("Must supply camera if config.doAddDistortionModel=True.")

        # Begin ISR processing.
        if self.config.doConvertIntToFloat:
            self.log.info("Converting exposure to floating point values")
            ccdExposure = self.convertIntToFloat(ccdExposure)

        # Amplifier level processing.
        overscans = []
        for amp in ccd:
            # if ccdExposure is one amp, check for coverage to prevent performing ops multiple times
            if ccdExposure.getBBox().contains(amp.getBBox()):
                # Check for fully masked bad amplifiers, and generate masks for SUSPECT and SATURATED values.
                badAmp = self.maskAmplifier(ccdExposure, amp, defects)

                if self.config.doOverscan and not badAmp:
                    # Overscan correction on amp-by-amp basis.
                    overscanResults = self.overscanCorrection(ccdExposure, amp)
                    self.log.debug("Corrected overscan for amplifier %s" % (amp.getName()))
                    if self.config.qa is not None and self.config.qa.saveStats is True:
                        if isinstance(overscanResults.overscanFit, float):
                            qaMedian = overscanResults.overscanFit
                            qaStdev = float("NaN")
                        else:
                            qaStats = afwMath.makeStatistics(overscanResults.overscanFit,
                                                             afwMath.MEDIAN | afwMath.STDEVCLIP)
                            qaMedian = qaStats.getValue(afwMath.MEDIAN)
                            qaStdev = qaStats.getValue(afwMath.STDEVCLIP)

                        self.metadata.set("ISR OSCAN {} MEDIAN".format(amp.getName()), qaMedian)
                        self.metadata.set("ISR OSCAN {} STDEV".format(amp.getName()), qaStdev)
                        self.log.debug("  Overscan stats for amplifer %s: %f +/- %f" %
                                       (amp.getName(), qaMedian, qaStdev))
                        ccdExposure.getMetadata().set('OVERSCAN', "Overscan corrected")
                else:
                    self.log.warn("Amplifier %s is bad." % (amp.getName()))
                    overscanResults = None

                overscans.append(overscanResults if overscanResults is not None else None)
            else:
                self.log.info("Skipped OSCAN")

        if self.config.doCrosstalk and self.config.doCrosstalkBeforeAssemble:
            self.log.info("Applying crosstalk correction.")
            self.crosstalk.run(ccdExposure, crosstalkSources=crosstalkSources)
            self.debugView(ccdExposure, "doCrosstalk")

        if self.config.doAssembleCcd:
            self.log.info("Assembling CCD from amplifiers")
            ccdExposure = self.assembleCcd.assembleCcd(ccdExposure)

            if self.config.expectWcs and not ccdExposure.getWcs():
                self.log.warn("No WCS found in input exposure")
            self.debugView(ccdExposure, "doAssembleCcd")

        ossThumb = None
        if self.config.qa.doThumbnailOss:
            ossThumb = isrQa.makeThumbnail(ccdExposure, isrQaConfig=self.config.qa)

        if self.config.doBias:
            self.log.info("Applying bias correction.")
            isrFunctions.biasCorrection(ccdExposure.getMaskedImage(), bias.getMaskedImage(),
                                        trimToFit=self.config.doTrimToMatchCalib)
            self.debugView(ccdExposure, "doBias")

        if self.config.doVariance:
            for amp, overscanResults in zip(ccd, overscans):
                if ccdExposure.getBBox().contains(amp.getBBox()):
                    self.log.debug("Constructing variance map for amplifer %s" % (amp.getName()))
                    ampExposure = ccdExposure.Factory(ccdExposure, amp.getBBox())
                    if overscanResults is not None:
                        self.updateVariance(ampExposure, amp,
                                            overscanImage=overscanResults.overscanImage)
                    else:
                        self.updateVariance(ampExposure, amp,
                                            overscanImage=None)
                    if self.config.qa is not None and self.config.qa.saveStats is True:
                        qaStats = afwMath.makeStatistics(ampExposure.getVariance(),
                                                         afwMath.MEDIAN | afwMath.STDEVCLIP)
                        self.metadata.set("ISR VARIANCE {} MEDIAN".format(amp.getName()),
                                          qaStats.getValue(afwMath.MEDIAN))
                        self.metadata.set("ISR VARIANCE {} STDEV".format(amp.getName()),
                                          qaStats.getValue(afwMath.STDEVCLIP))
                        self.log.debug("  Variance stats for amplifer %s: %f +/- %f" %
                                       (amp.getName(), qaStats.getValue(afwMath.MEDIAN),
                                        qaStats.getValue(afwMath.STDEVCLIP)))

        if self.doLinearize(ccd):
            self.log.info("Applying linearizer.")
            linearizer(image=ccdExposure.getMaskedImage().getImage(), detector=ccd, log=self.log)

        if self.config.doCrosstalk and not self.config.doCrosstalkBeforeAssemble:
            self.log.info("Applying crosstalk correction.")
            self.crosstalk.run(ccdExposure, crosstalkSources=crosstalkSources)
            self.debugView(ccdExposure, "doCrosstalk")

        if self.config.doWidenSaturationTrails:
            self.log.info("Widening saturation trails.")
            isrFunctions.widenSaturationTrails(ccdExposure.getMaskedImage().getMask())

        interpolationDone = False
        if self.config.doBrighterFatter:
            # We need to apply flats and darks before we can interpolate, and we
            # need to interpolate before we do B-F, but we do B-F without the
            # flats and darks applied so we can work in units of electrons or holes.
            # This context manager applies and then removes the darks and flats.
            with self.flatContext(ccdExposure, flat, dark):
                if self.config.doDefect:
                    self.maskAndInterpDefect(ccdExposure, defects)

                if self.config.doSaturationInterpolation:
                    self.saturationInterpolation(ccdExposure)

                self.maskAndInterpNan(ccdExposure)
                interpolationDone = True

            if self.config.brighterFatterLevel == 'DETECTOR':
                kernelElement = bfKernel
            else:
                # TODO: DM-15631 for implementing this
                raise NotImplementedError("per-amplifier brighter-fatter correction not yet implemented")
            self.log.info("Applying brighter fatter correction.")
            isrFunctions.brighterFatterCorrection(ccdExposure, kernelElement,
                                                  self.config.brighterFatterMaxIter,
                                                  self.config.brighterFatterThreshold,
                                                  self.config.brighterFatterApplyGain,
                                                  )
            self.debugView(ccdExposure, "doBrighterFatter")

        if self.config.doDark:
            self.log.info("Applying dark correction.")
            self.darkCorrection(ccdExposure, dark)
            self.debugView(ccdExposure, "doDark")

        if self.config.doFringe and not self.config.fringeAfterFlat:
            self.log.info("Applying fringe correction before flat.")
            self.fringe.run(ccdExposure, **fringes.getDict())
            self.debugView(ccdExposure, "doFringe")

        if self.config.doStrayLight:
            self.log.info("Applying stray light correction.")
            self.strayLight.run(ccdExposure)
            self.debugView(ccdExposure, "doStrayLight")

        if self.config.doFlat:
            self.log.info("Applying flat correction.")
            self.flatCorrection(ccdExposure, flat)
            self.debugView(ccdExposure, "doFlat")

        if self.config.doApplyGains:
            self.log.info("Applying gain correction instead of flat.")
            isrFunctions.applyGains(ccdExposure, self.config.normalizeGains)

        if self.config.doDefect and not interpolationDone:
            self.log.info("Masking and interpolating defects.")
            self.maskAndInterpDefect(ccdExposure, defects)

        if self.config.doSaturationInterpolation and not interpolationDone:
            self.log.info("Interpolating saturated pixels.")
            self.saturationInterpolation(ccdExposure)

        if self.config.doNanInterpAfterFlat or not interpolationDone:
            self.log.info("Masking and interpolating NAN value pixels.")
            self.maskAndInterpNan(ccdExposure)

        if self.config.doFringe and self.config.fringeAfterFlat:
            self.log.info("Applying fringe correction after flat.")
            self.fringe.run(ccdExposure, **fringes.getDict())

        if self.config.doSetBadRegions:
            badPixelCount, badPixelValue = isrFunctions.setBadRegions(ccdExposure)
            if badPixelCount > 0:
                self.log.info("Set %d BAD pixels to %f." % (badPixelCount, badPixelValue))

        flattenedThumb = None
        if self.config.qa.doThumbnailFlattened:
            flattenedThumb = isrQa.makeThumbnail(ccdExposure, isrQaConfig=self.config.qa)

        if self.config.doCameraSpecificMasking:
            self.log.info("Masking regions for camera specific reasons.")
            self.masking.run(ccdExposure)

        self.roughZeroPoint(ccdExposure)

        if self.config.doVignette:
            self.log.info("Constructing Vignette polygon.")
            self.vignettePolygon = self.vignette.run(ccdExposure)

            if self.config.vignette.doWriteVignettePolygon:
                self.setValidPolygonIntersect(ccdExposure, self.vignettePolygon)

        if self.config.doAttachTransmissionCurve:
            self.log.info("Adding transmission curves.")
            isrFunctions.attachTransmissionCurve(ccdExposure, opticsTransmission=opticsTransmission,
                                                 filterTransmission=filterTransmission,
                                                 sensorTransmission=sensorTransmission,
                                                 atmosphereTransmission=atmosphereTransmission)

        if self.config.doAddDistortionModel:
            self.log.info("Adding a distortion model to the WCS.")
            isrFunctions.addDistortionModel(exposure=ccdExposure, camera=camera)

        if self.config.doMeasureBackground:
            self.log.info("Measuring background level:")
            self.measureBackground(ccdExposure, self.config.qa)

            if self.config.qa is not None and self.config.qa.saveStats is True:
                for amp in ccd:
                    ampExposure = ccdExposure.Factory(ccdExposure, amp.getBBox())
                    qaStats = afwMath.makeStatistics(ampExposure.getImage(),
                                                     afwMath.MEDIAN | afwMath.STDEVCLIP)
                    self.metadata.set("ISR BACKGROUND {} MEDIAN".format(amp.getName()),
                                      qaStats.getValue(afwMath.MEDIAN))
                    self.metadata.set("ISR BACKGROUND {} STDEV".format(amp.getName()),
                                      qaStats.getValue(afwMath.STDEVCLIP))
                    self.log.debug("  Background stats for amplifer %s: %f +/- %f" %
                                   (amp.getName(), qaStats.getValue(afwMath.MEDIAN),
                                    qaStats.getValue(afwMath.STDEVCLIP)))

        self.debugView(ccdExposure, "postISRCCD")

        return pipeBase.Struct(
            exposure=ccdExposure,
            ossThumb=ossThumb,
            flattenedThumb=flattenedThumb,

            outputExposure=ccdExposure,
            outputOssThumbnail=ossThumb,
            outputFlattenedThumbnail=flattenedThumb,
        )

    @pipeBase.timeMethod
    def runDataRef(self, sensorRef):
        """Perform instrument signature removal on a ButlerDataRef of a Sensor.

        This method contains the `CmdLineTask` interface to the ISR
        processing.  All IO is handled here, freeing the `run()` method
        to manage only pixel-level calculations.  The steps performed
        are:
        - Read in necessary detrending/isr/calibration data.
        - Process raw exposure in `run()`.
        - Persist the ISR-corrected exposure as "postISRCCD" if
          config.doWrite=True.

        Parameters
        ----------
        sensorRef : `daf.persistence.butlerSubset.ButlerDataRef`
            DataRef of the detector data to be processed

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with component:
            - ``exposure`` : `afw.image.Exposure`
                The fully ISR corrected exposure.

        Raises
        ------
        RuntimeError
            Raised if a configuration option is set to True, but the
            required calibration data does not exist.

        """
        self.log.info("Performing ISR on sensor %s" % (sensorRef.dataId))

        ccdExposure = sensorRef.get(self.config.datasetType)

        camera = sensorRef.get("camera")
        if camera is None and self.config.doAddDistortionModel:
            raise RuntimeError("config.doAddDistortionModel is True "
                               "but could not get a camera from the butler")
        isrData = self.readIsrData(sensorRef, ccdExposure)

        result = self.run(ccdExposure, camera=camera, **isrData.getDict())

        if self.config.doWrite:
            sensorRef.put(result.exposure, "postISRCCD")
        if result.ossThumb is not None:
            isrQa.writeThumbnail(sensorRef, result.ossThumb, "ossThumb")
        if result.flattenedThumb is not None:
            isrQa.writeThumbnail(sensorRef, result.flattenedThumb, "flattenedThumb")

        return result

    def getIsrExposure(self, dataRef, datasetType, immediate=True):
        """!Retrieve a calibration dataset for removing instrument signature.

        Parameters
        ----------

        dataRef : `daf.persistence.butlerSubset.ButlerDataRef`
            DataRef of the detector data to find calibration datasets
            for.
        datasetType : `str`
            Type of dataset to retrieve (e.g. 'bias', 'flat', etc).
        immediate : `Bool`
            If True, disable butler proxies to enable error handling
            within this routine.

        Returns
        -------
        exposure : `lsst.afw.image.Exposure`
            Requested calibration frame.

        Raises
        ------
        RuntimeError
            Raised if no matching calibration frame can be found.
        """
        try:
            exp = dataRef.get(datasetType, immediate=immediate)
        except Exception as exc1:
            if not self.config.fallbackFilterName:
                raise RuntimeError("Unable to retrieve %s for %s: %s" % (datasetType, dataRef.dataId, exc1))
            try:
                exp = dataRef.get(datasetType, filter=self.config.fallbackFilterName, immediate=immediate)
            except Exception as exc2:
                raise RuntimeError("Unable to retrieve %s for %s, even with fallback filter %s: %s AND %s" %
                                   (datasetType, dataRef.dataId, self.config.fallbackFilterName, exc1, exc2))
            self.log.warn("Using fallback calibration from filter %s" % self.config.fallbackFilterName)

        if self.config.doAssembleIsrExposures:
            exp = self.assembleCcd.assembleCcd(exp)
        return exp

    def ensureExposure(self, inputExp, camera, detectorNum):
        """Ensure that the data returned by Butler is a fully constructed exposure.

        ISR requires exposure-level image data for historical reasons, so if we did
        not recieve that from Butler, construct it from what we have, modifying the
        input in place.

        Parameters
        ----------
        inputExp : `lsst.afw.image.Exposure`, `lsst.afw.image.DecoratedImageU`, or
                   `lsst.afw.image.ImageF`
            The input data structure obtained from Butler.
        camera : `lsst.afw.cameraGeom.camera`
            The camera associated with the image.  Used to find the appropriate
            detector.
        detectorNum : `int`
            The detector this exposure should match.

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
        else:
            raise TypeError(f"Input Exposure is not known type in isrTask.ensureExposure: {type(inputExp)}")

        if inputExp.getDetector() is None:
            inputExp.setDetector(camera[detectorNum])

        return inputExp

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
            return exposure
        if not hasattr(exposure, "convertF"):
            raise RuntimeError("Unable to convert exposure (%s) to float" % type(exposure))

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
        amp : `lsst.afw.table.AmpInfoCatalog`
            Catalog of parameters defining the amplifier on this
            exposure to mask.
        defects : `list`
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

        # Check if entire amp region is defined as a defect (need to use amp.getBBox() for correct
        # comparison with current defects definition.
        if defects is not None:
            badAmp = bool(sum([v.getBBox().contains(amp.getBBox()) for v in defects]))

        # In the case of a bad amp, we will set mask to "BAD" (here use amp.getRawBBox() for correct
        # association with pixels in current ccdExposure).
        if badAmp:
            dataView = afwImage.MaskedImageF(maskedImage, amp.getRawBBox(),
                                             afwImage.PARENT)
            maskView = dataView.getMask()
            maskView |= maskView.getPlaneBitMask("BAD")
            del maskView
            return badAmp

        # Mask remaining defects after assembleCcd() to allow for defects that cross amplifier boundaries.
        # Saturation and suspect pixels can be masked now, though.
        limits = dict()
        if self.config.doSaturation and not badAmp:
            limits.update({self.config.saturatedMaskName: amp.getSaturation()})
        if self.config.doSuspect and not badAmp:
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

        # Determine if we've fully masked this amplifier with SUSPECT and SAT pixels.
        maskView = afwImage.Mask(maskedImage.getMask(), amp.getRawDataBBox(),
                                 afwImage.PARENT)
        maskVal = maskView.getPlaneBitMask([self.config.saturatedMaskName,
                                            self.config.suspectMaskName])
        if numpy.all(maskView.getArray() & maskVal > 0):
            badAmp = True

        return badAmp

    def overscanCorrection(self, ccdExposure, amp):
        """Apply overscan correction in place.

        This method does initial pixel rejection of the overscan
        region.  The overscan can also be optionally segmented to
        allow for discontinuous overscan responses to be fit
        separately.  The actual overscan subtraction is performed by
        the `lsst.ip.isr.isrFunctions.overscanCorrection` function,
        which is called here after the amplifier is preprocessed.

        Parameters
        ----------
        ccdExposure : `lsst.afw.image.Exposure`
            Exposure to have overscan correction performed.
        amp : `lsst.afw.table.AmpInfoCatalog`
            The amplifier to consider while correcting the overscan.

        Returns
        -------
        overscanResults : `lsst.pipe.base.Struct`
            Result struct with components:
            - ``imageFit`` : scalar or `lsst.afw.image.Image`
                Value or fit subtracted from the amplifier image data.
            - ``overscanFit`` : scalar or `lsst.afw.image.Image`
                Value or fit subtracted from the overscan image data.
            - ``overscanImage`` : `lsst.afw.image.Image`
                Image of the overscan region with the overscan
                correction applied. This quantity is used to estimate
                the amplifier read noise empirically.

        Raises
        ------
        RuntimeError
            Raised if the ``amp`` does not contain raw pixel information.

        See Also
        --------
        lsst.ip.isr.isrFunctions.overscanCorrection
        """
        if not amp.getHasRawInfo():
            raise RuntimeError("This method must be executed on an amp with raw information.")

        if amp.getRawHorizontalOverscanBBox().isEmpty():
            self.log.info("ISR_OSCAN: No overscan region.  Not performing overscan correction.")
            return None

        statControl = afwMath.StatisticsControl()
        statControl.setAndMask(ccdExposure.mask.getPlaneBitMask("SAT"))

        # Determine the bounding boxes
        dataBBox = amp.getRawDataBBox()
        oscanBBox = amp.getRawHorizontalOverscanBBox()
        dx0 = 0
        dx1 = 0

        prescanBBox = amp.getRawPrescanBBox()
        if (oscanBBox.getBeginX() > prescanBBox.getBeginX()):  # amp is at the right
            dx0 += self.config.overscanNumLeadingColumnsToSkip
            dx1 -= self.config.overscanNumTrailingColumnsToSkip
        else:
            dx0 += self.config.overscanNumTrailingColumnsToSkip
            dx1 -= self.config.overscanNumLeadingColumnsToSkip

        # Determine if we need to work on subregions of the amplifier and overscan.
        imageBBoxes = []
        overscanBBoxes = []

        if ((self.config.overscanBiasJump and
             self.config.overscanBiasJumpLocation) and
            (ccdExposure.getMetadata().exists(self.config.overscanBiasJumpKeyword) and
             ccdExposure.getMetadata().getScalar(self.config.overscanBiasJumpKeyword) in
             self.config.overscanBiasJumpDevices)):
            if amp.getReadoutCorner() in (afwTable.LL, afwTable.LR):
                yLower = self.config.overscanBiasJumpLocation
                yUpper = dataBBox.getHeight() - yLower
            else:
                yUpper = self.config.overscanBiasJumpLocation
                yLower = dataBBox.getHeight() - yUpper

            imageBBoxes.append(afwGeom.Box2I(dataBBox.getBegin(),
                                             afwGeom.Extent2I(dataBBox.getWidth(), yLower)))
            overscanBBoxes.append(afwGeom.Box2I(oscanBBox.getBegin() +
                                                afwGeom.Extent2I(dx0, 0),
                                                afwGeom.Extent2I(oscanBBox.getWidth() - dx0 + dx1,
                                                                 yLower)))

            imageBBoxes.append(afwGeom.Box2I(dataBBox.getBegin() + afwGeom.Extent2I(0, yLower),
                                             afwGeom.Extent2I(dataBBox.getWidth(), yUpper)))
            overscanBBoxes.append(afwGeom.Box2I(oscanBBox.getBegin() + afwGeom.Extent2I(dx0, yLower),
                                                afwGeom.Extent2I(oscanBBox.getWidth() - dx0 + dx1,
                                                                 yUpper)))
        else:
            imageBBoxes.append(afwGeom.Box2I(dataBBox.getBegin(),
                                             afwGeom.Extent2I(dataBBox.getWidth(), dataBBox.getHeight())))
            overscanBBoxes.append(afwGeom.Box2I(oscanBBox.getBegin() + afwGeom.Extent2I(dx0, 0),
                                                afwGeom.Extent2I(oscanBBox.getWidth() - dx0 + dx1,
                                                                 oscanBBox.getHeight())))

        # Perform overscan correction on subregions, ensuring saturated pixels are masked.
        for imageBBox, overscanBBox in zip(imageBBoxes, overscanBBoxes):
            ampImage = ccdExposure.maskedImage[imageBBox]
            overscanImage = ccdExposure.maskedImage[overscanBBox]

            overscanArray = overscanImage.image.array
            median = numpy.ma.median(numpy.ma.masked_where(overscanImage.mask.array, overscanArray))
            bad = numpy.where(numpy.abs(overscanArray - median) > self.config.overscanMaxDev)
            overscanImage.mask.array[bad] = overscanImage.mask.getPlaneBitMask("SAT")

            statControl = afwMath.StatisticsControl()
            statControl.setAndMask(ccdExposure.mask.getPlaneBitMask("SAT"))

            overscanResults = isrFunctions.overscanCorrection(ampMaskedImage=ampImage,
                                                              overscanImage=overscanImage,
                                                              fitType=self.config.overscanFitType,
                                                              order=self.config.overscanOrder,
                                                              collapseRej=self.config.overscanNumSigmaClip,
                                                              statControl=statControl,
                                                              overscanIsInt=self.config.overscanIsInt
                                                              )

            # Measure average overscan levels and record them in the metadata
            levelStat = afwMath.MEDIAN
            sigmaStat = afwMath.STDEVCLIP

            sctrl = afwMath.StatisticsControl(self.config.qa.flatness.clipSigma,
                                              self.config.qa.flatness.nIter)
            metadata = ccdExposure.getMetadata()
            ampNum = amp.getName()
            if self.config.overscanFitType in ("MEDIAN", "MEAN", "MEANCLIP"):
                metadata.set("ISR_OSCAN_LEVEL%s" % ampNum, overscanResults.overscanFit)
                metadata.set("ISR_OSCAN_SIGMA%s" % ampNum, 0.0)
            else:
                stats = afwMath.makeStatistics(overscanResults.overscanFit, levelStat | sigmaStat, sctrl)
                metadata.set("ISR_OSCAN_LEVEL%s" % ampNum, stats.getValue(levelStat))
                metadata.set("ISR_OSCAN_SIGMA%s" % ampNum, stats.getValue(sigmaStat))

        return overscanResults

    def updateVariance(self, ampExposure, amp, overscanImage=None):
        """Set the variance plane using the amplifier gain and read noise

        The read noise is calculated from the ``overscanImage`` if the
        ``doEmpiricalReadNoise`` option is set in the configuration; otherwise
        the value from the amplifier data is used.

        Parameters
        ----------
        ampExposure : `lsst.afw.image.Exposure`
            Exposure to process.
        amp : `lsst.afw.table.AmpInfoRecord` or `FakeAmp`
            Amplifier detector data.
        overscanImage : `lsst.afw.image.MaskedImage`, optional.
            Image of overscan, required only for empirical read noise.

        See also
        --------
        lsst.ip.isr.isrFunctions.updateVariance
        """
        maskPlanes = [self.config.saturatedMaskName, self.config.suspectMaskName]
        gain = amp.getGain()

        if math.isnan(gain):
            gain = 1.0
            self.log.warn("Gain set to NAN!  Updating to 1.0 to generate Poisson variance.")
        elif gain <= 0:
            patchedGain = 1.0
            self.log.warn("Gain for amp %s == %g <= 0; setting to %f" %
                          (amp.getName(), gain, patchedGain))
            gain = patchedGain

        if self.config.doEmpiricalReadNoise and overscanImage is None:
            self.log.info("Overscan is none for EmpiricalReadNoise")

        if self.config.doEmpiricalReadNoise and overscanImage is not None:
            stats = afwMath.StatisticsControl()
            stats.setAndMask(overscanImage.mask.getPlaneBitMask(maskPlanes))
            readNoise = afwMath.makeStatistics(overscanImage, afwMath.STDEVCLIP, stats).getValue()
            self.log.info("Calculated empirical read noise for amp %s: %f", amp.getName(), readNoise)
        else:
            readNoise = amp.getReadNoise()

        isrFunctions.updateVariance(
            maskedImage=ampExposure.getMaskedImage(),
            gain=gain,
            readNoise=readNoise,
        )

    def darkCorrection(self, exposure, darkExposure, invert=False):
        """!Apply dark correction in place.

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
            raise RuntimeError("Exposure darktime is NAN")
        if darkExposure.getInfo().getVisitInfo() is not None:
            darkScale = darkExposure.getInfo().getVisitInfo().getDarkTime()
        else:
            # DM-17444: darkExposure.getInfo.getVisitInfo() is None
            #           so getDarkTime() does not exist.
            darkScale = 1.0

        if math.isnan(darkScale):
            raise RuntimeError("Dark calib darktime is NAN")
        isrFunctions.darkCorrection(
            maskedImage=exposure.getMaskedImage(),
            darkMaskedImage=darkExposure.getMaskedImage(),
            expScale=expScale,
            darkScale=darkScale,
            invert=invert,
            trimToFit=self.config.doTrimToMatchCalib
        )

    def doLinearize(self, detector):
        """!Check if linearization is needed for the detector cameraGeom.

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
            detector.getAmpInfoCatalog()[0].getLinearityType() != NullLinearityType

    def flatCorrection(self, exposure, flatExposure, invert=False):
        """!Apply flat correction in place.

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
        """!Detect saturated pixels and mask them using mask plane config.saturatedMaskName, in place.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to process.  Only the amplifier DataSec is processed.
        amp : `lsst.afw.table.AmpInfoCatalog`
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

    def saturationInterpolation(self, ccdExposure):
        """!Interpolate over saturated pixels, in place.

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
            maskedImage=ccdExposure.getMaskedImage(),
            fwhm=self.config.fwhm,
            growFootprints=self.config.growSaturationFootprintSize,
            maskName=self.config.saturatedMaskName,
        )

    def suspectDetection(self, exposure, amp):
        """!Detect suspect pixels and mask them using mask plane config.suspectMaskName, in place.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to process.  Only the amplifier DataSec is processed.
        amp : `lsst.afw.table.AmpInfoCatalog`
            Amplifier detector data.

        See Also
        --------
        lsst.ip.isr.isrFunctions.makeThresholdMask

        Notes
        -----
        Suspect pixels are pixels whose value is greater than amp.getSuspectLevel().
        This is intended to indicate pixels that may be affected by unknown systematics;
        for example if non-linearity corrections above a certain level are unstable
        then that would be a useful value for suspectLevel. A value of `nan` indicates
        that no such level exists and no pixels are to be masked as suspicious.
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

    def maskAndInterpDefect(self, ccdExposure, defectBaseList):
        """!Mask defects using mask plane "BAD" and interpolate over them, in place.

        Parameters
        ----------
        ccdExposure : `lsst.afw.image.Exposure`
            Exposure to process.
        defectBaseList : `List`
            List of defects to mask and interpolate.

        Notes
        -----
        Call this after CCD assembly, since defects may cross amplifier boundaries.
        """
        maskedImage = ccdExposure.getMaskedImage()
        defectList = []
        for d in defectBaseList:
            bbox = d.getBBox()
            nd = measAlg.Defect(bbox)
            defectList.append(nd)
        isrFunctions.maskPixelsFromDefectList(maskedImage, defectList, maskName='BAD')
        isrFunctions.interpolateDefectList(
            maskedImage=maskedImage,
            defectList=defectList,
            fwhm=self.config.fwhm,
        )

        if self.config.numEdgeSuspect > 0:
            goodBBox = maskedImage.getBBox()
            # This makes a bbox numEdgeSuspect pixels smaller than the image on each side
            goodBBox.grow(-self.config.numEdgeSuspect)
            # Mask pixels outside goodBBox as SUSPECT
            SourceDetectionTask.setEdgeBits(
                maskedImage,
                goodBBox,
                maskedImage.getMask().getPlaneBitMask("SUSPECT")
            )

    def maskAndInterpNan(self, exposure):
        """!Mask NaNs using mask plane "UNMASKEDNAN" and interpolate over them, in place.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to process.

        Notes
        -----
        We mask and interpolate over all NaNs, including those
        that are masked with other bits (because those may or may
        not be interpolated over later, and we want to remove all
        NaNs).  Despite this behaviour, the "UNMASKEDNAN" mask plane
        is used to preserve the historical name.
        """
        maskedImage = exposure.getMaskedImage()

        # Find and mask NaNs
        maskedImage.getMask().addMaskPlane("UNMASKEDNAN")
        maskVal = maskedImage.getMask().getPlaneBitMask("UNMASKEDNAN")
        numNans = maskNans(maskedImage, maskVal)
        self.metadata.set("NUMNANS", numNans)

        # Interpolate over these previously-unmasked NaNs
        if numNans > 0:
            self.log.warn("There were %i unmasked NaNs", numNans)
            nanDefectList = isrFunctions.getDefectListFromMask(
                maskedImage=maskedImage,
                maskName='UNMASKEDNAN',
            )
            isrFunctions.interpolateDefectList(
                maskedImage=exposure.getMaskedImage(),
                defectList=nanDefectList,
                fwhm=self.config.fwhm,
            )

    def measureBackground(self, exposure, IsrQaConfig=None):
        """Measure the image background in subgrids, for quality control purposes.

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
            self.log.info("Flattened sky level: %f +/- %f" % (skyLevel, skySigma))
            metadata = exposure.getMetadata()
            metadata.set('SKYLEVEL', skyLevel)
            metadata.set('SKYSIGMA', skySigma)

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

                    bbox = afwGeom.Box2I(afwGeom.Point2I(xLLC, yLLC), afwGeom.Point2I(xURC, yURC))
                    miMesh = maskedImage.Factory(exposure.getMaskedImage(), bbox, afwImage.LOCAL)

                    skyLevels[i, j] = afwMath.makeStatistics(miMesh, stat, statsControl).getValue()

            good = numpy.where(numpy.isfinite(skyLevels))
            skyMedian = numpy.median(skyLevels[good])
            flatness = (skyLevels[good] - skyMedian) / skyMedian
            flatness_rms = numpy.std(flatness)
            flatness_pp = flatness.max() - flatness.min() if len(flatness) > 0 else numpy.nan

            self.log.info("Measuring sky levels in %dx%d grids: %f" % (nX, nY, skyMedian))
            self.log.info("Sky flatness in %dx%d grids - pp: %f rms: %f" %
                          (nX, nY, flatness_pp, flatness_rms))

            metadata.set('FLATNESS_PP', float(flatness_pp))
            metadata.set('FLATNESS_RMS', float(flatness_rms))
            metadata.set('FLATNESS_NGRIDS', '%dx%d' % (nX, nY))
            metadata.set('FLATNESS_MESHX', IsrQaConfig.flatness.meshX)
            metadata.set('FLATNESS_MESHY', IsrQaConfig.flatness.meshY)

    def roughZeroPoint(self, exposure):
        """Set an approximate magnitude zero point for the exposure.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to process.
        """
        filterName = afwImage.Filter(exposure.getFilter().getId()).getName()  # Canonical name for filter
        if filterName in self.config.fluxMag0T1:
            fluxMag0 = self.config.fluxMag0T1[filterName]
        else:
            self.log.warn("No rough magnitude zero point set for filter %s" % filterName)
            fluxMag0 = self.config.defaultFluxMag0T1

        expTime = exposure.getInfo().getVisitInfo().getExposureTime()
        if not expTime > 0:  # handle NaN as well as <= 0
            self.log.warn("Non-positive exposure time; skipping rough zero point")
            return

        self.log.info("Setting rough magnitude zero point: %f" % (2.5*math.log10(fluxMag0*expTime),))
        exposure.getCalib().setFluxMag0(fluxMag0*expTime)

    def setValidPolygonIntersect(self, ccdExposure, fpPolygon):
        """!Set the valid polygon as the intersection of fpPolygon and the ccd corners.

        Parameters
        ----------
        ccdExposure : `lsst.afw.image.Exposure`
            Exposure to process.
        fpPolygon : `lsst.afw.geom.Polygon`
            Polygon in focal plane coordinates.
        """
        # Get ccd corners in focal plane coordinates
        ccd = ccdExposure.getDetector()
        fpCorners = ccd.getCorners(FOCAL_PLANE)
        ccdPolygon = Polygon(fpCorners)

        # Get intersection of ccd corners with fpPolygon
        intersect = ccdPolygon.intersectionSingle(fpPolygon)

        # Transform back to pixel positions and build new polygon
        ccdPoints = ccd.transform(intersect, FOCAL_PLANE, PIXELS)
        validPolygon = Polygon(ccdPoints)
        ccdExposure.getInfo().setValidPolygon(validPolygon)

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
        self._RawHorizontalOverscanBBox = afwGeom.Box2I()
        self._gain = config.gain
        self._readNoise = config.readNoise
        self._saturation = config.saturation

    def getBBox(self):
        return self._bbox

    def getRawBBox(self):
        return self._bbox

    def getHasRawInfo(self):
        return True                     # but see getRawHorizontalOverscanBBox()

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


class RunIsrConfig(pexConfig.Config):
    isr = pexConfig.ConfigurableField(target=IsrTask, doc="Instrument signature removal")


class RunIsrTask(pipeBase.CmdLineTask):
    """Task to wrap the default IsrTask to allow it to be retargeted.

    The standard IsrTask can be called directly from a command line
    program, but doing so removes the ability of the task to be
    retargeted.  As most cameras override some set of the IsrTask
    methods, this would remove those data-specific methods in the
    output post-ISR images.  This wrapping class fixes the issue,
    allowing identical post-ISR images to be generated by both the
    processCcd and isrTask code.
    """
    ConfigClass = RunIsrConfig
    _DefaultName = "runIsr"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.makeSubtask("isr")

    def runDataRef(self, dataRef):
        """
        Parameters
        ----------
        dataRef : `lsst.daf.persistence.ButlerDataRef`
            data reference of the detector data to be processed

        Returns
        -------
        result : `pipeBase.Struct`
            Result struct with component:

            - exposure : `lsst.afw.image.Exposure`
                Post-ISR processed exposure.
        """
        return self.isr.runDataRef(dataRef)
