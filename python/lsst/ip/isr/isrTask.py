#
# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import math
import numpy

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.meas.algorithms as measAlg
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.math as afwMath
from lsst.daf.persistence import ButlerDataRef
from lsstDebug import getDebugFrame
from lsst.afw.display import getDisplay
from . import isrFunctions
from .assembleCcdTask import AssembleCcdTask
from .fringe import FringeTask
from lsst.afw.geom import Polygon
from lsst.afw.geom.wcsUtils import makeDistortedTanWcs
from lsst.afw.cameraGeom import PIXELS, FOCAL_PLANE, FIELD_ANGLE, NullLinearityType
from contextlib import contextmanager
from .isr import maskNans
from .crosstalk import CrosstalkTask

__all__ = ('IsrTaskConfig', 'IsrTask', 'FakeAmp')

class IsrTaskConfig(pexConfig.Config):
    doBias = pexConfig.Field(
        dtype=bool,
        doc="Apply bias frame correction?",
        default=True,
    )
    doDark = pexConfig.Field(
        dtype=bool,
        doc="Apply dark frame correction?",
        default=True,
    )
    doFlat = pexConfig.Field(
        dtype=bool,
        doc="Apply flat field correction?",
        default=True,
    )
    doFringe = pexConfig.Field(
        dtype=bool,
        doc="Apply fringe correction?",
        default=True,
    )
    doDefect = pexConfig.Field(
        dtype=bool,
        doc="Apply correction for CCD defects, e.g. hot pixels?",
        default=True,
    )
    doAddDistortionModel = pexConfig.Field(
        dtype=bool,
        doc="Apply a distortion model based on camera geometry to the WCS?",
        default=True,
    )
    doWrite = pexConfig.Field(
        dtype=bool,
        doc="Persist postISRCCD?",
        default=True,
    )
    biasDataProductName = pexConfig.Field(
        dtype=str,
        doc="Name of the bias data product",
        default="bias",
    )
    darkDataProductName = pexConfig.Field(
        dtype=str,
        doc="Name of the dark data product",
        default="dark",
    )
    flatDataProductName = pexConfig.Field(
        dtype=str,
        doc="Name of the flat data product",
        default="flat",
    )
    assembleCcd = pexConfig.ConfigurableField(
        target=AssembleCcdTask,
        doc="CCD assembly task",
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
    saturation = pexConfig.Field(
        dtype=float,
        doc="The saturation level to use if no Detector is present in the Exposure (ignored if NaN)",
        default=float("NaN"),
    )
    fringeAfterFlat = pexConfig.Field(
        dtype=bool,
        doc="Do fringe subtraction after flat-fielding?",
        default=True,
    )
    fringe = pexConfig.ConfigurableField(
        target=FringeTask,
        doc="Fringe subtraction task",
    )
    fwhm = pexConfig.Field(
        dtype=float,
        doc="FWHM of PSF (arcsec)",
        default=1.0,
    )
    saturatedMaskName = pexConfig.Field(
        dtype=str,
        doc="Name of mask plane to use in saturation detection and interpolation",
        default="SAT",
    )
    suspectMaskName = pexConfig.Field(
        dtype=str,
        doc="Name of mask plane to use for suspect pixels",
        default="SUSPECT",
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
            "MEDIAN": "Correct using the median of the overscan region",
        },
    )
    overscanOrder = pexConfig.Field(
        dtype=int,
        doc=("Order of polynomial or to fit if overscan fit type is a polynomial, " +
             "or number of spline knots if overscan fit type is a spline."),
        default=1,
    )
    overscanRej = pexConfig.Field(
        dtype=float,
        doc="Rejection threshold (sigma) for collapsing overscan before fit",
        default=3.0,
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
    growSaturationFootprintSize = pexConfig.Field(
        dtype=int,
        doc="Number of pixels by which to grow the saturation footprints",
        default=1,
    )
    doSaturationInterpolation = pexConfig.Field(
        dtype=bool,
        doc="Perform interpolation over pixels masked as saturated?",
        default=True,
    )
    doNanInterpAfterFlat = pexConfig.Field(
        dtype=bool,
        doc=("If True, ensure we interpolate NaNs after flat-fielding, even if we "
             "also have to interpolate them before flat-fielding."),
        default=False,
    )
    fluxMag0T1 = pexConfig.Field(
        dtype=float,
        doc="The approximate flux of a zero-magnitude object in a one-second exposure",
        default=1e10,
    )
    keysToRemoveFromAssembledCcd = pexConfig.ListField(
        dtype=str,
        doc="fields to remove from the metadata of the assembled ccd.",
        default=[],
    )
    doAssembleIsrExposures = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Assemble amp-level calibration exposures into ccd-level exposure?"
    )
    doAssembleCcd = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Assemble amp-level exposures into a ccd-level exposure?"
    )
    expectWcs = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Expect input science images to have a WCS (set False for e.g. spectrographs)"
    )
    doLinearize = pexConfig.Field(
        dtype=bool,
        doc="Correct for nonlinearity of the detector's response?",
        default=True,
    )
    doCrosstalk = pexConfig.Field(
        dtype=bool,
        doc="Apply intra-CCD crosstalk correction?",
        default=False,
    )
    crosstalk = pexConfig.ConfigurableField(
        target=CrosstalkTask,
        doc="Intra-CCD crosstalk correction",
    )
    doBrighterFatter = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Apply the brighter fatter correction"
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
    datasetType = pexConfig.Field(
        dtype=str,
        doc="Dataset type for input data; users will typically leave this alone, "
        "but camera-specific ISR tasks will override it",
        default="raw",
    )
    fallbackFilterName = pexConfig.Field(dtype=str,
                                         doc="Fallback default filter name for calibrations", optional=True)
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
    doEmpiricalReadNoise = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Calculate empirical read noise instead of value from AmpInfo data?"
    )

## @addtogroup LSST_task_documentation
## @{
## @page IsrTask
## @ref IsrTask_ "IsrTask"
## @copybrief IsrTask
## @}


class IsrTask(pipeBase.CmdLineTask):
    """The process for correcting imaging data is very similar from camera to camera.
    This task provides a vanilla implementation of doing these corrections, including
    the ability to turn certain corrections off if they are not needed.
    The inputs to the primary method, run, are a raw exposure to be corrected and the
    calibration data products. The raw input is a single chip sized mosaic of all amps
    including overscans and other non-science pixels.
    The method runDataRef() is intended for use by a lsst.pipe.base.cmdLineTask.CmdLineTask
    and takes as input only a daf.persistence.butlerSubset.ButlerDataRef.
    This task may not meet all needs and it is expected that it will be subclassed for
    specific applications.

    Parameters
    ----------
    args :
        a list of positional arguments passed on to the Task constructor
    kwargs :
        a dictionary of keyword arguments passed on to the Task constructor
        Call the ``lsst.pipe.base.task.Task.__init__`` method
        Then setup the assembly and fringe correction subtasks

    Notes
    -----
    The available variables in IsrTask are:

    display :
        A dictionary containing debug point names as keys with frame number as value. Valid keys are:

    postISRCCD :
        display exposure after ISR has been applied

    Examples
    --------

    For example, put something like

    .. code-block :: none

        import lsstDebug
        def DebugInfo(name):
            di = lsstDebug.getInfo(name)        # N.b. lsstDebug.Info(name) would call us recursively
            if name == "lsst.ip.isrFunctions.isrTask":
                di.display = {'postISRCCD':2}
            return di
        lsstDebug.Info = DebugInfo

    into your debug.py file and run the commandline task with the @c --debug flag.

    """
    ConfigClass = IsrTaskConfig
    _DefaultName = "isr"

    def __init__(self, *args, **kwargs):
        '''Constructor for IsrTask
        '''
        pipeBase.Task.__init__(self, *args, **kwargs)
        self.makeSubtask("assembleCcd")
        self.makeSubtask("fringe")
        self.makeSubtask("crosstalk")

    def readIsrData(self, dataRef, rawExposure):
        """Retrieve necessary frames for instrument signature removal

        Parameters
        ----------
        dataRef : `daf.persistence.butlerSubset.ButlerDataRef`
            a daf.persistence.butlerSubset.ButlerDataRef
            of the detector data to be processed
        rawExposure :
            a reference raw exposure that will later be
            corrected with the retrieved calibration data;
            should not be modified in this method.

        Returns
        -------
        Struct : `pipeBase.Struct`
            a pipeBase.Struct with fields containing kwargs expected by run():

            - ``bias`` : exposure of bias frame
            - ``dark`` : exposure of dark frame
            - ``flat`` : exposure of flat field
            - ``defects`` : list of detects
            - ``fringeStruct`` : a pipeBase.Struct with field fringes containing
                exposure of fringe frame or list of fringe exposure

        """
        ccd = rawExposure.getDetector()

        biasExposure = self.getIsrExposure(dataRef, self.config.biasDataProductName) \
            if self.config.doBias else None
        # immediate=True required for functors and linearizers are functors; see ticket DM-6515
        linearizer = dataRef.get("linearizer", immediate=True) if self.doLinearize(ccd) else None
        darkExposure = self.getIsrExposure(dataRef, self.config.darkDataProductName) \
            if self.config.doDark else None
        flatExposure = self.getIsrExposure(dataRef, self.config.flatDataProductName) \
            if self.config.doFlat else None
        brighterFatterKernel = dataRef.get("brighterFatterKernel") if self.config.doBrighterFatter else None
        defectList = dataRef.get("defects") if self.config.doDefect else None

        if self.config.doCrosstalk:
            crosstalkSources = self.crosstalk.prepCrosstalk(dataRef)
        else:
            crosstalkSources = None

        if self.config.doFringe and self.fringe.checkFilter(rawExposure):
            fringeStruct = self.fringe.readFringes(dataRef, assembler=self.assembleCcd
                                                   if self.config.doAssembleIsrExposures else None)
        else:
            fringeStruct = pipeBase.Struct(fringes=None)

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
                               dark=darkExposure,
                               flat=flatExposure,
                               defects=defectList,
                               fringes=fringeStruct,
                               bfKernel=brighterFatterKernel,
                               opticsTransmission=opticsTransmission,
                               filterTransmission=filterTransmission,
                               sensorTransmission=sensorTransmission,
                               atmosphereTransmission=atmosphereTransmission,
                               crosstalkSources=crosstalkSources,
                               )

    @pipeBase.timeMethod
    def run(self, ccdExposure, bias=None, linearizer=None, dark=None, flat=None, defects=None,
            fringes=None, bfKernel=None, camera=None,
            opticsTransmission=None, filterTransmission=None,
            sensorTransmission=None, atmosphereTransmission=None,
            crosstalkSources=None):
        """Perform instrument signature removal on an exposure

        Steps include:

        - Detect saturation, apply overscan correction, bias, dark and flat
        - Perform CCD assembly
        - Interpolate over defects, saturated pixels and all NaNs

        Parameters
        ----------
        ccdExposure : `lsst.afw.image.exposure`
            lsst.afw.image.exposure of detector data
        bias :
            exposure of bias frame
        linearizer :
            linearizing functor; a subclass of lsst.ip.isrFunctions.LinearizeBase
        dark :
            exposure of dark frame
        flat :
            exposure of flatfield
        defects :
            list of detects
        fringes :
            a pipeBase.Struct with field fringes containing
            exposure of fringe frame or list of fringe exposure
        bfKernel :
            kernel for brighter-fatter correction
        camera : `lsst.afw.cameraGeom.Camera`
            camera geometry, used by addDistortionModel
        opticsTransmission :
            a TransmissionCurve for the optics
        filterTransmission :
            a TransmissionCurve for the filter
        sensorTransmission :
            a TransmissionCurve for the sensor
        atmosphereTransmission :
            a TransmissionCurve for the atmosphere
        crosstalkSources :
            a defaultdict used for DECam inter-CCD crosstalk

        Returns
        -------
        Struct :  ``pipeBase.Struct``
            with field:

            - ``exposure`` :

        """
        # parseAndRun expects to be able to call run() with a dataRef; see DM-6640
        if isinstance(ccdExposure, ButlerDataRef):
            return self.runDataRef(ccdExposure)

        ccd = ccdExposure.getDetector()

        # Validate Input
        if self.config.doBias and bias is None:
            raise RuntimeError("Must supply a bias exposure if config.doBias True")
        if self.doLinearize(ccd) and linearizer is None:
            raise RuntimeError("Must supply a linearizer if config.doBias True")
        if self.config.doDark and dark is None:
            raise RuntimeError("Must supply a dark exposure if config.doDark True")
        if self.config.doFlat and flat is None:
            raise RuntimeError("Must supply a flat exposure if config.doFlat True")
        if self.config.doBrighterFatter and bfKernel is None:
            raise RuntimeError("Must supply a kernel if config.doBrighterFatter True")
        if fringes is None:
            fringes = pipeBase.Struct(fringes=None)
        if self.config.doFringe and not isinstance(fringes, pipeBase.Struct):
            raise RuntimeError("Must supply fringe exposure as a pipeBase.Struct")
        if self.config.doDefect and defects is None:
            raise RuntimeError("Must supply defects if config.doDefect True")
        if self.config.doAddDistortionModel and camera is None:
            raise RuntimeError("Must supply camera if config.doAddDistortionModel True")

        ccdExposure = self.convertIntToFloat(ccdExposure)

        if not ccd:
            assert not self.config.doAssembleCcd, "You need a Detector to run assembleCcd"
            ccd = [FakeAmp(ccdExposure, self.config)]

        overscans = []
        for amp in ccd:
            # if ccdExposure is one amp, check for coverage to prevent performing ops multiple times
            if ccdExposure.getBBox().contains(amp.getBBox()):
                self.saturationDetection(ccdExposure, amp)
                self.suspectDetection(ccdExposure, amp)
                overscanResults = self.overscanCorrection(ccdExposure, amp)
                overscans.append(overscanResults.overscanImage if overscanResults is not None else None)
            else:
                overscans.append(None)

        if self.config.doCrosstalk:
            self.crosstalk.run(ccdExposure, crosstalkSources)

        if self.config.doAssembleCcd:
            ccdExposure = self.assembleCcd.assembleCcd(ccdExposure)
            if self.config.expectWcs and not ccdExposure.getWcs():
                self.log.warn("No WCS found in input exposure")

        if self.config.doBias:
            self.biasCorrection(ccdExposure, bias)

        if self.doLinearize(ccd):
            linearizer(image=ccdExposure.getMaskedImage().getImage(), detector=ccd, log=self.log)

        assert len(ccd) == len(overscans)
        for amp, overscanImage in zip(ccd, overscans):
            # if ccdExposure is one amp, check for coverage to prevent performing ops multiple times
            if ccdExposure.getBBox().contains(amp.getBBox()):
                ampExposure = ccdExposure.Factory(ccdExposure, amp.getBBox())
                self.updateVariance(ampExposure, amp, overscanImage)

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

            self.brighterFatterCorrection(ccdExposure, bfKernel,
                                          self.config.brighterFatterMaxIter,
                                          self.config.brighterFatterThreshold,
                                          self.config.brighterFatterApplyGain,
                                          )

        if self.config.doDark:
            self.darkCorrection(ccdExposure, dark)

        if self.config.doFringe and not self.config.fringeAfterFlat:
            self.fringe.run(ccdExposure, **fringes.getDict())

        if self.config.doFlat:
            self.flatCorrection(ccdExposure, flat)

        if not interpolationDone:
            if self.config.doDefect:
                self.maskAndInterpDefect(ccdExposure, defects)
            if self.config.doSaturationInterpolation:
                self.saturationInterpolation(ccdExposure)
        if not interpolationDone or self.config.doNanInterpAfterFlat:
            self.maskAndInterpNan(ccdExposure)

        if self.config.doFringe and self.config.fringeAfterFlat:
            self.fringe.run(ccdExposure, **fringes.getDict())

        exposureTime = ccdExposure.getInfo().getVisitInfo().getExposureTime()
        ccdExposure.getCalib().setFluxMag0(self.config.fluxMag0T1*exposureTime)

        if self.config.doAddDistortionModel:
            self.addDistortionModel(exposure=ccdExposure, camera=camera)

        if self.config.doAttachTransmissionCurve:
            self.attachTransmissionCurve(ccdExposure, opticsTransmission=opticsTransmission,
                                         filterTransmission=filterTransmission,
                                         sensorTransmission=sensorTransmission,
                                         atmosphereTransmission=atmosphereTransmission)

        frame = getDebugFrame(self._display, "postISRCCD")
        if frame:
            getDisplay(frame).mtv(ccdExposure)

        return pipeBase.Struct(
            exposure=ccdExposure,
        )

    @pipeBase.timeMethod
    def runDataRef(self, sensorRef):
        """Perform instrument signature removal on a ButlerDataRef of a Sensor

        - Read in necessary detrending/isr/calibration data
        - Process raw exposure in run()
        - Persist the ISR-corrected exposure as "postISRCCD" if config.doWrite is True

        Parameters
        ----------
        sensorRef : `daf.persistence.butlerSubset.ButlerDataRef`
            DataRef of the detector data to be processed

        Returns
        -------
        result : `pipeBase.Struct`
            Struct contains field "exposure," which is the exposure after application of ISR
        """
        self.log.info("Performing ISR on sensor %s" % (sensorRef.dataId))
        ccdExposure = sensorRef.get('raw')
        camera = sensorRef.get("camera")
        if camera is None and self.config.doAddDistortionModel:
            raise RuntimeError("config.doAddDistortionModel is True "
                               "but could not get a camera from the butler")
        isrData = self.readIsrData(sensorRef, ccdExposure)

        result = self.run(ccdExposure, camera=camera, **isrData.getDict())

        if self.config.doWrite:
            sensorRef.put(result.exposure, "postISRCCD")

        return result

    def convertIntToFloat(self, exposure):
        """Convert an exposure from uint16 to float, set variance plane to 1 and mask plane to 0
        """
        if isinstance(exposure, afwImage.ExposureF):
            # Nothing to be done
            return exposure
        if not hasattr(exposure, "convertF"):
            raise RuntimeError("Unable to convert exposure (%s) to float" % type(exposure))

        newexposure = exposure.convertF()
        maskedImage = newexposure.getMaskedImage()
        varArray = maskedImage.getVariance().getArray()
        varArray[:, :] = 1
        maskArray = maskedImage.getMask().getArray()
        maskArray[:, :] = 0
        return newexposure

    def biasCorrection(self, exposure, biasExposure):
        """Apply bias correction in place

        Parameters
        ----------
        exposure :
            exposure to process
        biasExposure :
            bias exposure of same size as exposure
        """
        isrFunctions.biasCorrection(exposure.getMaskedImage(), biasExposure.getMaskedImage())

    def darkCorrection(self, exposure, darkExposure, invert=False):
        """Apply dark correction in place

        Parameters
        ----------
        exposure :
            exposure to process
        darkExposure :
            dark exposure of same size as exposure
        invert :
            if True, remove the dark from an already-corrected image
        """
        expScale = exposure.getInfo().getVisitInfo().getDarkTime()
        if math.isnan(expScale):
            raise RuntimeError("Exposure darktime is NAN")
        darkScale = darkExposure.getInfo().getVisitInfo().getDarkTime()
        if math.isnan(darkScale):
            raise RuntimeError("Dark calib darktime is NAN")
        isrFunctions.darkCorrection(
            maskedImage=exposure.getMaskedImage(),
            darkMaskedImage=darkExposure.getMaskedImage(),
            expScale=expScale,
            darkScale=darkScale,
            invert=invert
        )

    def doLinearize(self, detector):
        """Is linearization wanted for this detector?

        Checks config.doLinearize and the linearity type of the first amplifier.

        Parameters
        ----------
        detector : `lsst.afw.cameraGeom.Detector`
            detector information (an lsst.afw.cameraGeom.Detector)
        """
        return self.config.doLinearize and \
            detector.getAmpInfoCatalog()[0].getLinearityType() != NullLinearityType

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
        """
        maskPlanes = [self.config.saturatedMaskName, self.config.suspectMaskName]
        gain = amp.getGain()
        if not math.isnan(gain):
            if gain <= 0:
                patchedGain = 1.0
                self.log.warn("Gain for amp %s == %g <= 0; setting to %f" %
                              (amp.getName(), gain, patchedGain))
                gain = patchedGain

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

    def flatCorrection(self, exposure, flatExposure, invert=False):
        """Apply flat correction in place

        Parameters
        ----------
        exposure :
            exposure to process
        flatExposure :
            flatfield exposure same size as exposure
        invert :
            if True, unflatten an already-flattened image instead.
        """
        isrFunctions.flatCorrection(
            maskedImage=exposure.getMaskedImage(),
            flatMaskedImage=flatExposure.getMaskedImage(),
            scalingType=self.config.flatScalingType,
            userScale=self.config.flatUserScale,
            invert=invert
        )

    def getIsrExposure(self, dataRef, datasetType, immediate=True):
        """Retrieve a calibration dataset for removing instrument signature

        Parameters
        ----------
        dataRef :
            data reference for exposure
        datasetType :
            type of dataset to retrieve (e.g. 'bias', 'flat')
        immediate :
            if True, disable butler proxies to enable error
                                        handling within this routine

        Returns
        -------
        exp :
            exposure
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

    def saturationDetection(self, exposure, amp):
        """Detect saturated pixels and mask them using mask plane config.saturatedMaskName, in place

        Parameters
        ----------
        exposure :
            exposure to process; only the amp DataSec is processed
        amp :
            amplifier device data
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
        """Interpolate over saturated pixels, in place

        Parameters
        ----------
        ccdExposure :
            exposure to process

        Notes
        -----
        Warnings :

            - Call saturationDetection first, so that saturated pixels have been identified in the "SAT" mask.
            - Call this after CCD assembly, since saturated regions may cross amplifier boundaries

        """
        isrFunctions.interpolateFromMask(
            maskedImage=ccdExposure.getMaskedImage(),
            fwhm=self.config.fwhm,
            growFootprints=self.config.growSaturationFootprintSize,
            maskName=self.config.saturatedMaskName,
        )

    def suspectDetection(self, exposure, amp):
        """Detect suspect pixels and mask them using mask plane config.suspectMaskName, in place

        Suspect pixels are pixels whose value is greater than amp.getSuspectLevel().
        This is intended to indicate pixels that may be affected by unknown systematics;
        for example if non-linearity corrections above a certain level are unstable
        then that would be a useful value for suspectLevel. A value of `nan` indicates
        that no such level exists and no pixels are to be masked as suspicious.

        Parameters
        ----------
        exposure :
            exposure to process; only the amp DataSec is processed
        amp :
            amplifier device data
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
        """Mask defects using mask plane "BAD" and interpolate over them, in place

        Parameters
        ----------
        ccdExposure :
            exposure to process
        defectBaseList :
            a list of defects to mask and interpolate

        Notes
        -----
        warning : call this after CCD assembly, since defects may cross amplifier boundaries
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

    def maskAndInterpNan(self, exposure):
        """Mask NaNs using mask plane "UNMASKEDNAN" and interpolate over them, in place

        We mask and interpolate over all NaNs, including those
        that are masked with other bits (because those may or may
        not be interpolated over later, and we want to remove all
        NaNs).  Despite this behaviour, the "UNMASKEDNAN" mask plane
        is used to preserve the historical name.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            exposure to process
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

    def overscanCorrection(self, exposure, amp):
        """Apply overscan correction, in-place

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to process; must include both data and bias regions.
        amp : `lsst.afw.table.AmpInfoRecord`
            Amplifier device data.

        Returns
        -------
        result : `lsst.pipe.base.Struct` or `NoneType`
            `None` if there is no overscan; otherwise, this is a
            result struct with components:

            - ``imageFit``: Value(s) removed from image (scalar or
                `lsst.afw.image.Image`).
            - ``overscanFit``: Value(s) removed from overscan (scalar or
                `lsst.afw.image.Image`).
            - ``overscanImage``: Image of the overscan, post-subtraction
                (`lsst.afw.image.Image`).
        """
        if not amp.getHasRawInfo():
            raise RuntimeError("This method must be executed on an amp with raw information.")

        if amp.getRawHorizontalOverscanBBox().isEmpty():
            self.log.info("No Overscan region. Not performing Overscan Correction.")
            return None

        oscanBBox = amp.getRawHorizontalOverscanBBox()

        # afw.cameraGeom.assembleImage.makeUpdatedDetector doesn't update readoutCorner; DM-15559
        x0, x1 = oscanBBox.getBeginX(), oscanBBox.getEndX()

        prescanBBox = amp.getRawPrescanBBox()
        if oscanBBox.getBeginX() > prescanBBox.getBeginX():  # amp is at the right
            x0 += self.config.overscanNumLeadingColumnsToSkip
            x1 -= self.config.overscanNumTrailingColumnsToSkip
        else:
            x0 += self.config.overscanNumTrailingColumnsToSkip
            x1 -= self.config.overscanNumLeadingColumnsToSkip

        oscanBBox = afwGeom.BoxI(afwGeom.PointI(x0, oscanBBox.getBeginY()),
                                 afwGeom.PointI(x1 - 1, oscanBBox.getEndY() - 1))

        maskedImage = exposure.maskedImage
        dataView = maskedImage[amp.getRawDataBBox()]
        overscanImage = maskedImage[oscanBBox]

        results = isrFunctions.overscanCorrection(
            ampMaskedImage=dataView,
            overscanImage=overscanImage,
            fitType=self.config.overscanFitType,
            order=self.config.overscanOrder,
            collapseRej=self.config.overscanRej,
        )
        results.overscanImage = overscanImage
        return results

    def addDistortionModel(self, exposure, camera):
        """Update the WCS in exposure with a distortion model based on camera geometry

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            exposure to process; must include a Detector and a WCS;
            the WCS of the exposure is modified in place
        camera : `lsst.afw.cameraGeom.Camera`
            camera geometry; an lsst.afw.cameraGeom.Camera

        Notes
        -----
        Add a model for optical distortion based on geometry found in `camera`
        and the `exposure`'s detector. The raw input exposure is assumed
        have a TAN WCS that has no compensation for optical distortion.
        Two other possibilities are:
        - The raw input exposure already has a model for optical distortion,
        as is the case for raw DECam data.
        In that case you should set config.doAddDistortionModel False.
        - The raw input exposure has a model for distortion, but it has known
        deficiencies severe enough to be worth fixing (e.g. because they
        cause problems for fitting a better WCS). In that case you should
        override this method with a version suitable for your raw data.

        """
        self.log.info("Adding a distortion model to the WCS")
        wcs = exposure.getWcs()
        if wcs is None:
            raise RuntimeError("exposure has no WCS")
        if camera is None:
            raise RuntimeError("camera is None")
        detector = exposure.getDetector()
        if detector is None:
            raise RuntimeError("exposure has no Detector")
        pixelToFocalPlane = detector.getTransform(PIXELS, FOCAL_PLANE)
        focalPlaneToFieldAngle = camera.getTransformMap().getTransform(FOCAL_PLANE, FIELD_ANGLE)
        distortedWcs = makeDistortedTanWcs(wcs, pixelToFocalPlane, focalPlaneToFieldAngle)
        exposure.setWcs(distortedWcs)

    def setValidPolygonIntersect(self, ccdExposure, fpPolygon):
        """Set the valid polygon as the intersection of fpPolygon and the ccd corners

        Parameters
        ----------
        ccdExposure :
            exposure to process
        fpPolygon :
            Polygon in focal plane coordinates
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

    def brighterFatterCorrection(self, exposure, kernel, maxIter, threshold, applyGain):
        """Apply brighter fatter correction in place for the image

        This correction takes a kernel that has been derived from flat field images to
        redistribute the charge.  The gradient of the kernel is the deflection
        field due to the accumulated charge.

        Notes
        -----
        Given the original image I(x) and the kernel K(x) we can compute the corrected image  Ic(x)
        using the following equation:

        Ic(x) = I(x) + 0.5*d/dx(I(x)*d/dx(int( dy*K(x-y)*I(y))))

        To evaluate the derivative term we expand it as follows:

        0.5 * ( d/dx(I(x))*d/dx(int(dy*K(x-y)*I(y))) + I(x)*d^2/dx^2(int(dy* K(x-y)*I(y))) )

        Because we use the measured counts instead of the incident counts we apply the correction
        iteratively to reconstruct the original counts and the correction.  We stop iterating when the
        summed difference between the current corrected image and the one from the previous iteration
        is below the threshold.  We do not require convergence because the number of iterations is
        too large a computational cost.  How we define the threshold still needs to be evaluated, the
        current default was shown to work reasonably well on a small set of images.  For more information
        on the method see DocuShare Document-19407.

        The edges as defined by the kernel are not corrected because they have spurious values
        due to the convolution.
        """
        self.log.info("Applying brighter fatter correction")

        image = exposure.getMaskedImage().getImage()

        # The image needs to be units of electrons/holes
        with self.gainContext(exposure, image, applyGain):

            kLx = numpy.shape(kernel)[0]
            kLy = numpy.shape(kernel)[1]
            kernelImage = afwImage.ImageD(kLx, kLy)
            kernelImage.getArray()[:, :] = kernel
            tempImage = image.clone()

            nanIndex = numpy.isnan(tempImage.getArray())
            tempImage.getArray()[nanIndex] = 0.

            outImage = afwImage.ImageF(image.getDimensions())
            corr = numpy.zeros_like(image.getArray())
            prev_image = numpy.zeros_like(image.getArray())
            convCntrl = afwMath.ConvolutionControl(False, True, 1)
            fixedKernel = afwMath.FixedKernel(kernelImage)

            # Define boundary by convolution region.  The region that the correction will be
            # calculated for is one fewer in each dimension because of the second derivative terms.
            # NOTE: these need to use integer math, as we're using start:end as numpy index ranges.
            startX = kLx//2
            endX = -kLx//2
            startY = kLy//2
            endY = -kLy//2

            for iteration in range(maxIter):

                afwMath.convolve(outImage, tempImage, fixedKernel, convCntrl)
                tmpArray = tempImage.getArray()
                outArray = outImage.getArray()

                with numpy.errstate(invalid="ignore", over="ignore"):
                    # First derivative term
                    gradTmp = numpy.gradient(tmpArray[startY:endY, startX:endX])
                    gradOut = numpy.gradient(outArray[startY:endY, startX:endX])
                    first = (gradTmp[0]*gradOut[0] + gradTmp[1]*gradOut[1])[1:-1, 1:-1]

                    # Second derivative term
                    diffOut20 = numpy.diff(outArray, 2, 0)[startY:endY, startX + 1:endX - 1]
                    diffOut21 = numpy.diff(outArray, 2, 1)[startY + 1:endY - 1, startX:endX]
                    second = tmpArray[startY + 1:endY - 1, startX + 1:endX - 1]*(diffOut20 + diffOut21)

                    corr[startY + 1:endY - 1, startX + 1:endX - 1] = 0.5*(first + second)

                    tmpArray[:, :] = image.getArray()[:, :]
                    tmpArray[nanIndex] = 0.
                    tmpArray[startY:endY, startX:endX] += corr[startY:endY, startX:endX]

                if iteration > 0:
                    diff = numpy.sum(numpy.abs(prev_image - tmpArray))

                    if diff < threshold:
                        break
                    prev_image[:, :] = tmpArray[:, :]

            if iteration == maxIter - 1:
                self.log.warn("Brighter fatter correction did not converge, final difference %f" % diff)

            self.log.info("Finished brighter fatter in %d iterations" % (iteration + 1))
            image.getArray()[startY + 1:endY - 1, startX + 1:endX - 1] += \
                corr[startY + 1:endY - 1, startX + 1:endX - 1]

    def attachTransmissionCurve(self, exposure, opticsTransmission=None, filterTransmission=None,
                                sensorTransmission=None, atmosphereTransmission=None):
        """Attach a TransmissionCurve to an Exposure, given separate curves for
        different components.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure object to modify by attaching the product of all given
            ``TransmissionCurves`` in post-assembly trimmed detector
            coordinates.  Must have a valid ``Detector`` attached that matches
            the detector associated with sensorTransmission.
        opticsTransmission : `lsst.afw.image.TransmissionCurve`
            A ``TransmissionCurve`` that represents the throughput of the
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

        Notes
        -----
        All ``TransmissionCurve`` arguments are optional; if none are provided,
        the attached ``TransmissionCurve`` will have unit transmission
        everywhere.

        Returns
        -------
        combined : ``lsst.afw.image.TransmissionCurve``
            The TransmissionCurve attached to the exposure.
        """
        return isrFunctions.attachTransmissionCurve(exposure, opticsTransmission=opticsTransmission,
                                                    filterTransmission=filterTransmission,
                                                    sensorTransmission=sensorTransmission,
                                                    atmosphereTransmission=atmosphereTransmission)

    @contextmanager
    def gainContext(self, exp, image, apply):
        """Context manager that applies and removes gain
        """
        if apply:
            ccd = exp.getDetector()
            for amp in ccd:
                sim = image.Factory(image, amp.getBBox())
                sim *= amp.getGain()

        try:
            yield exp
        finally:
            if apply:
                ccd = exp.getDetector()
                for amp in ccd:
                    sim = image.Factory(image, amp.getBBox())
                    sim /= amp.getGain()

    @contextmanager
    def flatContext(self, exp, flat, dark=None):
        """Context manager that applies and removes flats and darks,
        if the task is configured to apply them.
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


class FakeAmp(object):
    """A Detector-like object that supports returning gain and saturation level"""

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

## @addtogroup LSST_task_documentation
## @{
## @page RunIsrTask
## @ref RunIsrTask_ "RunIsrTask"
## @copybrief RunIsrTask
## @}


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
