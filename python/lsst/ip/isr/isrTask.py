#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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
import lsst.afw.image as afwImage
import lsst.meas.algorithms as measAlg
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from . import isr
from .isrLib import maskNans
from .assembleCcdTask import AssembleCcdTask
from .fringe import FringeTask

class IsrTaskConfig(pexConfig.Config):
    doBias = pexConfig.Field(
        dtype = bool,
        doc = "Apply bias frame correction?",
        default = True,
    )
    doDark = pexConfig.Field(
        dtype = bool,
        doc = "Apply dark frame correction?",
        default = True,
    )
    doFlat = pexConfig.Field(
        dtype = bool,
        doc = "Apply flat field correction?",
        default = True,
    )
    doFringe = pexConfig.Field(
        dtype = bool,
        doc = "Apply fringe correction?",
        default = True,
        )
    doWrite = pexConfig.Field(
        dtype = bool,
        doc = "Persist postISRCCD?",
        default = True,
    )
    assembleCcd = pexConfig.ConfigurableField(
        target = AssembleCcdTask,
        doc = "CCD assembly task",
    )
    fringeAfterFlat = pexConfig.Field(
        dtype = bool,
        doc = "Do fringe subtraction after flat-fielding?",
        default = True,
        )
    fringe = pexConfig.ConfigurableField(
        target = FringeTask,
        doc = "Fringe subtraction task",
        )
    fwhm = pexConfig.Field(
        dtype = float,
        doc = "FWHM of PSF (arcsec)",
        default = 1.0,
    )
    saturatedMaskName = pexConfig.Field(
        dtype = str,
        doc = "Name of mask plane to use in saturation detection and interpolation",
        default = "SAT",
    )
    flatScalingType = pexConfig.ChoiceField(
        dtype = str,
        doc = "The method for scaling the flat on the fly.",
        default = 'USER',
        allowed = {
            "USER":   "Scale by flatUserScale",
            "MEAN":   "Scale by the inverse of the mean",
            "MEDIAN": "Scale by the inverse of the median",
        },
    )
    flatUserScale = pexConfig.Field(
        dtype = float,
        doc = "If flatScalingType is 'USER' then scale flat by this amount; ignored otherwise",
        default = 1.0,
    )
    overscanFitType = pexConfig.ChoiceField(
        dtype = str,
        doc = "The method for fitting the overscan bias level.",
        default = 'MEDIAN',
        allowed = {
            "POLY": "Fit ordinary polynomial to the longest axis of the overscan region",
            "CHEB": "Fit Chebyshev polynomial to the longest axis of the overscan region",
            "LEG":  "Fit Legendre polynomial to the longest axis of the overscan region",
            "NATURAL_SPLINE": "Fit natural spline to the longest axis of the overscan region",
            "CUBIC_SPLINE": "Fit cubic spline to the longest axis of the overscan region",
            "AKIMA_SPLINE": "Fit Akima spline to the longest axis of the overscan region",
            "MEAN": "Correct using the mean of the overscan region",
            "MEDIAN": "Correct using the median of the overscan region",
        },
    )
    overscanOrder = pexConfig.Field(
        dtype = int,
        doc = ("Order of polynomial or to fit if overscan fit type is a polynomial, " +
               "or number of spline knots if overscan fit type is a spline."),
        default = 1,
    )
    overscanRej = pexConfig.Field(
        dtype = float,
        doc = "Rejection threshold (sigma) for collapsing overscan before fit",
        default = 3.0,
    )
    growSaturationFootprintSize = pexConfig.Field(
        dtype = int,
        doc = "Number of pixels by which to grow the saturation footprints",
        default = 1,
    )
    fluxMag0T1 = pexConfig.Field(
        dtype = float,
        doc = "The approximate flux of a zero-magnitude object in a one-second exposure",
        default = 1e10,
    )
    setGainAssembledCcd = pexConfig.Field(
        dtype = bool,
        doc = "update exposure metadata in the assembled ccd to reflect the effective gain of the assembled chip",
        default = True,
    )
    keysToRemoveFromAssembledCcd = pexConfig.ListField(
        dtype = str,
        doc = "fields to remove from the metadata of the assembled ccd.",
        default = [],
    )
    doAssembleDetrends = pexConfig.Field(
        dtype = bool,
        default = False,
        doc = "Assemble detrend/calibration frames?"
        )

## \addtogroup LSST_task_documentation
## \{
## \page IsrTask
## \ref IsrTask_ "IsrTask"
## \copybrief IsrTask
## \}

class IsrTask(pipeBase.CmdLineTask):
    """!
    \anchor IsrTask_

    \brief Apply common instrument signature correction algorithms to a raw frame.

    \section ip_isr_isr_Contents Contents

     - \ref ip_isr_isr_Purpose
     - \ref ip_isr_isr_Initialize
     - \ref ip_isr_isr_IO
     - \ref ip_isr_isr_Config
     - \ref ip_isr_isr_Debug
     - \ref ip_isr_isr_Example

    \section ip_isr_isr_Purpose Description

    The process for correcting imaging data is very similar from camera to camera.
    This task provides a vanilla implementation of doing these corrections, including 
    the ability to turn certain corrections off if they are not needed.  The input
    is a daf.persistence.butlerSubset.ButlerDataRef.  The data reference can return
    the raw input and all the calibration products.  The raw input is a single chip
    sized mosaic of all amps including overscans and other non-science pixels. This
    task may not meet all needs and it is expected that it will be subclassed for specific
    applications.

    \section ip_isr_isr_Initialize Task initialization

    \copydoc \_\_init\_\_

    \section ip_isr_isr_IO Inputs/Outputs to the assembleCcd method

    \copydoc run

    \section ip_isr_isr_Config Configuration parameters

    See \ref IsrTaskConfig

    \section ip_isr_isr_Debug Debug variables

    The \link lsst.pipe.base.cmdLineTask.CmdLineTask command line task\endlink interface supports a
    flag \c -d to import \b debug.py from your \c PYTHONPATH; see <a
    href="http://lsst-web.ncsa.illinois.edu/~buildbot/doxygen/x_masterDoxyDoc/base_debug.html">
    Using lsstDebug to control debugging output</a> for more about \b debug.py files.

    The available variables in AssembleCcdTask are:
    <DL>
      <DT> \c display
      <DD> A dictionary containing debug point names as keys with frame number as value. Valid keys are:
        <DL>
          <DT> postISRCCD 
          <DD> display exposure after ISR has been applied 
        </DL>
    </DL>  

    \section ip_isr_isr_Example A complete example of using IsrTask

    This code is in runIsrTask.py in the examples directory, and can be run as \em e.g.
    \code
    python examples/runIsrTask.py
    \endcode
<HR>
    \dontinclude runIsrTask.py
    Import the task.  There are other imports.  Read the source file for more info.
    \skipline IsrTask

    \dontinclude exampleUtils.py
    Create the input data reference with the help of some utilities in examples/exampleUtils.py.  This
    is a mock data reference that has all the right methods to run through ISR.  We will only
    do overscan, dark and flat correction, so it only needs to know how to get those products (and an
    empty list of defects).
    \skip FakeDataRef
    \until writeFits
    The above numbers can be changed to modify the gradient in the flat, for example.  

    \dontinclude exampleUtils.py
    The data are constructed by hand so that all effects will be corrected for essentially perfectly
    \skip makeRaw
    \until return flatExposure


    \dontinclude runIsrTask.py
    Construct the task and set some config parameters.  Specifically, we don't want to
    do zero or fringe correction.  We also don't want the assembler messing with the gain.
    \skip runIsr
    \until config=isrConfig

    Now make the fake data reference and run it through ISR.
    \skip sensorRef
    \until isrTask.run

    <HR>
    To investigate the \ref ip_isr_isr_Debug, put something like
    \code{.py}
    import lsstDebug
    def DebugInfo(name):
        di = lsstDebug.getInfo(name)        # N.b. lsstDebug.Info(name) would call us recursively
        if name == "lsst.ip.isr.isrTask":
            di.display = {'postISRCCD':2} 
        return di

    lsstDebug.Info = DebugInfo
    \endcode
    into your debug.py file and run runAssembleTask.py with the \c --debug flag.


    Conversion notes:
        Display code should be updated once we settle on a standard way of controlling what is displayed.
    """
    ConfigClass = IsrTaskConfig
    _DefaultName = "isr"

    def __init__(self, *args, **kwargs):
        '''!Constructor for IsrTask
        \param[in] *args -- a list of positional arguments passed on to the Task constructor
        \param[in] **kwargs -- a dictionary of keyword arguments passed on to the Task constructor
        Call the lsst.pipe.base.task.Task.__init__ method
        Then setup the assembly and fringe correction subtasks
        '''
        pipeBase.Task.__init__(self, *args, **kwargs)
        self.makeSubtask("assembleCcd")
        self.makeSubtask("fringe")

    @pipeBase.timeMethod
    def run(self, sensorRef):
        """!Perform instrument signature removal on an exposure

        Steps include:
        - Detect saturation, apply overscan correction, bias, dark and flat
        - Perform CCD assembly
        - Interpolate over defects, saturated pixels and all NaNs
        - Persist the ISR-corrected exposure as "postISRCCD" if config.doWrite is True

        \param[in] sensorRef -- daf.persistence.butlerSubset.ButlerDataRef of the detector data to be processed
        \return a pipeBase.Struct with fields:
        - exposure: the exposure after application of ISR
        """
        self.log.log(self.log.INFO, "Performing ISR on sensor %s" % (sensorRef.dataId))
        ccdExposure = sensorRef.get('raw')
        ccd = ccdExposure.getDetector()

        ccdExposure = self.convertIntToFloat(ccdExposure)
        for amp in ccd:
            self.saturationDetection(ccdExposure, amp)

            self.overscanCorrection(ccdExposure, amp)

        ccdExposure = self.assembleCcd.assembleCcd(ccdExposure)

        ccd = ccdExposure.getDetector()

        if self.config.doBias:
            biasExposure = self.getDetrend(sensorRef, "bias")
            self.biasCorrection(ccdExposure, biasExposure)

        if self.config.doDark:
            darkExposure = self.getDetrend(sensorRef, "dark")
            self.darkCorrection(ccdExposure, darkExposure)

        for amp in ccd:
            ampExposure = ccdExposure.Factory(ccdExposure, amp.getBBox())

            self.updateVariance(ampExposure, amp)

        if self.config.doFringe and not self.config.fringeAfterFlat:
            self.fringe.run(ccdExposure, sensorRef,
                            assembler=self.assembleCcd if self.config.doAssembleDetrends else None)

        if self.config.doFlat:
            flatExposure = self.getDetrend(sensorRef, "flat")
            self.flatCorrection(ccdExposure, flatExposure)

        defects = sensorRef.get('defects')
        self.maskAndInterpDefect(ccdExposure, defects)

        self.saturationInterpolation(ccdExposure)

        self.maskAndInterpNan(ccdExposure)

        if self.config.doFringe and self.config.fringeAfterFlat:
            self.fringe.run(ccdExposure, sensorRef,
                            assembler=self.assembleCcd if self.config.doAssembleDetrends else None)

        ccdExposure.getCalib().setFluxMag0(self.config.fluxMag0T1 * ccdExposure.getCalib().getExptime())

        if self.config.doWrite:
            sensorRef.put(ccdExposure, "postISRCCD")

        self.display("postISRCCD", ccdExposure)

        return pipeBase.Struct(
            exposure = ccdExposure,
        )

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
        varArray[:,:] = 1
        maskArray = maskedImage.getMask().getArray()
        maskArray[:,:] = 0
        return newexposure

    def biasCorrection(self, exposure, biasExposure):
        """!Apply bias correction in place

        \param[in,out]  exposure        exposure to process
        \param[in]      biasExposure    bias exposure of same size as exposure
        """
        isr.biasCorrection(exposure.getMaskedImage(), biasExposure.getMaskedImage())

    def darkCorrection(self, exposure, darkExposure):
        """!Apply dark correction in place

        \param[in,out]  exposure        exposure to process
        \param[in]      darkExposure    dark exposure of same size as exposure
        """
        darkCalib = darkExposure.getCalib()
        isr.darkCorrection(
            maskedImage = exposure.getMaskedImage(),
            darkMaskedImage = darkExposure.getMaskedImage(),
            expScale = exposure.getCalib().getExptime(),
            darkScale = darkCalib.getExptime(),
        )

    def updateVariance(self, ampExposure, amp):
        """!Set the variance plane based on the image plane, plus amplifier gain and read noise

        \param[in,out]  ampExposure     exposure to process
        \param[in]      amp             amplifier detector information
        """
        isr.updateVariance(
            maskedImage = ampExposure.getMaskedImage(),
            gain = amp.getGain(),
            readNoise = amp.getReadNoise(),
        )

    def flatCorrection(self, exposure, flatExposure):
        """!Apply flat correction in place

        \param[in,out]  exposure        exposure to process
        \param[in]      flatExposure    flatfield exposure same size as exposure
        """
        isr.flatCorrection(
            maskedImage = exposure.getMaskedImage(),
            flatMaskedImage = flatExposure.getMaskedImage(),
            scalingType = self.config.flatScalingType,
            userScale = self.config.flatUserScale,
        )

    def getDetrend(self, dataRef, detrend, immediate=True):
        """!Get a detrend exposure

        \param[in]      dataRef         data reference for exposure
        \param[in]      detrend         detrend/calibration to read
        \param[in]      immediate       if True, disable butler proxies to enable error
                                        handling within this routine
        \return Detrend exposure
        """
        try:
            exp = dataRef.get(detrend, immediate=immediate)
        except Exception, e:
            raise RuntimeError("Unable to retrieve %s for %s: %s" % (detrend, dataRef.dataId, e))
        if self.config.doAssembleDetrends:
            exp = self.assembleCcd.assembleCcd(exp)
        return exp

    def saturationDetection(self, exposure, amp):
        """!Detect saturated pixels and mask them using mask plane "SAT", in place

        \param[in,out]  exposure    exposure to process; only the amp DataSec is processed
        \param[in]      amp         amplifier device data
        """
        maskedImage = exposure.getMaskedImage()
        dataView = maskedImage.Factory(maskedImage, amp.getRawBBox())
        isr.makeThresholdMask(
            maskedImage = dataView,
            threshold = amp.getSaturation(),
            growFootprints = 0,
            maskName = self.config.saturatedMaskName,
        )

    def saturationInterpolation(self, ccdExposure):
        """!Interpolate over saturated pixels, in place

        \param[in,out]  ccdExposure     exposure to process

        \warning:
        - Call saturationDetection first, so that saturated pixels have been identified in the "SAT" mask.
        - Call this after CCD assembly, since saturated regions may cross amplifier boundaries
        """
        isr.interpolateFromMask(
            maskedImage = ccdExposure.getMaskedImage(),
            fwhm = self.config.fwhm,
            growFootprints = self.config.growSaturationFootprintSize,
            maskName = self.config.saturatedMaskName,
        )

    def maskAndInterpDefect(self, ccdExposure, defectBaseList):
        """!Mask defects using mask plane "BAD" and interpolate over them, in place

        \param[in,out]  ccdExposure     exposure to process
        \param[in] defectBaseList a list of defects to mask and interpolate

        \warning: call this after CCD assembly, since defects may cross amplifier boundaries
        """
        maskedImage = ccdExposure.getMaskedImage()
        defectList = measAlg.DefectListT()
        for d in defectBaseList:
            bbox = d.getBBox()
            nd = measAlg.Defect(bbox)
            defectList.append(nd)
        isr.maskPixelsFromDefectList(maskedImage, defectList, maskName='BAD')
        isr.interpolateDefectList(
            maskedImage = maskedImage,
            defectList = defectList,
            fwhm = self.config.fwhm,
        )

    def maskAndInterpNan(self, exposure):
        """!Mask NaNs using mask plane "UNMASKEDNAN" and interpolate over them, in place

        We mask and interpolate over all NaNs, including those
        that are masked with other bits (because those may or may
        not be interpolated over later, and we want to remove all
        NaNs).  Despite this behaviour, the "UNMASKEDNAN" mask plane
        is used to preserve the historical name.

        \param[in,out]  exposure        exposure to process
        """
        maskedImage = exposure.getMaskedImage()

        # Find and mask NaNs
        maskedImage.getMask().addMaskPlane("UNMASKEDNAN") 
        maskVal = maskedImage.getMask().getPlaneBitMask("UNMASKEDNAN")
        numNans = maskNans(maskedImage, maskVal)
        self.metadata.set("NUMNANS", numNans)

        # Interpolate over these previously-unmasked NaNs
        if numNans > 0:
            self.log.log(self.log.WARN, "There were %i unmasked NaNs" % (numNans,))
            nanDefectList = isr.getDefectListFromMask(
                maskedImage = maskedImage,
                maskName = 'UNMASKEDNAN',
                growFootprints = 0,
            )
            isr.interpolateDefectList(
                maskedImage = exposure.getMaskedImage(),
                defectList = nanDefectList,
                fwhm = self.config.fwhm,
            )

    def overscanCorrection(self, exposure, amp):
        """!Apply overscan correction, in place

        \param[in,out]  exposure    exposure to process; must include both DataSec and BiasSec pixels
        \param[in]      amp         amplifier device data
        """
        if not amp.getHasRawInfo():
            raise RuntimeError("This method must be executed on an amp with raw information.")
        maskedImage = exposure.getMaskedImage()
        dataView = maskedImage.Factory(maskedImage, amp.getRawDataBBox())

        expImage = exposure.getMaskedImage().getImage()
        overscanImage = expImage.Factory(expImage, amp.getRawHorizontalOverscanBBox())

        isr.overscanCorrection(
            ampMaskedImage = dataView,
            overscanImage = overscanImage,
            fitType = self.config.overscanFitType,
            order = self.config.overscanOrder,
            collapseRej = self.config.overscanRej,
        )
