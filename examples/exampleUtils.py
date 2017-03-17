from builtins import range
from builtins import object
import numpy

from lsst.afw.cameraGeom import DetectorConfig, PIXELS
from lsst.afw.cameraGeom.cameraFactory import makeDetector
import lsst.afw.cameraGeom.utils as cameraGeomUtils
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage


def getReadCorner(flipx, flipy):
    """Get the read corner from flips of the pixel grid
    \param[in] flipx: Flip the x-axis?
    \param[in] flipy: Flip the y-axis?
    \return The read corner in assembled coordinates.
    """
    cornerMap = {(True, True): afwTable.UR,
                 (True, False): afwTable.LR,
                 (False, True): afwTable.UL,
                 (False, False): afwTable.LL}
    return cornerMap[(bool(flipx), bool(flipy))]


def populateAmpBoxes(nx, ny, nprescan, nhoverscan, nvoverscan, nextended, flipx, flipy, ix, iy,
                     isPerAmp, record):
    '''!Fill ampInfo tables
    \param[in] isPerAmp -- If True, return a dictionary of amp exposures keyed by amp name.
                           If False, return a single exposure with amps mosaiced preserving non-science pixels
                           (e.g. overscan)
    \param[in] nx -- number of pixels in the serial register
    \param[in] ny -- number of rows in the parallel direction
    \param[in] nprescan -- number of prescan rows
    \param[in] nhoverscan -- number of horizonatal overscan columns
    \param[in] nvoverscan -- number of vertical overscan rows
    \param[in] nextended -- number of pixels in the extended register
    \param[in] flipx -- should the amp be flipped about the x axis when inserted in the chip mosaic?
    \param[in] flipy -- should the amp be flipped about the y axis when inserted int he chip mosaic?
    \param[in] ix -- index in x direction of the amp in the chip
    \param[in] iy -- index in y direction of the amp in the chip
    \param[in] isPerAmp -- are the raw data per amp or assembled into a mosaiced image
    \param[in, out] record -- record to add this amp to
    '''
    def makeBbox(x0, y0, x_extent, y_extent):
        return afwGeom.BoxI(afwGeom.PointI(x0, y0), afwGeom.ExtentI(x_extent, y_extent))

    bbox = makeBbox(0, 0, nx, ny)

    dataBox = makeBbox(0, 0, nx, ny)
    dataBox.shift(afwGeom.ExtentI(nextended, nprescan))

    allBox = afwGeom.BoxI()

    preBox = makeBbox(0, 0, nx, nprescan)
    preBox.shift(afwGeom.ExtentI(nextended, 0))

    extBox = makeBbox(0, 0, nextended, ny)
    extBox.shift(afwGeom.ExtentI(0, nprescan))

    hOscanBox = makeBbox(0, 0, nhoverscan, ny)
    hOscanBox.shift(afwGeom.ExtentI(nextended+nx, nprescan))

    vOscanBox = makeBbox(0, 0, nx, nvoverscan)
    vOscanBox.shift(afwGeom.ExtentI(nextended, nprescan+ny))

    allBox.include(dataBox)
    allBox.include(preBox)
    allBox.include(extBox)
    allBox.include(hOscanBox)
    allBox.include(vOscanBox)

    bbox.shift(afwGeom.ExtentI(ix*nx, iy*ny))
    xtot = allBox.getDimensions().getX()
    ytot = allBox.getDimensions().getY()
    rShiftExt = afwGeom.ExtentI(ix*xtot, iy*ytot)
    if not isPerAmp:
        #Set read corner in assembled coordinates
        record.setReadoutCorner(getReadCorner(flipx, flipy))

        allBox.shift(rShiftExt)

        if flipx:
            dataBox.flipLR(xtot)
            preBox.flipLR(xtot)
            extBox.flipLR(xtot)
            hOscanBox.flipLR(xtot)
            vOscanBox.flipLR(xtot)
            flipx = False
        if flipy:
            dataBox.flipTB(ytot)
            preBox.flipTB(ytot)
            extBox.flipTB(ytot)
            hOscanBox.flipTB(ytot)
            vOscanBox.flipTB(ytot)
            flipy = False

        dataBox.shift(rShiftExt)
        preBox.shift(rShiftExt)
        extBox.shift(rShiftExt)
        hOscanBox.shift(rShiftExt)
        vOscanBox.shift(rShiftExt)
        rawXoff = 0
        rawYoff = 0

    else:
        #We assume that single amp images have the first pixel read in the
        #lower left and that the pixels are arrange such that the
        #serial is along the x-axis.
        record.setReadoutCorner(afwTable.LL)
        rawXoff = rShiftExt.getX()
        rawYoff = rShiftExt.getY()

    record.setBBox(bbox)
    record.setName("A:%i,%i"%(ix, iy))

    record.setGain(1.)
    record.setSaturation(100000)
    record.setReadNoise(1.)
    record.setLinearityCoeffs((0., 1., 0., 0.))
    record.setLinearityType('Polynomial')
    record.setHasRawInfo(True)
    record.setRawFlipX(flipx)
    record.setRawFlipY(flipy)
    record.setRawBBox(allBox)
    record.setRawXYOffset(afwGeom.ExtentI(rawXoff, rawYoff))
    record.setRawDataBBox(dataBox)
    record.setRawHorizontalOverscanBBox(hOscanBox)
    record.setRawVerticalOverscanBBox(vOscanBox)
    record.setRawPrescanBBox(preBox)


def createDetector(nAmpX, nAmpY, nPixX, nPixY, pre, hOscan, vOscan, ext, isPerAmp):
    '''!Fill ampInfo tables
    \param[in] nAmpX -- Number of amps in the x direction
    \param[in] nAmpY -- Number of amps in the y direction
    \param[in] nPixX -- Number of pixels in the amp in the x direction
    \param[in] nPixY -- Number of pixels in the amp in the y direction
    \param[in] pre -- Number of prescan rows
    \param[in] hOscan -- Number of horizontal overscan columns
    \param[in] vOscan -- Number of vertical overscan rows
    \param[in] ext -- Number of pixels in the extended register
    \param[in] isPerAmp -- Are the raw amp data in separate images?
    \return an lsst.afw.cameraGeom.Detector object
    '''
    schema = afwTable.AmpInfoTable.makeMinimalSchema()
    ampCatalog = afwTable.AmpInfoCatalog(schema)
    flipy = True
    for iy in range(nAmpY):
        flipy = not flipy
        flipx = True
        for ix in range(nAmpX):
            flipx = not flipx
            record = ampCatalog.addNew()
            populateAmpBoxes(nPixX, nPixY, pre, hOscan, vOscan, ext, flipx, flipy, ix, iy,
                             isPerAmp, record)
            record.setGain(ix+iy*nAmpX+1.)

    detConfig = DetectorConfig()
    detConfig.name = 'TestDetector'
    detConfig.id = 0
    detConfig.bbox_x0 = 0
    detConfig.bbox_y0 = 0
    detConfig.bbox_x1 = nAmpX*nPixX - 1
    detConfig.bbox_y1 = nAmpY*nPixY - 1
    detConfig.detectorType = 0 #Science type
    detConfig.serial = 'THX1138'
    detConfig.offset_x = 0.
    detConfig.offset_y = 0.
    detConfig.refpos_x = nAmpX*nPixX*0.5 - 0.5
    detConfig.refpos_y = nAmpY*nPixY*0.5 - 0.5
    detConfig.yawDeg = 0.
    detConfig.pitchDeg = 0.
    detConfig.rollDeg = 0.
    detConfig.pixelSize_x = 10./1000. #in mm
    detConfig.pixelSize_y = 10./1000. #in mm
    detConfig.transposeDetector = False
    detConfig.transformDict.nativeSys = PIXELS.getSysName()

    fpTransform = afwGeom.xyTransformRegistry['identity']()
    return makeDetector(detConfig, ampCatalog, fpTransform)


def makeFakeWcs():
    '''!Make a wcs to put in an exposure
    \return a Wcs object
    '''
    return afwImage.makeWcs(afwCoord.IcrsCoord(45.0*afwGeom.degrees, 45.0*afwGeom.degrees),
                            afwGeom.Point2D(0.0, 0.0), 1.0, 0.0, 0.0, 1.0)


def makeExpFromIm(im, detector):
    wcs = makeFakeWcs()
    var = afwImage.ImageF(im)
    mask = afwImage.MaskU(im.getDimensions())
    mi = afwImage.makeMaskedImage(im, mask, var)
    exp = afwImage.makeExposure(mi)
    exp.setDetector(detector)
    exp.setWcs(wcs)
    return exp


def makeAmpInput(detector):
    '''!Make a dictionary of amp images for assembly
    \param[in] detector -- An lsst.afw.cameraGeom.Detector describing the detector to create
    \return a dictionary of amp exposures keyed on the amp names
    '''
    inputData = {}
    for amp in detector:
        im = cameraGeomUtils.makeImageFromAmp(amp, imageFactory=afwImage.ImageF)
        exp = makeExpFromIm(im, detector)
        inputData[amp.getName()] = exp
    return inputData


def makeAssemblyInput(isPerAmp, doTrim=False):
    '''!Make the input to pass to the assembly task
    \param[in] isPerAmp -- If True, return a dictionary of amp exposures keyed by amp name.
                           If False, return a single exposure with amps mosaiced preserving non-science pixels
                           (e.g. overscan)
    \param[in doTrim -- Trim out the non-data pixels (e.g. overscan)?  Ignored if isPerAmp is True.
    \return Either a dictionary of amp exposures or an exposure contining the mosaiced amps.
    '''

    #number of amps in x and y
    nAmpX = 3
    nAmpY = 2

    #number of science pixels in each amp in x and y
    nPixX = 512
    nPixY = 1024

    #number of prescan rows
    pre = 4

    #number of horizontal overscan columns
    hOscan = 10

    #number of vertical overscan rows
    vOscan = 15

    #number of pixels in the extended register
    ext = 1

    detector = createDetector(nAmpX, nAmpY, nPixX, nPixY, pre, hOscan, vOscan, ext, isPerAmp)
    if isPerAmp:
        return makeAmpInput(detector)
    else:
        im = cameraGeomUtils.makeImageFromCcd(detector, isTrimmed=doTrim, imageFactory=afwImage.ImageF)
        ccdAssemblyInput = makeExpFromIm(im, detector)
        ccdAssemblyInput.setDetector(detector)
        return ccdAssemblyInput


def makeRaw(darkval, oscan, gradient, exptime):
    '''!Make a raw image for input to ISR
    \param[in] darkval -- dark current e-/sec
    \param[in] oscan -- overscan value
    \param[in] gradient -- fractional gradient in the flat
    \param[in] exptime -- exposure time of this observation
    \return a raw exposure containing mosaiced raw amps
    '''
    rawExposure = makeAssemblyInput(False)
    detector = rawExposure.getDetector()
    visitInfo = afwImage.VisitInfo(exposureTime=exptime)
    rawExposure.getInfo().setVisitInfo(visitInfo)
    im = rawExposure.getMaskedImage().getImage()
    for amp in detector:
        subim = im.Factory(im, amp.getRawDataBBox())
        subim.set(5000.)
        subArr = subim.getArray()
        nPixY = subim.getDimensions().getY()
        grad = numpy.interp(range(nPixY), (0, nPixY-1), (1., 1.-gradient))
        subArr *= grad[:, numpy.newaxis]
        subArr += darkval*exptime
        subArr /= amp.getGain()
        subArr += oscan
        suboscanim = im.Factory(im, amp.getRawHorizontalOverscanBBox())
        suboscanim.set(oscan)
    return rawExposure


def makeDark(darkval, exptime):
    '''!Make a dark exposure in DN
    \param[in] darkval -- dark current in e-/sec
    \param[in] exptime -- exposure time of the dark frame
    \return an assembled dark exposure
    '''
    darkExposure = makeAssemblyInput(False, doTrim=True)
    detector = darkExposure.getDetector()
    visitInfo = afwImage.VisitInfo(exposureTime=exptime)
    darkExposure.getInfo().setVisitInfo(visitInfo)
    im = darkExposure.getMaskedImage().getImage()
    for amp in detector:
        subim = im.Factory(im, amp.getBBox())
        subim.set(darkval*exptime/amp.getGain())
    return darkExposure


def makeFlat(gradient):
    '''!Make a flat exposure including gain variation
    \param[in] gradient -- fractional gradient in the flat from bottom to top
    \return an assembled flat exposure
    '''
    flatExposure = makeAssemblyInput(False, doTrim=True)
    detector = flatExposure.getDetector()
    im = flatExposure.getMaskedImage().getImage()
    for amp in detector:
        subim = im.Factory(im, amp.getBBox())
        subim.set(1.)
        subArr = subim.getArray()
        nPixY = subim.getDimensions().getY()
        grad = numpy.interp(range(nPixY), (0, nPixY-1), (1., 1.-gradient))
        subArr *= grad[:, numpy.newaxis]
        subArr /= amp.getGain()
    return flatExposure


class FakeDataRef(object):
    '''!A mock data reference to use in the example IsrTask runner
    The main thing is to define the get method with the relevant datatypes.
    This can be extended to mimic other getters (fringe, e.g.) if needed.
    '''
    darkval = 2. #e-/sec
    oscan = 1000. #DN
    gradient = .10
    exptime = 15 #seconds
    darkexptime = 40. #seconds
    dataId = "My Fake Data"

    def get(self, dataType, **kwargs):
        if dataType == 'raw':
            return makeRaw(self.darkval, self.oscan, self.gradient, self.exptime)
        if dataType == 'dark':
            return makeDark(self.darkval, self.darkexptime)
        if dataType == 'flat':
            return makeFlat(self.gradient)
        if dataType == 'defects':
            return []

    def put(self, exposure, filename):
        exposure.writeFits(filename+".fits")
