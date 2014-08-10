from lsst.afw.cameraGeom import DetectorConfig, PIXELS
from lsst.afw.cameraGeom.cameraFactory import makeDetector
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage

def populateAmpBoxes(nx, ny, nprescan, nhoverscan, nvoverscan, nextended, flipx, flipy, ix, iy,
                      isPerAmp, record):
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
        rawXoff = rShiftExt.getX()
        rawYoff = rShiftExt.getY()

    record.setBBox(bbox)
    record.setName("A:%i,%i"%(ix, iy))

    #The readout corner is in the assemble coordinates
    if flipx and flipy:
        record.setReadoutCorner(afwTable.UR)
    elif flipx and not flipy:
        record.setReadoutCorner(afwTable.LR)
    elif not flipx and flipy:
        record.setReadoutCorner(afwTable.UL)
    elif not flipx and not flipy:
        record.setReadoutCorner(afwTable.LL)
    else:
        raise ValueError("Cannont determine read corner given flipx: %s, flipy: %s"%(flipx, flipy))

    record.setGain(1.)
    record.setSaturation(1)
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
    plateScale = 1.
    return makeDetector(detConfig, ampCatalog, fpTransform, plateScale)

def makeFakeAmp(amp):
    im = afwImage.ImageF(amp.getRawBBox().getDimensions())
    im.set(amp.getGain())
    markBox = afwGeom.BoxI(amp.getRawDataBBox().getMin(), afwGeom.ExtentI(10, 10))
    subim = afwImage.ImageF(im, markBox)
    subim.set(0)
    return im

def makeAmpInput(detector):
    inputData = {}
    for amp in detector:
        im = makeFakeAmp(amp)
        var = afwImage.ImageF(im)
        mask = afwImage.MaskU(im.getDimensions())
        mi = afwImage.makeMaskedImage(im, mask, var)
        exp = afwImage.makeExposure(mi)
        exp.setDetector(detector)
        inputData[amp.getName()] = exp
    return inputData
