import os
import re
import lsst.pex.policy as pexPolicy
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.cameraGeom.utils as cameraGeomUtils
import lsst.afw.display.ds9 as ds9

class ccdImageFactory(cameraGeomUtils.GetCcdImage):
    """A class to return an Image of a given Simulated CCD"""
    
    def __init__(self, frameId, isTrimmed=True, imageFactory=afwImage.ImageU):
        self.frameId = frameId
        self.isTrimmed = isTrimmed
        self.imageFactory = imageFactory

    def getFilename(self, ccd, amp, expType=None):
        """Return the filename of specified Ccd"""
        
        mat = re.search(r"^R:(\d),(\d)\s+S:(\d),(\d)\s*$", ccd.getId().getName())
        if mat:
            raftId = "R%s%s" % (mat.group(1), mat.group(2))
            sensorId = "S%s%s" % (mat.group(3), mat.group(4))

        ampTuple = amp.getId().getIndex()
        ampName = "C%d%d" % (ampTuple[1], ampTuple[0])

        fileName = "%s/ImSim/processed/imsim_%08d_%s_%s_%s_E000.fits" % (
                os.environ['AFWDATA_DIR'],
                self.frameId, raftId, sensorId, ampName)

        return fileName

    def getImage(self, ccd, amp, expType=None, imageFactory=None):
        """Return an image of the specified Amp in the specified Ccd"""

        if not imageFactory:
            imageFactory = self.imageFactory

        return self.getImageFromFilename(self.getFilename(ccd, amp), ccd, amp, imageFactory=imageFactory)

def getGeomPolicy(cameraGeomPolicyFile):
    """Return a Policy describing a Camera's geometry"""

    policyFile = pexPolicy.DefaultPolicyFile("afw", "CameraGeomDictionary.paf", "policy")
    defPolicy = pexPolicy.Policy.createPolicy(policyFile, policyFile.getRepositoryPath(), True)

    geomPolicy = pexPolicy.Policy.createPolicy(cameraGeomPolicyFile, True)

    geomPolicy.mergeDefaults(defPolicy.getDictionary())

    return geomPolicy

def getCamera(cameraGeomPolicyFile):
    """Return a Camera object given the camera description PAF file

If oneAmpPerFile is True, assume that each Amp is written to a separate file
    """
    geomPolicy = getGeomPolicy(cameraGeomPolicyFile)

    return cameraGeomUtils.makeCamera(geomPolicy)

def getCcdId(ccdId):
    if isinstance(ccdId, cameraGeom.Id):
        if ccdId.getSerial() > 0:
            ccdId = ccdId.getSerial()
        else:
            ccdId = ccdId.getName()

    return ccdId

def getCcd(camera, ccdId):
    """Return a Camera object given a Camera"""
    return cameraGeomUtils.findCcd(camera, cameraGeomUtils.cameraGeom.Id(getCcdId(ccdId)))

def showAmps(camera, ccdId, perFile=False):
    ccd = cameraGeomUtils.findCcd(camera, cameraGeomUtils.cameraGeom.Id(getCcdId(ccdId)))

    for a in ccd:
        print "%-10s %s" % (a.getId(), a.getDataSec(False))

def foo(frameId=85751839, ccdName="R:2,3 S:1,1", geomPolicyFile="tests/Full_STA_geom.paf",
        isTrimmed=False, display=True):
    cif = ccdImageFactory(frameId)
    camera = getCamera(geomPolicyFile)
    raft = cameraGeomUtils.findRaft(camera, cameraGeomUtils.cameraGeom.Id(23,"R:2,3"))
    ccd = getCcd(raft, ccdName)
    if False:
        showAmps(camera, ccdName)
    ccd.setTrimmed(isTrimmed)
    ccdImage = cameraGeomUtils.makeImageFromCcd(ccd, cif, isTrimmed=isTrimmed)
    for a in ccd:
        biasbbox = a.getBiasSec()
        print a.getAllPixels()
        print biasbbox.getX0(), biasbbox.getY0, biasbbox.shift(-biasbbox.getX0(), -biasbbox.getY0()), a.getDataSec()

    if display:
        cameraGeomUtils.showCcd(ccd, ccdImage, isTrimmed=isTrimmed, frame=0)

    return ccd, ccdImage

if __name__ == "__main__":
    foo(display=False)
