import sys
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath

def getTemplateSubexposure(scienceExposure, templateFile, growSize=20, remap=True):
    # Template information
    templateMetadata = afwImage.readMetadata(templateFile)
    templateWcs      = afwImage.Wcs(templateMetadata)

    # Find overlap of Science exposure with Template
    templateBBox     = findOverlapBBox(scienceExposure, templateWcs, growSize)

    # Clip the BBox to the exposure size
    LLC = templateBBox.getLLC()
    LLC.setX( max(0, LLC.getX()) )
    LLC.setY( max(0, LLC.getY()) ) 

    URC = templateBBox.getURC()
    URC.setX( 10 ) #min(1055, URC.getX()) )
    URC.setY( 10 ) #min(1153, URC.getY()) )

    clippedBBox = afwImage.BBox(LLC, URC)
    print clippedBBox.getX0(), clippedBBox.getY0(), clippedBBox.getX1(), clippedBBox.getY1()
    
    # Get the Template Subexposure
    templateImage       = afwImage.ImageF(templateFile, 0, None, clippedBBox)
    templateMaskedImage = afwImage.MaskedImageF(templateImage)
    mask = templateMaskedImage.getMask()
    mask = 0
    var  = templateMaskedImage.getVariance()
    var  = 0.01
    
    templateWcs         = afwImage.Wcs(templateMetadata)
    templateExposure    = afwImage.ExposureF(templateMaskedImage, templateWcs)

    # Remap to match science exposure
    if remap:
        afwMath.warpExposure(templateExposure,
                             scienceExposure,
                             afwMath.LanczosWarpingKernel(2)
                             )
                             
    
def findOverlapBBox(scienceExposure, templateWcs, growSize):
    """Find the overlap of a Template's Wcs with the extent of a
    Science Exposure"""

    # Get the extent of the science exposure on the sky
    scienceMi  = scienceExposure.getMaskedImage()
    scienceWcs = scienceExposure.getWcs()
    scienceOriginSky = scienceWcs.xyToRaDec(0, 0)
    scienceLimitSky  = scienceWcs.xyToRaDec( scienceMi.getWidth(),
                                             scienceMi.getHeight() )


    # Get the x,y location of the template overlap
    templateOriginOverlap = templateWcs.raDecToXY( scienceOriginSky )
    templateLimitOverlap  = templateWcs.raDecToXY( scienceLimitSky  )

    # Which is the minimium?
    if (templateOriginOverlap.getX() < templateLimitOverlap.getX()):
        lower = templateOriginOverlap
        upper = templateLimitOverlap
    else:
        upper = templateOriginOverlap
        lower = templateLimitOverlap

    # Grow the 2 points
    growL = afwImage.PointI( int(lower.getX() - growSize + 0.5),
                             int(lower.getY() - growSize + 0.5) )

    growU = afwImage.PointI( int(upper.getX() + growSize + 0.5),
                             int(upper.getY() + growSize + 0.5) )


    # Turn into BBox
    templateBBox = afwImage.BBox( growL, growU )

    return templateBBox

if __name__ == '__main__':
    scienceFile        = sys.argv[1]
    scienceImage       = afwImage.ImageF(scienceFile)
    scienceMetadata    = afwImage.readMetadata(scienceFile)
    scienceWcs         = afwImage.Wcs(scienceMetadata)
    
    scienceMaskedImage = afwImage.MaskedImageF(scienceImage)
    mask = scienceMaskedImage.getMask()
    mask = 0
    var  = scienceMaskedImage.getVariance()
    var  = 0.01

    scienceExposure    = afwImage.ExposureF(scienceMaskedImage, scienceWcs)
    
    templateFile       = sys.argv[2]
    getTemplateSubexposure(scienceExposure, templateFile)
