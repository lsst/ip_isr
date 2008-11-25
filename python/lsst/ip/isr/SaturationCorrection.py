"""

@brief: Implementation of the stage, Saturation Correction, for the
 nightly Instrument Signature Removal Pipeline

@author: Nicole M. Silvestri
 contact: nms@astro.washington.edu
 file created: 

@file
"""

import lsst.afw.image as afwImage
import lsst.daf.base as dafBase
import lsst.detection as det
import lsst.pex.exceptions as pexEx
import lsst.pex.logging as pexLog
import lsst.pex.policy as pexPolicy
import lsst.ip.isr as ipIsr

def saturationCorrection(chunkExposure, isrPolicy, lookupTable):

    """Saturation Correction for the Chunk Exposure

    Detect and mask pixels that are saturated in the A/D converter or have
    excessive non-linear response in the Chunk Exposure.  Uses the
    lsst::detection::DetectionSet to create footprints of saturated pixels
    above a threshold as given in the ISR Policy.

    Grow by additional pixels (as given in the ISR Policy) to mask charge
    spillover.  Set appropriate bits in the Mask and interpolate over
    masked pixels using the 'InterpolateOverMaskedPixels' utility
    function.

    @return chunkExposure with saturated pixels masked and interpolated
  
    @throw NotFound if any metadata or policy parameter can not be obtained
   
    TO DO (as of Tue 11/20/08):
    - delineate between A/D saturated pixels and other?
    - Calculate additional SDQA metrics as requested by SDQA team
    - inplement option for using satLimit from LookupTable??
    - sigma clipping?
    """

    stage = "lsst.ip.isr.saturationCorrection"   
    pexLog.Trace("Entering ISR Stage: ", 4, "%s" % (stage,))

    # Parse the Policy File
    try:
        pexLog.Trace("%s" % (stage,), 4, "Parsing the ISR Policy File.")
        satPolicy = isrPolicy.getPolicy("saturationPolicy")
        
        satTable = satPolicy.getBool("satTable")
        useDefSat = satPolicy.getBool("useDefSat")
        satGrow = satPolicy.getInt("grow")
        satThresh = satPolicy.getInt("threshold")
        #nSatPixMin = satPolicy.getInt("nSatPixMin")
    except pexEx.LsstExceptionStack, e:
        print "Cannot parse the ISR Policy File: %s" % e   
    
    chunkMaskedImage = chunkExposure.getMaskedImage()
    chunkMetadata = chunkMaskedImage.getImage().getMetaData()

    try:
        pexLog.Trace("%s" % (stage,), 4, "Obtaining additional parameters from chunkMetadata.")
        satField = chunkMetadata.findUnique("SATURATE")
        satLimit = satField.getValue()
        print "satLimit: %s" % (satLimit,)
    except pexEx.LsstExceptionStack, e:
        print "Cannot get the saturation limit from the chunkExposure: %s" % e      

    chunkMask = chunkMaskedImage.getMask()
    satMaskBit = chunkMask.getPlaneBitMask("SAT")

    # QQQ: Can we distinguish between pixels saturated in the A/D
    # converter and those just saurated on chip?  If so, we don't want
    # to grow around the A/D saturated pixels (in unbinned data).

    pexLog.Trace("%s" % (stage,), 4, "Finding saturated footprints.")
    satFootprintSet = det.DetectionSetD(chunkMaskedImage, det.Threshold(satThresh))
    satFootprintList = satFootprintSet.getFootprints()
    numSatFootprints = len(satFootprintList)
    print "Found %s saturated footprints." % (numSatFootprints,)

    numSatPix = 0
    for feet in satFootprintList:
        numSatPix = numSatPix + feet.getNpix()

    print "Found %s saturated pixels." % (numSatPix,)
    
    pexLog.Trace("%s" % (stage,), 4, "Growing around all saturated footprints.")
    grownSatFootprintSet = det.DetectionSetD(satFootprintSet, satGrow)
    grownSatFootprintList = grownSatFootprintSet.getFootprints()
    numGrownSatFootprints = len(grownSatFootprintList)
    print "Found %s grown saturated footprints." % (numGrownSatFootprints,)

    numGrownSatPix = 0
    for bigFeet in grownSatFootprintList:
        numGrownSatPix = numGrownSatPix + bigFeet.getNpix()

    print "Found %s grown saturated pixels." % (numGrownSatPix,)
    
    # Mask all of the grown saturated pixel footprints.  Using "SAT"
    # bitmask for all pixels in the footprint.

    # QQQ: Do we want to distinguish between grown pixels and those
    # that were actually saturated?  What bitmask would that be
    # ("GROW")?
    
    det.setMaskFromFootprintList(chunkMask, grownSatFootprintList, satMaskBit)

    pexLog.Trace("%s" % (stage,), 4, "Interpolating over all saturated footprints.")
    ipIsr.interpolateOverMaskedPixels(chunkExposure, isrPolicy)
    
    # Record the stage provenance to the Image Metadata (for now) this
    # will eventually go into an ISR_info object.

    pexLog.Trace("%s" % (stage,), 4, "Recording ISR provenance information.")
    chunkMetadata.addProperty(dafBase.DataProperty("SATU_TH", satThresh))
    chunkMetadata.addProperty(dafBase.DataProperty("SATU_PIX", numSatPix))    
    chunkMetadata.addProperty(dafBase.DataProperty("SATU_FP", numSatFootprints))
    chunkMetadata.addProperty(dafBase.DataProperty("SATU_END", "Completed Successfully")) 
    chunkMaskedImage.setMetadata(chunkMetadata)

    # Calculate any additional SDQA Metrics and write all metrics to
    # the SDQA object (or directly to the clipboard)
                               
    # pexLog.Trace("%s" % (stage,), 4, "Recording SDQA metric information." )
                              
    """ Return the following for SDQA:
    - numSatFootprints
    - numSatPix
    - numGrownSatFootprints
    - numGrownSatPix                          
    """                       
    pexLog.Trace("%s" % (stage,), 4, "Completed Successfully" )
    pexLog.Trace("Leaving ISR Stage: ", 4, "%s" % (stage,))
                              
