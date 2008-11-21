"""

@brief: Implementation of the stage, Overscan Correct and Trim, for the
 nightly Instrument Signature Removal Pipeline

@author: Nicole M. Silvestri
 contact: nms@astro.washington.edu
 file created: Thu Nov 20, 2008  

@file
"""

import lsst.afw.image as afwImage
import lsst.daf.base as dafBase
import lsst.detection as det
import lsst.pex.exceptions as pexEx
import lsst.pex.logging as pexLog
import lsst.pex.policy as pexPolicy
import lsst.ip.isr as ipIsr

def overscanCorrectAndTrim(chunkExposure, isrPolicy):

    """Overscan Correct and Trim
    
    Subtract a single constant (using all input pixels for the
    statistic) or one-dimensional function from the Chunk Exposure
    which varies along the overscan.  Use spline or polynomial fit
    along the overscan region (want a slowly varying function).  For
    the 1-D representation, Pan-STARRS determines the input values to
    the fit as representations on the coordinate along the overscan,
    with the statistic derived from the pixels in the perpendicular
    direction at each location.  Sigma-clipping on input data is a
    necessary option.
  
    The Image is also trimmed to remove the overscan region and any
    unexposed pixels.
  
    @return chunkExposure corrected for overscan and trimmed of all
    non-illuminated pixels
  
    @throw Runtime if sage has already been run on the image
    @throw Runtime if any Policy or metadata value can not be obtained
    @throw InvalidParameter if functional form for the overscan fit is invalid 
   
    TO DO (as of 11/20/08):
    - Implement sigma-clipping
    - Trim ramp-up on overscan (do not include in fits)
    - Implement additional constant modes and spline fit
    - add a smoothing option to the functional fit? 
    - add any additional SDQA statistics requested by SDQA team
    """

    satStage = "lsst.ip.isr.overscanCorrectAndTrim"   
    pexLog.Trace("Entering ISR Stage: ", 4, "%s" % (satStage,))

    # Parse the Policy File
    try:
     overscanPolicy = isrPolicy.getPolicy("overscanPolicy")
     method = overPolicy.getString("method")
     smooth = overPolicy.getBool("smooth")
     sigClip = overPolicy.getBool("sigClip")
     datasecKey = overPolicy.getString("datasec")
     biassecKey = overPolicy.getString("biassec")
     trimsecKey = overPolicy.getString("trimsec")
    except pexEx.LsstExceptionStack, e:
        print "Cannot parse the ISR Policy File: %s" % e   

    try:
        if method == "function":
            funcForm = overPolicy.getString("funcForm")
            funcOrder = overPolicy.getInt("funcOrder")
            stepSize = overPolicy.getDouble("stepSize")
            numParams = overPolicy.getDouble("numParams")
        else: 
            constantMeth = overPolicy.getString("constantMeth")
    except pexEx.LsstExceptionStack, e:
        print "Cannot parse the Overscan Policy File for 'method': %s" % e

     try:
         if smooth = true:
             smoothOpt = overPolicy.getString("smoothOpt")
             if smoothOpt == "gaussian":
                 smoothSigma = overPolicy.getInt("smoothSigma")
             else:
                 smoothOrder = overPolicy.getInt("smoothOrder")
     except pexEx.LsstExceptionStack, e:
         print "Cannot parse the Overscan Policy File for 'smooth': %s" % e

     try:
         if sigClip = true:
             sigClipVal = overPolicy.getInt("sigClipVal") 
     except pexEx.LsstExceptionStack, e:
         print "Cannot parse the Overscan Policy File for 'sigClip': %s" % e
    
    chunkMaskedImage = chunkExposure.getMaskedImage()
    chunkMetadata = chunkMaskedImage.getImage().getMetaData()

    try:
        pexLog.Trace("%s" % (satStage,), 4, "Obtaining additional parameters from chunkMetadata.")
        datasecProp = chunkMetadata.findUnique(datasecKey)
        datasec = datasecProp.getValue()
        print "DATASEC: ", datasec
        trimsecProp = chunkMetadata.findUnique(trimsecKey)
    except pexEx.LsstExceptionStack, e:
        print "Cannot get the requested metadata from the chunkExposure: %s" % e      

    chunkMask = chunkMaskedImage.getMask()
    satMaskBit = chunkMask.getPlaneBitMask("SAT")

    # QQQ: Can we distinguish between pixels saturated in the A/D
    # converter and those just saurated on chip?  If so, we don't want
    # to grow around the A/D saturated pixels (in unbinned data).

    pexLog.Trace("%s" % (satStage,), 4, "Finding saturated footprints.")
    satFootprintSet = det.DetectionSetD(chunkMaskedImage, det.Threshold(satThresh))
    satFootprintList = satFootprintSet.getFootprints()
    numSatFootprints = len(satFootprintList)
    print "Found %s saturated footprints." % (numSatFootprints,)

    numSatPix = 0
    for feet in satFootprintList:
        numSatPix = numSatPix + feet.getNpix()

    print "Found %s saturated pixels." % (numSatPix,)
    
    pexLog.Trace("%s" % (satStage,), 4, "Growing around all saturated footprints.")
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

    pexLog.Trace("%s" % (satStage,), 4, "Interpolating over all saturated footprints.")
    ipIsr.interpolateOverMaskedPixels(chunkExposure, isrPolicy)
    
    # Record the stage provenance to the Image Metadata (for now) this
    # will eventually go into an ISR_info object.

    pexLog.Trace("%s" % (satStage,), 4, "Recording ISR provenance information.")
    chunkMetadata.addProperty(dafBase.DataProperty("SATU_TH", satThresh))
    chunkMetadata.addProperty(dafBase.DataProperty("SATU_PIX", numSatPix))    
    chunkMetadata.addProperty(dafBase.DataProperty("SATU_FP", numSatFootprints))
    chunkMetadata.addProperty(dafBase.DataProperty("SATU_END", "Completed Successfully")) 
    chunkMaskedImage.setMetadata(chunkMetadata)

    # Calculate any additional SDQA Metrics and write all metrics to
    # the SDQA object (or directly to the clipboard)
                               
    # pexLog.Trace("%s" % (satStage,), 4, "Recording SDQA metric information." )
                              
    """ Return the following for SDQA:
    - numSatFootprints
    - numSatPix
    - numGrownSatFootprints
    - numGrownSatPix                          
    """                       
    pexLog.Trace("%s" % (satStage,), 4, "Completed Successfully" )
    pexLog.Trace("Leaving ISR Stage: ", 4, "%s" % (satStage,))
                              
