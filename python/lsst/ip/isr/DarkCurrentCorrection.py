"""

@brief Implementation of the stage, Dark Currecnt Correction, for the
       nightly Instrument Signature Removal Pipeline

@author Nicole M. Silvestri
        University of Washington
        nms@astro.washington.edu
        created: Mon Nov 24, 2008  

@file
"""
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.daf.base as dafBase
import lsst.pex.exceptions as pexEx
import lsst.pex.logging as pexLog
import lsst.pex.policy as pexPolicy
import lsst.ip.isr as ipIsr

def darkCurrentCorrection(chunkExposure, masterExposure, isrPolicy):

    """Dark Current Correction

    @brief The appropriate Master Dark Curent Exposure is
            retrieved from the Clipboard, scaled, and subtracted from
            the Chunk Exposure to correct for the thermal noise
            contribution of the electronics.

    @return chunkExposure corrected for thermal noise 

    @throw LengthError if chunk and master Exposures are different sizes
    @throw RangeError if chunk and master Exposures are derived from different pixels 
    @throw NotFound if any Policy or metadata value can not be obtained
   
    TO DO (as of 1/08/09):
    - add any additional SDQA statistics requested by SDQA team
    - scale darktime properly
    """

    stage = "lsst.ip.isr.darkCurrentCorrection"   
    pexLog.Trace("%s" % (stage,), 4,"Entering ISR Dark Correction Stage" )

    # Parse the Policy File
    pexLog.Trace("%s" % (stage,), 4, "Parsing the ISR Policy File." )
    try:
        darkPolicy = isrPolicy.getPolicy("darkPolicy")
        chunkType = isrPolicy.getString("chunkType")
        darkScale = darkPolicy.getDouble("darkScale")
        chunkField = darkPolicy.getString("chunkField")
        fileNameField = darkPolicy.getString("fileName")
        expField = darkPolicy.getString("exptimeField")
        darkField = darkPolicy.getString("darktimeField")
        sigClip = darkPolicy.getBool("sigClip")
        if sigClip == True:
            sigClipVal = darkPolicy.getDouble("sigClipVal")
            # add sigClipping policy info here - not yet implemented
    except pexEx.LsstExceptionStack, e:
        pexLog.Trace("%s" % (stage,), 4, "Can not parse the ISR Policy File." % (e,))   
        raise pexExcept.NotFound, "In %s: Can not obtain policy parameters from the ISR Policy File." % (stage,)
   
    chunkMaskedImage = chunkExposure.getMaskedImage()
    chunkMetadata = chunkMaskedImage.getImage().getMetaData()
    masterMaskedImage = masterExposure.getMaskedImage()
    masterMetadata = chunkMaskedImage.getImage().getMetaData()

    # Check that the Master Bias Chunk Exposure and Chunk Exposure are
    # the same size.
    
    pexLog.Trace("%s" % (stage,), 4, "Verifying Master and Chunk Exposures are the same size." )
    numWidth = chunkMaskedImage.getWidth()
    numHeight = chunkMaskedImage.getHeight() 
    pexLog.Trace("%s" % (stage,), 4, "Chunk Exposure NumWidth, NumHeight: %s, %s" % numWidth, numHeight )

    mnumWidth = masterMaskedImage.getWidth()
    mnumHeight = masterMaskedImage.getHeight() 
    pexLog.Trace("%s" % (stage,), 4, "Master Exposure NumWidth, NumHeight: %s, %s" % mnumWidth, mnumHeight )

    if (numWidth != mnumWidth) or (numHeight != mnumHeight):
        raise pexExcept.LengthError, "In %s: Chunk Exposure and Master Dark Current Exposure are not the same size." % (stage,)
    else:
        pexLog.Trace("%s" % (stage,), 4, "Success: Master and Chunk Exposures are the same size." )

    # Check that the Master Bias Chunk Exposure and Chunk Exposure are
    # derived from the same pixels (eg. both are from the same 'chunk').

    pexLog.Trace("%s" % (stage,), 4, "Verifying Master and Chunk Exposures are derived from the same pixels." )
        
    chunkField = chunkMetadata.findUnique(chunkField)
    mchunkField = masterMetadata.findUnique(chunkField);
    if chunkField and mampidField:
        chunkId = chunkField.getValueString()
        print "ID chunk: ", chunkId
        mchunkId = mchunkField.getValueString()
        print "ID master: ", mchunkId
    else:
        raise pexExcept.NotFound, "In %s: Could not get %s from the Metadata." % (stage, chunkType)

    if chunkId != mchunkId:
        raise pexExcept.RangeError, "In %s: Chunk Exposure and Master Dark Current Exposure are not derived from the same pixels." % (stage,)
    else:
        pexLog.Trace("%s" % (stage,), 4, "Success: Master and Chunk Exposures are derived from the same pixels." )

    # Get additional metadata
    
    pexLog.Trace("%s" % (stage,), 4, "Obtaining additional parameters from the metadata." )
    fileNameKey = masterMetadata.findUnique(fileNameField)
    if fileNameKey:
        fileName = fileNameKey.getValueString()
        pexLog.Trace("%s" % (stage,), 4, "Master Dark Current Filename: %s" % (fileName,))
    else:
        raise pexExcept.NotFound,"In %s: Could not get FILENAME from the Metadata." %(stage,) 

#    meanDarkField = masterMetadata.findUnique("MEAN");
#    if meanDarkField:
#        meanDark = meanDarkField.getValue()
#    else:
#        raise pexExcept.NotFound,"In %s: Could not get MEAN from the master Metadata." %(stage,)

    exptimeField = chunkMetadata.findUnique(expField)
    darkTimeField = chunkMetadata.findUnique(darkField)
    mexptimeField = masterMetadata.findUnique(expField)
    mdarkTimeField = masterMetadata.findUnique(darkField)
    if exptimeField and mexptimeField and darkTimeField and mdarkTimeField:
        exptime = exptimeField.getValueDouble()
        pexLog.Trace("%s" % (stage,), 4,"EXPTIME for Chunk: %s" % (exptime,))
        mexptime = mexptimeField.getValueDouble()
        pexLog.Trace("%s" % (stage,), 4,"EXPTIME for Master Dark: %s" % (mexptime,))
        darkTime = darkTimeField.getValueDouble()
        pexLog.Trace("%s" % (stage,), 4,"DARKTIME for Chunk: %s" % (darktime,))
        mdarkTime = mdarkTimeField.getValueDouble()
        pexLog.Trace("%s" % (stage,), 4,"DARKTIME for Master Dark: %s" % (mdarktime,))
    else:
        raise pexExcept.NotFound, "In %s: Could not get EXPTIME or DARKTIME from Chunk Metadata." % (stage,)

    pexLog.Trace("%s" % (stage,), 4, "Normalizing Exposures to DARKTIME." )
  
  # The darktime scaling needs some work - need to normalize the dark
  # times by fitting polynomial to the exposure times for surrounding
  # chunks (?) as PS does it?  Is this because of the shutter delay
  # (dark and exptimes will be different by fractions of a second) 
  #      if darktime and mdarkTime:  
  #      scale = exptime/darkTime
  #      masterMaskedImage *= scale
    pexLog.Trace("%s" % (stage,), 4, "Scaling Master Dark Exposure to Chunk Exposure's EXPTIME." )
    if exptime != mexptime:
        scale = exptime/mexptime 
        masterMaskedImage *= scale 
    else:
        # What is this supposed to do? Early return? It's not in a loop
        continue

    # subtract the Master Dark Current MaskedImage form the Chunk
    # MaskedImage.  Allow for additional scaling if desired.
    
    pexLog.Trace("%s" % (stage,), 4, "Subtracting Master Dark from Chunk Exposure." )
    if darkScale:
        pexLog.Trace("%s" % (stage,), 4, "Additional scale factor applied to Master Dark: %s" %(darkScale,))
        masterMaskedImage *= darkScale
        chunkMaskedImage -= masterMaskedImage
    else:
        chunkMaskedImage -= masterMaskedImage
        
    # Record final stage provenance to the Image Metadata
    pexLog.Trace("%s" % (stage,), 4, "Recording final provenance information." )
    #dateTime = dafBase.DateTime.utc2mjd()
    #chunkMetadata.addProperty(dafBase.Dataproperty("DARK_MJD", dateTime))
    chunkMetadata.addProperty(dafBase.DataProperty("DARK_MC", fileName))
    #chunkMetadata.addProperty(dafBase.DataProperty("DARK_MU", meanDark))
    chunkMetadata.addProperty(dafBase.DataProperty("DARK_SCA", scale))
    chunkMetadata.addProperty(dafBase.DataProperty("DARK_END", "Completed Successfully")) 
    chunkMaskedImage.setMetadata(chunkMetadata);
  
    # Calculate any additional SDQA Metrics and write all metrics to
    # the SDQA object (or directly to the clipboard)
                               
    pexLog.Trace("%s" % (stage,), 4, "Recording SDQA metric information." )
                              
    """ Return the following for SDQA:
    - ?
    - ?
    - ?
    """
    
    pexLog.Trace("%s" % (stage,), 4, "Completed Successfully" )
    pexLog.Trace("Leaving ISR Stage: ", 4, "%s" % (stage,))
