"""

@brief: Implementation of the stage, Dark Currecnt Correction, for the
 nightly Instrument Signature Removal Pipeline

@author: Nicole M. Silvestri
 contact: nms@astro.washington.edu
 file created: Mon Nov 24, 2008  

@file
"""
import eups
import os

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.daf.base as dafBase
import lsst.pex.exceptions as pexEx
import lsst.pex.logging as pexLog
import lsst.pex.policy as pexPolicy
import lsst.ip.isr as ipIsr

def darkCurrentCorrection(chunkExposure, masterChunkExposure, isrPolicy):

    """Dark Current Correction

    @brief The appropriate Master Dark Curent Chunk Exposure is
            retrieved from the Clipboard, scaled, and subtracted from
            the Chunk Exposure to correct for the thermal noise
            contribution of the electronics.

    @return chunkExposure corrected for thermal noise 

    @throw LengthError if chunk and master Exposures are different sizes
    @throw RangeError if chunk and master are derived from different pixels 
    @throw NotFound if any Policy or metadata value can not be obtained
   
    TO DO (as of 11/24/08):
    - Implement sigma-clipping
    - add any additional SDQA statistics requested by SDQA team
    """

    stage = "lsst.ip.isr.darkCurrentCorrection"   
    pexLog.Trace("Entering ISR Stage: ", 4, "%s" % (stage,))

    # Parse the Policy File
    try:
     darkPolicy = isrPolicy.getPolicy("darkPolicy")
     chunkType = isrPolicy.getString("chunkType");
     darkScale = darkPolicy.getDouble("darkScale");
     sigClip = darkPolicy.getBool("sigClip");
     if sigClip = true:
         sigClipVal = darkPolicy.getDouble("sigClipVal");
         # add sigClipping policy info here - not yet implemented
     except pexEx.LsstExceptionStack, e:
         print "Cannot parse the ISR Policy File: %s" % e   
         raise pexExcept.NotFound, "Can not obtain %s Policy Parameters." % (stage,)
   
    chunkMaskedImage = chunkExposure.getMaskedImage()
    chunkMetadata = chunkMaskedImage.getImage().getMetaData()
    masterChunkMaskedImage = masterChunkExposure.getMaskedImage()
    masterChunkMetadata = chunkMaskedImage.getImage().getMetaData()

    # Check that the Master Bias Chunk Exposure and Chunk Exposure are
    # the same size.

    numCols = chunkMaskedImage.getCols()
    numRows = chunkMaskedImage.getRows() 

    mnumCols = masterChunkMaskedImage.getCols()
    mnumRows = masterChunkMaskedImage.getRows() 

    if numCols != mnumCols || numRows != mnumRows:
        raise pexExcept.LengthError, "In %s: Chunk Exposure and Master Bias Chunk Exposure are not the same size." % (stage,)

    # Check that the Master Bias Chunk Exposure and Chunk Exposure are
    # derived from the same pixels (eg. both are from the same amp,
    # CCD, or raft).

    if chunkType == "amp":
        
        ampidField = chunkMetadata.findUnique("AMPID")
        mampidField = masterChunkMetadata.findUnique("AMPID");
        if ampField:
            ampid = ampidField.getValue()
            mampid = mampidField.getValue()
        else:
            raise pexExcept.NotFound, "In %s: Could not get AMPID from the Metadata." %(stage,)

    if chunkType == "ccd":
        
        ccdidField = chunkMetadata.findUnique("CCDID")
        mccdidField = masterChunkMetadata.findUnique("CCDID");
        if ccdField:
            ccdid = ccdidField.getValue()
            mccdid = mccdidField.getValue()
        else:
            raise pexExcept.NotFound, "In %s: Could not get CCDID from the Metadata." %(stage,)

    if chunkType == "raft":
        
        raftidField = chunkMetadata.findUnique("RAFTID")
        mraftidField = masterChunkMetadata.findUnique("RAFTID");
        if raftField:
            raftid = raftidField.getValue()
            mraftid = mraftidField.getValue()
        else:
            raise pexExcept.NotFound, "In %s: Could not get RAFTID from the Metadata." %(stage,)
        
    if ampid != mampid || ccdid != mccdid || raftid != mraftid:
        raise pexExcept.RangeError, "In %s: Chunk Exposure and Master Bias Chunk Exposure are not derived from the same pixels." % (stage,)
        }

    # Get the relevant metadata 

    fileNameField = masterChunkMetadata.findUnique("FILENAME")
    if fileNameField:
        fileName = fileNameField.getValue()
    else:
        raise pexExcept.NotFound,"In %s: Could not get FILENAME from the Metadata." %(stage,) 

    meanDarkField = masterChunkMetadata.findUnique("MEAN");
    if meanDarkField:
        meanDark = meanDarkField.getValue()
    else:
        raise pexExcept.NotFound,"In %s: Could not get MEAN from the master Metadata." %(stage,)

    exptimeField = chunkMetadata.findUnique("EXPTIME")
    mexptimeField = masterChunkMetadata.findUnique("EXPTIME")
    if exptimeField:
        exptime = exptimeField.getValue()
        mexptime = mexptimeField.getValue()
    else:
        raise pexExcept.NotFound, "In %s: Could not get EXPTIME from Chunk Metadata." % (stage,)

    scale= exptime/mexptime 
    if exptime != mexptime:
        masterChunkMaskedImage *= scale 

    # subtract the master Bias MaskedImage form the Chunk MaskedImage.
    # Allow for additional scaling if desired.

    if darkScale:
        masterChunkMaskedImage *= darkScale
        chunkMaskedImage -= masterChunkMaskedImage
    else:
        chunkMaskedImage -= masterChunkMaskedImage
        
    # Record final stage provenance to the Image Metadata
    dateTime = dafBase.DateTime.utc2mjd()
    chunkMetadata.addProperty(dafBase.Dataproperty("DARK_MJD", dateTime))
    chunkMetadata.addProperty(dafBase.DataProperty("DARK_MC", fileName))
    chunkMetadata.addProperty(dafBase.DataProperty("DARK_MU", meanDark))
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
                              
