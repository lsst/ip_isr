"""

@brief: Implementation of the stage, Bias Correction, for the
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

dataDir = eups.productDir("afwdata")
    if not dataDir:
        raise RuntimeError("Must set up afwdata to run this program.")
   
def biasCorrection(chunkExposure, masterChunkExposure, isrPolicy):

    """Bias Correction

    @brief The appropriate Master Bias Chunk Exposure is retrieved
    from the Clipboard and is subtracted from the Chunk Exposure to
    correct for structure in the bias (offsets from the overscan
    level).

    @return chunkExposure corrected for bias

    @throw LengthError if chunk and master Exposures are different sizes
    @throw RangeError if chunk and master are derived from different pixels 
    @throw NotFound if any Policy or metadata value can not be obtained
   
    TO DO (as of 11/24/08):
    - Implement sigma-clipping
    - implement verification that images are derived from same raft if chunk=raft?
    - add any additional SDQA statistics requested by SDQA team
    """

    stage = "lsst.ip.isr.biasCorrection"   
    pexLog.Trace("Entering ISR Stage: ", 4, "%s" % (stage,))

    # Parse the Policy File
    try:
     biasPolicy = isrPolicy.getPolicy("biasPolicy")
     chunkType = biasPolicy.getString("chunkType");
     biasScale = biasPolicy.getDouble("biasScale");
     sigClip = biasPolicy.getBool("sigClip");
     if sigClip = true:
         sigClipVal = biasPolicy.getDouble("sigClipVal");
         # add sigClipping policy info here - not yet implemented
     except pexEx.LsstExceptionStack, e:
         print "Cannot parse the ISR Policy File: %s" % e   
         raise pexExcept.NotFound, "Can not obtain Bias Correction Policy Parameters."
   
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

    meanBiasField = masterChunkMetadata.findUnique("MEAN");
    if meanBiasField:
        meanBias = meanBiasField.getValue()
    else:
        raise pexExcept.NotFound,"In %s: Could not get MEAN from the master Metadata." %(stage,) 

    # subtract the master Bias MaskedImage form the Chunk MaskedImage.
    # Allow for aditional scaling if desired.

    if biasScale:
        masterChunkMaskedImage *= biasScale
        chunkMaskedImage -= masterChunkMaskedImage
    else:
        chunkMaskedImage -= masterChunkMaskedImage
        
    # Record final stage provenance to the Image Metadata
    dateTime = dafBase.DateTime.utc2mjd()
    chunkMetadata.addProperty(dafBase.Dataproperty("BIAS_MJD", dateTime))
    chunkMetadata.addProperty(dafBase.DataProperty("BIAS_MC", fileName))
    chunkMetadata.addProperty(dafBase.DataProperty("BIAS_MU", meanBias))
    chunkMetadata.addProperty(dafBase.DataProperty("BIAS_END", "Completed Successfully")) 
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
                              
