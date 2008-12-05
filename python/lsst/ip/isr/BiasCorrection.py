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
    pexLog.Trace("%s" % (stage,), 4, "Entering ISR Bias Correction stage." )

    # Parse the Policy File
    pexLog.Trace("%s" % (stage,), 4, "Parsing the ISR Policy File." )
    try:
        biasPolicy = isrPolicy.getPolicy("biasPolicy")
        chunkType = isrPolicy.getString("chunkType");
        biasScale = biasPolicy.getDouble("biasScale");
        sigClip = biasPolicy.getBool("sigClip");
        if sigClip == "true":
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

    pexLog.Trace("%s" % (stage,), 4, "Verifying Master and Chunk Exposures are the same size." )
    numCols = chunkMaskedImage.getCols()
    numRows = chunkMaskedImage.getRows() 

    mnumCols = masterChunkMaskedImage.getCols()
    mnumRows = masterChunkMaskedImage.getRows() 

    if numCols != mnumCols or numRows != mnumRows:
        raise pexExcept.LengthError, "In %s: Chunk Exposure and Master Bias Chunk Exposure are not the same size." % (stage,)
    else:
        pexLog.Trace("%s" % (stage,), 4, "Success: Master and Chunk are the same size." )

    # Check that the Master Bias Chunk Exposure and Chunk Exposure are
    # derived from the same pixels (eg. both are from the same amp,
    # CCD, or raft).

    pexLog.Trace("%s" % (stage,), 4, "Verifying Master and Chunk Exposures are derived from the same pixels." )
    if chunkType == "amp":

        ampidField = chunkMetadata.findUnique("AMPLIST")
        mampidField = masterChunkMetadata.findUnique("AMPLIST");
        if ampidField:
            ampid = ampidField.getValueString()
            mampid = mampidField.getValueString()
        else:
            raise pexExcept.NotFound, "In %s: Could not get AMPID from the Metadata." %(stage,)

        if ampid != mampid:
            raise pexExcept.RangeError, "In %s: Chunk Exposure and Master Bias Chunk Exposure are not derived from the same pixels." % (stage,)
        else:
            pexLog.Trace("%s" % (stage,), 4, "Success: Master and Chunk Exposures are derived from the same pixels." )

    if chunkType == "ccd":
        
        ccdidField = chunkMetadata.findUnique("CCDID")
        mccdidField = masterChunkMetadata.findUnique("CCDID");
        if ccdField:
            ccdid = ccdidField.getValueString()
            mccdid = mccdidField.getValueString()
        else:
            raise pexExcept.NotFound, "In %s: Could not get CCDID from the Metadata." %(stage,)

        if ccdid != mccdid:
            raise pexExcept.RangeError, "In %s: Chunk Exposure and Master Bias Chunk Exposure are not derived from the same pixels." % (stage,)
        else:
            pexLog.Trace("%s" % (stage,), 4, "Success: Master and Chunk Exposures are derived from the same pixels." )

    if chunkType == "raft":
        
        raftidField = chunkMetadata.findUnique("RAFTID")
        mraftidField = masterChunkMetadata.findUnique("RAFTID");
        if raftField:
            raftid = raftidField.getValueString()
            mraftid = mraftidField.getValueString()
        else:
            raise pexExcept.NotFound, "In %s: Could not get RAFTID from the Metadata." %(stage,)
        
        if raftid != mraftid:
            raise pexExcept.RangeError, "In %s: Chunk Exposure and Master Bias Chunk Exposure are not derived from the same pixels." % (stage,)
        else:
            pexLog.Trace("%s" % (stage,), 4, "Success: Master and Chunk Exposures are derived from the same pixels." )

    # Get the relevant metadata 

    pexLog.Trace("%s" % (stage,), 4, "Obtaining additional parameters from the metadata." )
    fileNameField = masterChunkMetadata.findUnique("FILENAME")
    if fileNameField:
        fileName = fileNameField.getValueString()
        print "Filename: ", fileName
    else:
        raise pexExcept.NotFound,"In %s: Could not get FILENAME from the Metadata." %(stage,) 

#    meanBiasField = masterChunkMetadata.findUnique("MEAN");
#    if meanBiasField:
#        meanBias = meanBiasField.getValueDouble()
#        print "Mean of Master Bias: ", meanBias
#    else:
#        raise pexExcept.NotFound,"In %s: Could not get MEAN from the master Metadata." %(stage,) 

    # subtract the master Bias MaskedImage form the Chunk MaskedImage.
    # Allow for aditional scaling if desired.

    pexLog.Trace("%s" % (stage,), 4, "Subtracting Master Bias from Chunk Exposure." )
    if biasScale:
        masterChunkMaskedImage *= biasScale
        chunkMaskedImage -= masterChunkMaskedImage
    else:
        chunkMaskedImage -= masterChunkMaskedImage
        
    # Record final stage provenance to the Image Metadata
    pexLog.Trace("%s" % (stage,), 4, "Recording final provenance information." )
#    dateTime = dafBase.DateTime.utc2mjd()
 # chunkMetadata.addProperty(dafBase.Dataproperty("BIAS_MJD", dateTime))
    chunkMetadata.addProperty(dafBase.DataProperty("BIAS_MC", fileName))
 #   chunkMetadata.addProperty(dafBase.DataProperty("BIAS_MU", meanBias))
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
                              
