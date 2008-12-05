"""
@brief: Implementation of the stage, Fringe Correction, for the
 nightly Instrument Signature Removal Pipeline

@author: Nicole M. Silvestri
 contact: nms@astro.washington.edu
 file created: Tue Nov 25, 2008  

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
   
def fringeCorrection(chunkExposure, masterChunkExposure, isrPolicy):

    """Fringe Correction

    @brief Remove fringe pattern associated with night sky emission lines.  The intensity of the night sky lines varies with time.  Additional 
scaling is necessary.  Fringing varies with filter/bandpass so master Fringe and chunk Exposure should be taken through the same filter.
    
    @return chunkExposure corrected for fringe

    @throw LengthError if chunk and master Exposures are different sizes
    @throw RangeError if chunk and master are derived from different pixels 
    @throw NotFound if any Policy or metadata value can not be obtained
   
    TO DO (as of 11/24/08):
    - deal with fringe scaling - this needs to be done using inter-slice communication
    """

    stage = "lsst.ip.isr.fringeCorrection"   
    pexLog.Trace("Entering ISR Stage: ", 4, "%s" % (stage,))

    # Parse the Policy File
    try:
        fringePolicy = isrPolicy.getPolicy("fringePolicy")
        chunkType = isrPolicy.getString("chunkType")
        sigClip = fringePolicy.getBool("sigClip")
        sigClipVal = fringePolicy.getDouble("sigClipVal")
        # THIS IS DERIVED FROM SEVERAL PLACES ON THE MASTER
        # FRINGE. THIS NEEDS MORE WORK THAN JUST GRABBING A SCALE
        # FACTOR
        fringeScale = fringePolicy.getDouble("fringeScale")
     except pexEx.LsstExceptionStack, e:
         print "Cannot parse the ISR Policy File: %s" % e   
         raise pexExcept.NotFound, "Can not obtain Fringe Correction Policy Parameters."
   
    chunkMaskedImage = chunkExposure.getMaskedImage()
    chunkMetadata = chunkMaskedImage.getImage().getMetaData()
    masterChunkMaskedImage = masterChunkExposure.getMaskedImage()
    masterChunkMetadata = chunkMaskedImage.getImage().getMetaData()

    # Check that the Master Fringe Chunk Exposure and Chunk Exposure are
    # the same size.

    numCols = chunkMaskedImage.getCols()
    numRows = chunkMaskedImage.getRows() 

    mnumCols = masterChunkMaskedImage.getCols()
    mnumRows = masterChunkMaskedImage.getRows() 

    if numCols != mnumCols || numRows != mnumRows:
        raise pexExcept.LengthError, "In %s: Chunk Exposure and Master Fringe Chunk Exposure are not the same size." % (stage,)

    # Check that the Master Fringe Chunk Exposure and Chunk Exposure are
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
        raise pexExcept.RangeError, "In %s: Chunk Exposure and Master Fringe Chunk Exposure are not derived from the same pixels." % (stage,)
        }

    # Get the relevant metadata 

    fileNameField = masterChunkMetadata.findUnique("FILENAME")
    if fileNameField:
        fileName = fileNameField.getValue()
    else:
        raise pexExcept.NotFound,"In %s: Could not get FILENAME from the Metadata." %(stage,) 

    meanFlatField = masterChunkMetadata.findUnique("MEAN");
    if meanFlatField:
        meanFlat = meanFlatField.getValue()
    else:
        raise pexExcept.NotFound,"In %s: Could not get MEAN from the master Metadata." %(stage,)

    flterField = chnkMetadata.findUnique("FILTER")
    mfilterField = masterChunkMetadata.findUnique("FILTER")
    if filterField:
        filter = filterField.getValue()
        mfilter = mfilterField.getValue()
    else:
        raise pexExcept.NotFound,"In %s: Could not get FILTER from the Metadata." %(stage,)

    # Make sure the chunk and master have the same filter designation 
    if filter != mfilter:
        raise pexExcept.NotFound,"In %s, Chunk Exposure and Master Fringe Chunk Exposure are not from the same FILTER." % (stage,) 

    #  Divide the Chunk Exposure by the scaled Fringe Chunk
    #  Exposure. Allow for additional scaling if requested.

    if fringeScale:
        masterChunkMaskedImage *= fringeScale
        chunkMaskedImage /= masterChunkMaskedImage
    else:
        chunkMaskedImage /= masterChunkMaskedImage
        
    # Record final stage provenance to the Image Metadata
    dateTime = dafBase .DateTime.utc2mjd()
    chunkMetadata.addProperty(dafBase.Dataproperty("FRIN_MJD", dateTime))
    chunkMetadata.addProperty(dafBase.DataProperty("FRIN_MC", fileName))
    chunkMetadata.addProperty(dafBase.DataProperty("FRIN_END", "Completed Successfully")) 
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
                              
