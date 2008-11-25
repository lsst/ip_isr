"""
@brief: Implementation of the stage, Flat Field Correction, for the
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
import lsst.ip.isr.IlluminationCorrection as ipIsrIllum


dataDir = eups.productDir("afwdata")
    if not dataDir:
        raise RuntimeError("Must set up afwdata to run this program.")
   
def flatFieldCorrection(chunkExposure, masterChunkExposure, isrPolicy):

    """Flat Field Correction

    @brief Divide the Chunk Exposure by the (normalized) Master Flat Field Chunk
    Exposure(s) to correct for pixel-to-pixel variations (eg. optics, vignetting,
    thickness variations, gain, etc.). The Master Flat Field Chunk Exposure can
    be one of potentially three different types of flats (dome, twilight, night
    sky) with further sub-divisions into LSST filters (ugrizy) or bandpasses.
   
    Calls 'Illumination Correction' stage (DR or nightly) to perform the
    illumination correction.
   
    Dome Flats: correct for the pixel-to-pixel variations in the response og the
    CCD.  These will be the 'Stubb's' tunable laser flats.
   
    Twilight Flats: correct for the large-scale illumination of the Chunk
    Exposure (compensates for any brightness gradients in the dome flats).  These
    will be more rare as the time to take them in asronomical twilight may not be
    enough to get these in all filters slated for observing for an evening.
   
    Night Sky Flats: correct for large-scale illumination effects.  These will be
    derived from the Science Chunk Exposures.
    
    NOTE: The bias subtraction sub-stage of the ISR must be run BEFORE this sub-stage.
    
    @return chunkExposure corrected for large-scale illumination

    @throw LengthError if chunk and master Exposures are different sizes
    @throw RangeError if chunk and master are derived from different pixels 
    @throw NotFound if any Policy or metadata value can not be obtained
   
    TO DO (as of 11/24/08):
    - handle stretch and scale factors better, if needed
    - QQQ: do we need to sig-clip here?
    - QQQ: stretchFactor - do we need it??
    """

    stage = "lsst.ip.isr.flatFieldCorrection"   
    pexLog.Trace("Entering ISR Stage: ", 4, "%s" % (stage,))

    # Parse the Policy File
    try:
        flatPolicy = isrPolicy.getPolicy("flatPolicy")
        chunkType = flatPolicy.getString("chunkType")
        flatFieldScale = flatPolicy.getDouble("flatFieldScale")
        stretchFactor = flatPolicy.getDouble("stretchFactor")
        sigClip = flatPolicy.getBool("sigClip")
        sigClipVal = flatPolicy->getDouble("sigClipVal")
        illumPolicy = isrPolicy.getPolicy("illumPolicy")
        run = isrPolicy.getString("run")
     except pexEx.LsstExceptionStack, e:
         print "Cannot parse the ISR Policy File: %s" % e   
         raise pexExcept.NotFound, "Can not obtain Flat Field Correction Policy Parameters."
   
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
        raise pexExcept.LengthError, "In %s: Chunk Exposure and Master Flat Field Chunk Exposure are not the same size." % (stage,)

    # Check that the Master Flat Field Chunk Exposure and Chunk Exposure are
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
        raise pexExcept.NotFound,"In %s, Chunk Exposure and Master Flat Field Chunk Exposure are not from the same FILTER." % (stage,) 

    # Has the Master Flat Field Chunk Exposure been normalized?

    # CFHT data lists all image processing flags as
    # 'IMRED_processingStep' eg. 'IMRED_NF' = elixir normalized the
    # master flat field.  Will need to put these processing flags into
    # the original exposure.

    isrNormalize = masterChunkMetadata.findUnique("ISR_NC")
    if isrNormalize:
        pexLog.Trace("In %s:" % (stage,), 4, "Master Flat Field Chunk Exposure has been normalized by the ISR.")
    else:
        
        # Normalize the Master Flat Field Chunk Exposure by dividing
        # the Master Flat Field Chunk Exposure by the mean value of
        # the entire Master Flat Field Chunk Exposure.

        ipIsr.easyMean(n, mu, sigma, masterChunkMaskedImage)
        masterChunkMaskedImage /= mu
        
        masterChunkMetadata.addProperty(dafBase.Dataproperty("NORM_MU", mu))
        masterChunkMetadata.addProperty(dafBase.DataProperty("NORM_SD", sigma))
        chunkMetadata.addProperty(dafBase.DataProperty("NORM_MU", mu))
        chunkMetadata.addProperty(dafBase.DataProperty("NORM_SD", sigma))
        masterChunkMetadata.addProperty(dafBase.DataProperty("NORM_END", "Completed Successfully")) 
        masterChunkMaskedImage.setMetadata(masterChunkMetadata)
        chunkMaskedImage.setMetadata(chunkMetadata)

    # Has an Illumination Correction been previously applied to the Chunk Exposure?  
  
    isrIllumination = masterChunkMetadata.findUnique("ISR_ILL");
    if isrIllumination:
        pexLog.Trace("In %s:" % (stage,), 4, "Master Flat Field Chunk Exposure has been illumination corrected.")
    else:

        # Correct the Master Dome (or Twilight) Flat Field Chunk
        # Exposure for scattered light.

        # This correction is different depending on which processing
        # pipelie the ISR stage is being run in.  So first, determine
        # if the current ISR pipeline being run for Data Release or
        # Nightly Processing?

        if run == "DR":
            masterSfChunkExposure = afwImage.ExposureD() # Master Night Sky Flat Field Chunk Exposure
            sfCurrent = illumPolicy.getString("sfCurrent")
            masterSfChunkExposure.readFits(sfCurrent)
            ipIsrIllum.illuminationCorrectionDR(masterChunkExposure, masterSfChunkExposure, isrPolicy)
  
        if run == "nightly":
            masterIcpChunkMaskedImage = afwImage.MaskedImageD() # Master Night Sky Flat Field Chunk Exposure from a previous night
            icPrevious = illumPolicy.getString("icPrevious")
            masterIcpChunkMaskedImage.readFits(icPrevious)
            masterDfpChunkExposure = afwImage.ExposureD() # Master Dome (or Twilight) Flat Field Chunk Exposure from a previous night
            dfPrevious = illumPolicy.getString("dfPrevious")
            masterDfpChunkExposure.readFits(dfPrevious)
            ipIsrIllum.illuminationCorrection(masterChunkExposure, masterDfpChunkExposure, masterIcpChunkMaskedImage, isrPolicy)    

    # Divide the Chunk Exposure by the normalized Master Flat Field
    # Chunk Exposure.  Allow for aditional scaling if desired.

    # QQQ: Do we need to preserve dynamic range by stretching 65K ADU by some factor??
    masterChunkMaskedImage *= stretchFactor



    if flatFieldScale:
        masterChunkMaskedImage *= flatFieldScale
        chunkMaskedImage /= masterChunkMaskedImage
    else:
        chunkMaskedImage /= masterChunkMaskedImage
        
    # Record final stage provenance to the Image Metadata
    dateTime = dafBase .DateTime.utc2mjd()
    chunkMetadata.addProperty(dafBase.Dataproperty("FLAT_MJD", dateTime))
    chunkMetadata.addProperty(dafBase.DataProperty("FLAT_MC", fileName))
    chunkMetadata.addProperty(dafBase.DataProperty("FLAT_END", "Completed Successfully")) 
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
                              
