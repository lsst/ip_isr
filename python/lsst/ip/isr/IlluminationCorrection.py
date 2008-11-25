"""
@brief: Implementation of the stage, Illumination Correction, for the
 Data Release Instrument Signature Removal Pipeline

@author: Nicole M. Silvestri
         University of Washington
         nms@astro.washington.edu

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
   
def illuminationCorrectionDR(masterChunkExposure, masterSFChunkExposure, isrPolicy):

    """Illumination Correction (for the Data Release Pipeline)

    @brief Create a Master Illumination Correction Chunk Exposure from the
    Master Dome and Night Sky Flat Field Chunk Exposures to account for the
    differences in dome vs. sky scattered-light illumination.  This corrects for
    structure in the system response for large-scale surveys.  This method is to
    be employed in the Data Release Pipeline (not for nightly processing)
    because it requires the night sky flat produced from a median-combine of the
    science images from the evening.  This method follows the proposed
    illumination correction described by A. Rest et al. (10/19/2007).
    
    @return masterChunkExposure corrected for illumination

    @throw LengthError if chunk and master Exposures are different sizes
    @throw RangeError if chunk and master are derived from different pixels 
    @throw NotFound if any Policy or metadata value can not be obtained

     The assumption is that the Master Night Sky Flat Field Chunk Exposure has
    had the following performed on it before it arrives at this stage: 
   
    (1) Each individual image was calibrated with the Master Dome (or Twilight)
        Flat Field Chunk Exposure.
    (2) All stars have been masked out of each individual image
    (3) All images have been normalized to teh same average sky value.
    These images are then combined to produce the Master Night Sky Flat Field
    Chunk Exposure used in this sub-stage.
   
    The illumination correction (I) is described as follows:
    I = smoothK((sum(Fs)/sum(Fd))^(-1))
    where Fs = Master Night Sky Flat Field Chunk Exposure
          Fd = Master Dome (or Twilight) Flat Field Chunk Exposure
          smoothK = smoothing kernel  
   
    The final illumination corrected Master Dome (or Twilight) Flat
    Field Chunk Exposure (Fi) is described as follows:    
    Fi = Fd * I 
    where Fd and I are normalized to 1.0.
   
    QQQ: filter dependence of ilumination corrections??
   
    TO DO (as of 11/25/08):
    - need to add code for a smoothing kernel
    - needs to be cleaned of bugs after scons & instantiation
    - need to add this use case version to the ISR's EA model
    """

    stage = "lsst.ip.isr.illuminationCorrectionDR"   
    pexLog.Trace("Entering ISR Stage: ", 4, "%s" % (stage,))

    # Parse the Policy File
    try:
        illumPolicy = isrPolicy.getPolicy("illuminationPolicy") 
        chunkType = illumPolicy.getString("chunkType")
        binSize = illumPolicy.getInt("binSize")
        kernel = illumPolicy.getString("kernel")
        kernelSize = illumPolicy.getInt("kernelSize") 
     except pexEx.LsstExceptionStack, e:
         print "Cannot parse the ISR Policy File: %s" % e   
         raise pexExcept.NotFound, "Can not obtain Illumination Correction Policy Parameters."
   
    masterChunkMaskedImage = masterChunkExposure.getMaskedImage()
    masterChunkMetadata = masterChunkMaskedImage.getImage().getMetaData()
    masterSFChunkMaskedImage = masterSFChunkExposure.getMaskedImage()
    masterSFChunkMetadata = masterSFChunkMaskedImage.getImage().getMetaData()

    # Check that the Master Illumination Chunk Exposure and Chunk Exposure are
    # the same size.

    numCols = masterChunkMaskedImage.getCols()
    numRows = masterChunkMaskedImage.getRows() 

    mnumCols = masterSFChunkMaskedImage.getCols()
    mnumRows = masterSFChunkMaskedImage.getRows() 

    if numCols != mnumCols || numRows != mnumRows:
        raise pexExcept.LengthError, "In %s: Chunk Exposure and Master Flat Field Chunk Exposure are not the same size." % (stage,)

    # Check that the Master Ilumination Chunk Exposure and Chunk Exposure are
    # derived from the same pixels (eg. both are from the same amp,
    # CCD, or raft).

    if chunkType == "amp":
        
        ampidField = masterChunkMetadata.findUnique("AMPID")
        mampidField = masterSFChunkMetadata.findUnique("AMPID")
        if ampField:
            ampid = ampidField.getValue()
            mampid = mampidField.getValue()
        else:
            raise pexExcept.NotFound, "In %s: Could not get AMPID from the Metadata." %(stage,)

    if chunkType == "ccd":
        
        ccdidField = masterChunkMetadata.findUnique("CCDID")
        mccdidField = masterSFChunkMetadata.findUnique("CCDID")
        if ccdField:
            ccdid = ccdidField.getValue()
            mccdid = mccdidField.getValue()
        else:
            raise pexExcept.NotFound, "In %s: Could not get CCDID from the Metadata." %(stage,)

    if chunkType == "raft":
        
        raftidField = masterChunkMetadata.findUnique("RAFTID")
        mraftidField = masterSFChunkMetadata.findUnique("RAFTID")
        if raftField:
            raftid = raftidField.getValue()
            mraftid = mraftidField.getValue()
        else:
            raise pexExcept.NotFound, "In %s: Could not get RAFTID from the Metadata." %(stage,)
        
    if ampid != mampid || ccdid != mccdid || raftid != mraftid:
        raise pexExcept.RangeError, "In %s: Chunk Exposure and Master Bias Chunk Exposure are not derived from the same pixels." % (stage,)

    # Get the relevant metadata 

    fileNameField = masterSFChunkMetadata.findUnique("FILENAME")
    if fileNameField:
        fileName = fileNameField.getValue()
    else:
        raise pexExcept.NotFound,"In %s: Could not get FILENAME from the Metadata." %(stage,) 

    filterField = chunkMetadata.findUnique("FILTER")
    mfilterField = masterSFChunkMetadata.findUnique("FILTER")
    if filterField:
        filter = filterField.getValue()
        mfilter = mfilterField.getValue()
    else:
        raise pexExcept.NotFound,"In %s: Could not get FILTER from the Metadata." %(stage,)

    # Make sure the chunk and master have the same filter designation 
    if filter != mfilter:
        raise pexExcept.NotFound,"In %s, Chunk Exposure and Master Flat Field Chunk Exposure are not from the same FILTER." % (stage,) 

    # Apply the Illumination Correction
    # assuming that masterSFChunkExposure has been normalized, Fs/Fd, and smoothed with some kernel
    masterChunkMaskedImage *= masterSFChunkMaskedImage

    # Record final stage provenance to the Image Metadata
    dateTime = dafBase .DateTime.utc2mjd()
    masterChunkMetadata.addProperty(dafBase.Dataproperty("ILDR_MJD", dateTime))
    masterChunkMetadata.addProperty(dafBase.DataProperty("ILDR_END", "Completed Successfully")) 
    masterChunkMaskedImage.setMetadata(masterChunkMetadata);
  
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


                              

def illuminationCorrection(masterChunkExposure, masterDfpChunkExposure, masterIcpChunkExposure, isrPolicy):

    """Illumination Correction (for the Nightly Processing Pipeline)

    @brief Correct Master Flat Field Chunk Exposures for the differences in dome
    vs. sky illumination.  This corrects for large scale structure in the system
    * response for large-scale surveys.  This method is for the nightly IPP only
    * because it derives an estimate of the illumination correction (assuming
    that * the global changes in illumination are negligible (or very small)
    from the * Master Dome Flat Field Chunk Exposure from the current night
    (Df(t2), the * Master Dome Flat Field Chunk Exposure from the provious night
    (or one close * in time; Df(t1)) and the illumination correction from the
    previous night (or * one close in time.  I(t2)) as follows:
   
    I(t2) = smoothK(Df(t1)/Df(t2))*I(t2) 
    where smoothK is the smoothing kernel
   
    @return masterChunkExposure corrected for illumination

    NOTE: masterChunkExposure = Dome (or Twilight) Flat Field Chunk Exposure (current night)
          masterDfpChunkExposure = the Master Dome (or Twilight) Flat Field Chunk Exposure (previous night)
          masterIcpChunkExposure = the Master Illumination Correction Chunk Exposure (previous night)  
   
    @throw NotFound if any Policy or metadata value can not be obtained
    @throw InvalidParameter if functional form for the lineaization fit is invalid
    
    QQQ: filter dependence of ilumination corrections.
    
    TO DO (as of  11/25/2008):
    - need to add this version of the code to the ISR's EA model
    """

    stage = "lsst.ip.isr.illuminationCorrection"   
    pexLog.Trace("Entering ISR Stage: ", 4, "%s" % (stage,))

    # Parse the Policy File
    try:
        illumPolicy = isrPolicy.getPolicy("illuminationPolicy") 
        chunkType = illumPolicy.getString("chunkType")
        binSize = illumPolicy.getInt("binSize")
        kernel = illumPolicy.getString("kernel")
        kernelSize = illumPolicy.getInt("kernelSize")
        illumMiName = illumPolicy.getString("illumMiName")
     except pexEx.LsstExceptionStack, e:
         print "Cannot parse the ISR Policy File: %s" % e   
         raise pexExcept.NotFound, "Can not obtain Illumination Correction Policy Parameters."
   
    masterChunkMaskedImage = masterChunkExposure.getMaskedImage()
    masterChunkMetadata = masterChunkMaskedImage.getImage().getMetaData()
    masterDfpChunkMaskedImage = masterDfpChunkExposure.getMaskedImage()
    masterDfpChunkMetadata = masterDfpChunkMaskedImage.getImage().getMetaData()
    masterIcpChunkMaskedImage = masterIcpChunkExposure.getMaskedImage()
    masterIcpChunkMetadata = masterIcpChunkMaskedImage.getImage().getMetaData()
    
    # Check that the Master Illumination Chunk Exposure and Master Dome Chunk Exposures are
    # the same size.

    numCols = masterChunkMaskedImage.getCols()
    numRows = masterChunkMaskedImage.getRows() 

    mnumCols = masterDfpChunkMaskedImage.getCols()
    mnumRows = masterDfpChunkMaskedImage.getRows()

    inumCols = masterIcpChunkMaskedImage.getCols()
    inumRows = masterIcpChunkMaskedImage.getRows() 
    

    if numCols != mnumCols || numRows != mnumRows || numCols != inumRows || numRows != inumRows:
        raise pexExcept.LengthError, "In %s: Master Illumination Chunk Exposure and Master Flat Field Chunk Exposures are not the same size." % (stage,)

    # Check that the Master Ilumination Chunk Exposure and Master Dome
    # Chunk Exposures are derived from the same pixels (eg. all are
    # from the same amp, CCD, or raft).

    if chunkType == "amp":
        
        ampidField = masterChunkMetadata.findUnique("AMPID")
        mampidField = masterDfpChunkMetadata.findUnique("AMPID")
        iampidField = masterIcpChunkMetadata.findUnique("AMPID")
        if ampField:
            ampid = ampidField.getValue()
            mampid = mampidField.getValue()
            iampid = iampidField.getValue()
        else:
            raise pexExcept.NotFound, "In %s: Could not get AMPID from the Metadata." %(stage,)

    if chunkType == "ccd":
        
        ccdidField = masterChunkMetadata.findUnique("CCDID")
        mccdidField = masterDfpChunkMetadata.findUnique("CCDID")
        iccdidField = masterIcpChunkMetadata.findUnique("CCDID")
        if ccdField:
            ccdid = ccdidField.getValue()
            mccdid = mccdidField.getValue()
            iccdid = iccdidField.getValue()
        else:
            raise pexExcept.NotFound, "In %s: Could not get CCDID from the Metadata." %(stage,)

    if chunkType == "raft":
        
        raftidField = masterChunkMetadata.findUnique("RAFTID")
        mraftidField = masterDfpChunkMetadata.findUnique("RAFTID")
        iraftidField = masterIcpChunkMetadata.findUnique("RAFTID")
        if raftField && mraftidField && iraftidField:
            raftid = raftidField.getValue()
            mraftid = mraftidField.getValue()
            iraftid = iraftidField.getValue()
        else:
            raise pexExcept.NotFound, "In %s: Could not get RAFTID from the Metadata." %(stage,)
        
    if ampid != mampid || ccdid != mccdid || raftid != mraftid || ampid != iampid || ccdid != iccdid || raftid != iraftid:
        raise pexExcept.RangeError, "In %s: Master Illumination Chunk Exposure and Master Dome Flat Field Chunk Exposures are not derived from the same pixels." % (stage,)

    # Get the relevant metadata 

    fileNameField = masterDfpChunkMetadata.findUnique("FILENAME")
    mfileNameField = masterDfpChunkMetadata.findUnique("FILENAME")
    ifileNameField = masterIcpChunkMetadata.findUnique("FILENAME")
    if fileNameField && ifileNameField:
        fileName = fileNameField.getValue()
        mfilename = mfileNameField.getValue()
        ifileName = ifileNameField.getValue()
    else:
        raise pexExcept.NotFound,"In %s: Could not get FILENAME from the Metadata." %(stage,) 

    filterField = chunkMetadata.findUnique("FILTER")
    mfilterField = masterDfpChunkMetadata.findUnique("FILTER")
    ifilterField = masterIcpChunkMetadata.findUnique("FILTER")
    if filterField && mfilterField && ifilterField:
        filter = filterField.getValue()
        mfilter = mfilterField.getValue()
        ifilter = ifilterField.getValue()
    else:
        raise pexExcept.NotFound,"In %s: Could not get FILTER from the Metadata." %(stage,)

    # Make sure the chunk and master have the same filter designation 
    if filter != mfilter != ifilter:
        raise pexExcept.NotFound,"In %s, Master Illumination Chunk Exposure and Master Dome Flat Field Chunk Exposures are not from the same FILTER." % (stage,) 

    # Have all of the Master Chunk Exposures been normalized?  The
    # 'Flat Field Correction' stage checks the current night's
    # masterChunkExposure
    
    misrNormalize = masterDfpChunkMetadata.findUnique("ISR_NC")
    if misrNormalize:
        pexLog.Trace("In %s:" % (stage,), 4, "Master Flat Field Chunk Exposure has been normalized by the ISR.")
    else:
        
        # Normalize the Master Dome (or twilight) Flat Field Exposure
        # from a previous night

        ipIsr.easyMean(n, mu, sigma, masterDfpChunkMaskedImage)
        masterDfpChunkMaskedImage /= mu
        
        masterDfpChunkMetadata.addProperty(dafBase.Dataproperty("NORM_MU", mu))
        masterDfpChunkMetadata.addProperty(dafBase.DataProperty("NORM_SD", sigma))      
        masterDfpChunkMetadata.addProperty(dafBase.DataProperty("NORM_END", "Completed Successfully")) 
        masterDfpChunkMaskedImage.setMetadata(masterDfpChunkMetadata)
      
    iisrNormalize = masterIcpChunkMetadata.findUnique("ISR_NC")
    if iisrNormalize:
        pexLog.Trace("In %s:" % (stage,), 4, "Master Flat Field Chunk Exposure has been normalized by the ISR.")
    else:
        
        # Normalize the Master Illumination Correction Exposure from a
        # previous night

        ipIsr.easyMean(n, mu, sigma, masterIcpChunkMaskedImage)
        masterIcpChunkMaskedImage /= mu
        
        masterIcpChunkMetadata.addProperty(dafBase.Dataproperty("NORM_MU", mu))
        masterIcpChunkMetadata.addProperty(dafBase.DataProperty("NORM_SD", sigma))      
        masterIcpChunkMetadata.addProperty(dafBase.DataProperty("NORM_END", "Completed Successfully")) 
        masterIcpChunkMaskedImage.setMetadata(masterIcpChunkMetadata)

    # Divide the previous night's Dome Flat with the current night's Dome Flat.
  
    masterDpfChunkMaskedImage /= masterChunkMaskedImage
    
    # Smooth the new masterDpfChunkMaskedImage with a kernel
    # NOTE: this is goingto be slow.  Bin the data and this will go faster.
    # NEED TO ADD THIS STEP IN (11/25/08)
    smoothedMasterChunkMaskedImage = afwImage.MaskedImageD()
    
    
    # Construct the final illumination correction
    #RETURN THIS TO THE CLIPBOARD...DON'T WRITE IT OUT
    
    smoothedMasterChunkMaskedImage *= masterIcpChunkMaskedImage
    smoothedMetadata = smoothedMasterChunkMaskedImage.getImage().getMetaData() 

    smoothedMetadata.addProperty(dafBase.DataProperty("ILNP_BS", binSize))
    smoothedMetadata.addProperty(dafBase.DataProperty("ILNP_KS", kernelSize))
    smoothedMetadata.addProperty(dafBase.DataProperty("ILNP_KT", kernelType))
    smoothedMasterChunkMaskedImage.setMetadata(smoothedMetadata)
    
    # Apply the Illumination Correction
  
    masterChunkMaskedImage *= smoothedMasterChunkMaskedImage

    # Record final stage provenance to the Image Metadata
    dateTime = dafBase .DateTime.utc2mjd()
    masterChunkMetadata.addProperty(dafBase.Dataproperty("ILNP_MJD", dateTime))
    masterChunkMetadata.addProperty(dafBase.DataProperty("ILNP_END", "Completed Successfully")) 
    masterChunkMaskedImage.setMetadata(masterChunkMetadata)
  
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
