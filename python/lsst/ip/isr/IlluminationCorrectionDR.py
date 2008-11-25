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
   
def flatFieldCorrection(masterChunkExposure, masterSFChunkExposure, isrPolicy):

    """Illumination Correction

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

    stage = "lsst.ip.isr.illuminationCorrection"   
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
    masterSFChunkMetadata = masterChunkMaskedImage.getImage().getMetaData()

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
        mampidField = masterSFChunkMetadata.findUnique("AMPID");
        if ampField:
            ampid = ampidField.getValue()
            mampid = mampidField.getValue()
        else:
            raise pexExcept.NotFound, "In %s: Could not get AMPID from the Metadata." %(stage,)

    if chunkType == "ccd":
        
        ccdidField = masterChunkMetadata.findUnique("CCDID")
        mccdidField = masterSFChunkMetadata.findUnique("CCDID");
        if ccdField:
            ccdid = ccdidField.getValue()
            mccdid = mccdidField.getValue()
        else:
            raise pexExcept.NotFound, "In %s: Could not get CCDID from the Metadata." %(stage,)

    if chunkType == "raft":
        
        raftidField = masterChunkMetadata.findUnique("RAFTID")
        mraftidField = masterSFChunkMetadata.findUnique("RAFTID");
        if raftField:
            raftid = raftidField.getValue()
            mraftid = mraftidField.getValue()
        else:
            raise pexExcept.NotFound, "In %s: Could not get RAFTID from the Metadata." %(stage,)
        
    if ampid != mampid || ccdid != mccdid || raftid != mraftid:
        raise pexExcept.RangeError, "In %s: Chunk Exposure and Master Bias Chunk Exposure are not derived from the same pixels." % (stage,)
        }

    # Get the relevant metadata 

    fileNameField = masterSFChunkMetadata.findUnique("FILENAME")
    if fileNameField:
        fileName = fileNameField.getValue()
    else:
        raise pexExcept.NotFound,"In %s: Could not get FILENAME from the Metadata." %(stage,) 

    meanFlatField = masterSFChunkMetadata.findUnique("MEAN");
    if meanFlatField:
        meanFlat = meanFlatField.getValue()
    else:
        raise pexExcept.NotFound,"In %s: Could not get MEAN from the master Metadata." %(stage,)

    flterField = chnkMetadata.findUnique("FILTER")
    mfilterField = masterSFChunkMetadata.findUnique("FILTER")
    if filterField:
        filter = filterField.getValue()
        mfilter = mfilterField.getValue()
    else:
        raise pexExcept.NotFound,"In %s: Could not get FILTER from the Metadata." %(stage,)

    # Make sure the chunk and master have the same filter designation 
    if filter != mfilter:
        raise pexExcept.NotFound,"In %s, Chunk Exposure and Master Flat Field Chunk Exposure are not from the same FILTER." % (stage,) 

    # Has the Master Night Sky Flat Field Chunk Exposure been normalized?

    isrNormalize = masterSFChunkMetadata.findUnique("ISR_NC")
    if isrNormalize:
        pexLog.Trace("In %s:" % (stage,), 4, "Master Flat Field Chunk Exposure has been normalized by the ISR.")
    else:
        
        # Normalize the Master Night Sky Flat Field Chunk Exposure by
        # dividing the Master NS Flat Field Chunk Exposure by the mean
        # value of the entire Master NS Flat Field Chunk Exposure.

        ipIsr.easyMean(n, mu, sigma, masterSFChunkMaskedImage)
        masterSFChunkMaskedImage /= mu
        
        masterSFChunkMetadata.addProperty(dafBase.Dataproperty("NORM_MU", mu))
        masterSFChunkMetadata.addProperty(dafBase.DataProperty("NORM_SD", sigma))
        masterChunkMetadata.addProperty(dafBase.DataProperty("NORM_MU", mu))
        masterChunkMetadata.addProperty(dafBase.DataProperty("NORM_SD", sigma))
        masterSFChunkMetadata.addProperty(dafBase.DataProperty("NORM_END", "Completed Successfully")) 
        masterSFChunkMaskedImage.setMetadata(masterSFChunkMetadata)
        masterChunkMaskedImage.setMetadata(masterChunkMetadata)

    # multiplicative inverse...
    masterTempChunkMaskedImage = afwImage.MaskedImageD()
    masterSfChunkMaskedImage /= masterChunkMaskedImage
    masterTempChunkMaskedImage = 1/masterSfChunkMaskedImage
    
    # Smooth the temporary MaskedImage with a kernel
    # NEED TO WRITE/FINISH THIS!

    masterIcChunkMaskedImage = afwImage.MaskedImageD()
    

    # Construct the final illumination corrected Master Dome (or Twilight) Flat
    # Field Chunk Exposure 

    masterChunkMaskedImage *= masterIcChunkMaskedImage

    #RETURN THIS TO TH CLIPBOARD...DON'T WRITE IT OUT
    # Write the Illumination Correction to Fits Storage
    #std::string illumMiName = illumPolicy->getString("illumMiName");
    #masterIcChunkMaskedImage.writeFits(illumMiName);

    # Record final stage provenance to the Image Metadata
    dateTime = dafBase .DateTime.utc2mjd()
    masterChunkMetadata.addProperty(dafBase.Dataproperty("ILDR_MJD", dateTime))
    masterChunkMetadata.addProperty(dafBase.DataProperty("ILDR_BS", binSize))
    masterChunkMetadata.addProperty(dafBase.DataProperty("ILDR_KS", kernelSize))
    masterChunkMetadata.addProperty(dafBase.DataProperty("ILDR_KT", kernelType))
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
                              
