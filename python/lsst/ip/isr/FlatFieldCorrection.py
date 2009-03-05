"""
@brief: Implementation of the stage, Flat Field Correction, for the
 nightly Instrument Signature Removal Pipeline

@author: Nicole M. Silvestri / ACB
 contact: nms@astro.washington.edu
 file created: Tue Nov 25, 2008  

@file
"""
import eups
import os

import lsst.afw.image as afwImage
import lsst.daf.base as dafBase
import lsst.pex.exceptions as pexEx
import lsst.pex.logging as pexLog
import lsst.pex.policy as pexPolicy
import lsst.ip.isr as ipIsr
import lsst.ip.isr.IlluminationCorrection as ipIsrIllum
import time

# global variables
STAGE_SIGNATURE = 'ISR_FLAT'
   
def flatFieldCorrection(chunkExposure, masterChunkExposure, isrPolicy):

    """Flat Field Correction

    @brief Divide the Chunk Exposure by the (normalized, illumination corrcted) Master Flat Field Chunk
    Exposure(s) to correct for pixel-to-pixel variations (eg. optics, vignetting,
    thickness variations, gain, etc.). The Master Flat Field Chunk Exposure can
    be one of several different types of flats (broadband dome, mono-chromatic dome, 
    night sky, twilight) with further sub-divisions into LSST filters (ugrizy) or bandpasses.
   
    Calls 'Illumination Correction' stage (DR or nightly) to perform the
    illumination correction.  
   
    Dome Flats: correct for the pixel-to-pixel variations in the response og the
    CCD.  These will be the 'Stubb's' tunable laser flats or broadband quartz lamp dome flats.
   
    Twilight Flats: correct for the large-scale illumination of the Chunk
    Exposure (compensates for any brightness gradients in the dome flats).  These
    will be more rare as the time to take them in asronomical twilight may not be
    enough to get these in all filters slated for observing for an evening.
   
    Night Sky Flats: correct for large-scale illumination effects.  These will be
    derived from the Science Chunk Exposures.  These have vignetting.
    
    NOTE: The bias subtraction stage of the ISR must be run BEFORE this stage.
    
    @return chunkExposure corrected for large-scale illumination

    @throw LengthError if chunk and master Exposures are different sizes
    @throw RangeError if chunk and master are derived from different pixels 
    @throw NotFound if any Policy or metadata value can not be obtained
    @throw Runtime if chunk and master Exposures have different filter designations
   
    TO DO (as of 12/02/08):
    - handle stretch and scale factors better, if needed
    - QQQ: do we need to sig-clip here?
    - QQQ: stretchFactor - do we need it??
    - generalize header info (eg. FILTER) - in Illum. Corr. stage as well
    """

    stage = "lsst.ip.isr.flatFieldCorrection"   
    pexLog.Trace("%s" % (stage,), 4, "Entering ISR Flat Field Correction Stage.")

    # Parse the Policy File
    pexLog.Trace("%s" % (stage,), 4, "Parsing the ISR Policy File." )
    try:
        chunkType      = isrPolicy.getString("chunkType")
        run            = isrPolicy.getString("run")
        illumPolicy    = isrPolicy.getPolicy("illumPolicy")

        flatPolicy     = isrPolicy.getPolicy("flatPolicy")
        fileNameField  = flatPolicy.getString("fileName")
        flatFieldScale = flatPolicy.getDouble("flatFieldScale")
        stretchFactor  = flatPolicy.getDouble("stretchFactor")
        sigClip        = flatPolicy.getBool("sigClip")
        sigClipVal     = flatPolicy.getDouble("sigClipVal")
    except pexEx.LsstExceptionStack, e:
        pexLog.Trace("%s" % (stage,), 4, "Cannot parse the ISR Policy File: %s" % e )
        raise pexExcept.LsstException, "%s: Can not obtain policy parameters from the ISR Policy File." % (stage,)
   
    chunkMaskedImage = chunkExposure.getMaskedImage()
    chunkMetadata    = chunkExposure.getMetaData()
    if chunkMetadata == None:
        chunkMetadata = dafBase.PropertySet()
    if chunkMetadata.exists(STAGE_SIGNATURE):
        pexLog.Trace("%s" % (stage,), 4, "BiasCorrection has already been run")
        return

    masterChunkMaskedImage = masterChunkExposure.getMaskedImage()
    masterChunkMetadata    = masterChunkExposure.getMetaData()

    # Get additional metadata 
    pexLog.Trace("%s" % (stage,), 4, "Obtaining additional parameters from the metadata." )
    fileName = masterMetadata.getString(fileNameField, 'None')

    # Has the Master Flat Field Chunk Exposure been normalized?
    pexLog.Trace("%s" % (stage,), 4, "Checking Master Flat Field Exposure for normalization." )
    
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

        pexLog.Trace("In %s:" % (stage,), 4, "Normalizing the Master Flat Field Exposure.")
        mu = afwMath.make_Statistics(masterChunkMaskedImage.getImage(), afwMath.MEAN).getValue(afwMath.MEAN)
        masterChunkMaskedImage /= mu

        pexLog.Trace("%s" % (stage,), 4, "Recording normalization provenance information." )
        masterChunkMetadata.addProperty(dafBase.Dataproperty("NORM_MU", mu))
        # masterChunkMetadata.addProperty(dafBase.DataProperty("NORM_SD", sigma))      
        masterChunkMetadata.addProperty(dafBase.DataProperty("NORM_END", "Completed Successfully")) 
        masterChunkMaskedImage.setMetadata(masterChunkMetadata)
      
    # Has an Illumination Correction been previously applied to the Chunk Exposure?  

    pexLog.Trace("%s" % (stage,), 4, "Checking Master Flat for Illumination Correction." )
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
            pexLog.Trace("In %s:" % (stage,), 4, "Applying Illumination Correction for Data Release processing.")

            # Master Night Sky Flat Field Chunk Exposure
            sfCurrent = illumPolicy.getString("sfCurrent")
            masterSfChunkExposure = afwImage.ExposureD(sfCurrent) 
            
            ipIsrIllum.illuminationCorrectionDR(masterChunkExposure, masterSfChunkExposure, isrPolicy)
  
        else:
            pexLog.Trace("In %s:" % (stage,), 4, "Applying Illumination Correction for nightly processing.")
            # Master Night Sky Flat Field Chunk Exposure from a previous night
            icPrevious = illumPolicy.getString("icPrevious")
            masterIcpChunkMaskedImage = afwImage.MaskedImageD(icPrevious) 
            
            # Master Dome (or Twilight) Flat Field Chunk Exposure from a previous night
            dfPrevious = illumPolicy.getString("dfPrevious")
            masterDfpChunkExposure = afwImage.ExposureD(dfPrevious) 
            
            ipIsrIllum.illuminationCorrection(masterChunkExposure, masterDfpChunkExposure, masterIcpChunkMaskedImage, isrPolicy)    

    # QQQ: Do we need to preserve dynamic range by stretching 65K ADU by some factor??
    # pexLog.Trace("In %s:" % (stage,), 4, "Applying stretch factor: %s." % (stretchFactor,)
    # masterChunkMaskedImage *= stretchFactor

    # Divide the Chunk Exposure by the normalized Master Flat Field
    # Chunk Exposure.  Allow for additional scaling if desired.

    pexLog.Trace("In %s:" % (stage,), 4, "Dividing MaskedImage by Master Flat Field MaskedImage.")
    
    if flatFieldScale:
        pexLog.Trace("In %s:" % (stage,), 4, "Scaling Flat Field MaskedImage by %s." % (flatFieldScale,))
        masterChunkMaskedImage *= flatFieldScale
        chunkMaskedImage /= masterChunkMaskedImage
    else:
        chunkMaskedImage /= masterChunkMaskedImage
        
    # Record final stage provenance to the Image Metadata
    pexLog.Trace("%s" % (stage,), 4, "Recording final provenance information." )
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
                              
