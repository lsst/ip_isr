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
import lsst.daf.base as dafBase
import lsst.pex.exceptions as pexEx
import lsst.pex.logging as pexLog
import lsst.pex.policy as pexPolicy
import lsst.ip.isr as ipIsr
import lsst.ip.isr.IlluminationCorrection as ipIsrIllum
   
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
        chunkType = isrPolicy.getString("chunkType")
        run = isrPolicy.getString("run")
        flatPolicy = isrPolicy.getPolicy("flatPolicy")
        illumPolicy = isrPolicy.getPolicy("illumPolicy")
        flatFieldScale = flatPolicy.getDouble("flatFieldScale")
        stretchFactor = flatPolicy.getDouble("stretchFactor")
        sigClip = flatPolicy.getBool("sigClip")
        sigClipVal = flatPolicy.getDouble("sigClipVal")
    except pexEx.LsstExceptionStack, e:
        pexLog.Trace("%s" % (stage,), 4, "Cannot parse the ISR Policy File: %s" % e )
        raise pexEx.NotFound, "Can not obtain Flat Field Correction Policy Parameters."
   
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
        raise pexEx.LengthError, "In %s: Chunk Exposure and Master Flat Field Chunk Exposure are not the same size." % (stage,)
    else:
        pexLog.Trace("%s" % (stage,), 4, "Success: Master and Chunk Exposures are the same size." )

    # Check that the Master Flat Field Chunk Exposure and Chunk Exposure are
    # derived from the same pixels (eg. both are from the same amp,
    # CCD, or raft).
    
    pexLog.Trace("%s" % (stage,), 4, "Verifying Master and Chunk Exposures are derived from the same pixels." )
    if chunkType == "AMP":
        
        ampidField = chunkMetadata.findUnique("AMPLIST")
        mampidField = masterChunkMetadata.findUnique("AMPLIST");
        if ampField:
            ampid = ampidField.getValueString()
            print "AMPID chunk: ", ampid
            mampid = mampidField.getValueString()
            print "AMPID master: ", mampid
        else:
            raise pexEx.NotFound, "In %s: Could not get AMPID from the Metadata." %(stage,)

        if ampid != mampid:
            raise pexEx.RangeError, "In %s: Chunk Exposure and Master Chunk Exposure are not derived from the same pixels." % (stage,)
        else:
            pexLog.Trace("%s" % (stage,), 4, "Success: Master and Chunk Exposures are derived from the same pixels.")

    elif chunkType == "CCD":
        
        ccdidField = chunkMetadata.findUnique("CCDID")
        mccdidField = masterChunkMetadata.findUnique("CCDID");
        if ccdField:
            ccdid = ccdidField.getValue()
            mccdid = mccdidField.getValue()
        else:
            raise pexEx.NotFound, "In %s: Could not get CCDID from the Metadata." %(stage,)

        if ccdid != mccdid:
            raise pexEx.RangeError, "In %s: Chunk Exposure and Master Chunk Exposure are not derived from the same pixels." % (stage,)
        else:
            pexLog.Trace("%s" % (stage,), 4, "Success: Master and Chunk Exposures are derived from the same pixels.")

    elif chunkType == "RAFT":
        
        raftidField = chunkMetadata.findUnique("RAFTID")
        mraftidField = masterChunkMetadata.findUnique("RAFTID");
        if raftField:
            raftid = raftidField.getValueString()
            mraftid = mraftidField.getValueString()
        else:
            raise pexEx.NotFound, "In %s: Could not get RAFTID from the Metadata." %(stage,)

        if raftid != mraftid:
            raise pexEx.RangeError, "In %s: Chunk Exposure and Master Chunk Exposure are not derived from the same pixels." % (stage,)
        else:
            pexLog.Trace("%s" % (stage,), 4, "Success: Master and Chunk Exposures are derived from the same pixels." )

    else:
        raise pexExcept.NotFound, "In %s: Chunk Type Not Implemented. Use 'AMP', 'CCD', or 'RAFT'." % (stage,)

    # Get the relevant metadata  
    pexLog.Trace("%s" % (stage,), 4, "Obtaining additional parameters from the metadata." )
    fileNameField = masterChunkMetadata.findUnique("FILENAME")
    if fileNameField:
        fileName = fileNameField.getValue()
    else:
        raise pexEx.NotFound,"In %s: Could not get FILENAME from the Metadata." %(stage,) 

#    meanFlatField = masterChunkMetadata.findUnique("MEAN");
#    if meanFlatField:
#        meanFlat = meanFlatField.getValue()
#    else:
#        raise pexEx.NotFound,"In %s: Could not get MEAN from the master Metadata." %(stage,)

    flterField = chunkMetadata.findUnique("FILTER")
    mfilterField = masterChunkMetadata.findUnique("FILTER")
    if filterField:
        #QQQ: this will likely be a numerical value for LSST so need to generalize here.
        filter = filterField.getValueString()
        print "FILTER chunk: ", filter
        mfilter = mfilterField.getValueString()
        print "FILTER master: ", mfilter
    else:
        raise pexEx.NotFound,"In %s: Could not get FILTER from the Metadata." % (stage,)

    # Make sure the chunk and master have the same filter designation
    pexLog.Trace("%s" % (stage,), 4, "Verifying Master and Chunk Exposures are the same FILTER." )
    if filter != mfilter:
        raise pexEx.Runtime,"In %s, Chunk Exposure and Master Flat Field Chunk Exposure are not from the same FILTER." % (stage,) 
    else:
        pexLog.Trace("%s" % (stage,), 4, "Success: Master and Chunk Exposures are from the same filter." )

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
        mu = ipIsr.easyMean(masterChunkMaskedImage)
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
            masterSfChunkExposure = afwImage.ExposureD() 
            sfCurrent = illumPolicy.getString("sfCurrent")
            masterSfChunkExposure.readFits(sfCurrent)
            
            ipIsrIllum.illuminationCorrectionDR(masterChunkExposure, masterSfChunkExposure, isrPolicy)
  
        else:
            pexLog.Trace("In %s:" % (stage,), 4, "Applying Illumination Correction for nightly processing.")
            # Master Night Sky Flat Field Chunk Exposure from a previous night
            masterIcpChunkMaskedImage = afwImage.MaskedImageD() 
            icPrevious = illumPolicy.getString("icPrevious")
            masterIcpChunkMaskedImage.readFits(icPrevious)
            
            # Master Dome (or Twilight) Flat Field Chunk Exposure from a previous night
            masterDfpChunkExposure = afwImage.ExposureD() 
            dfPrevious = illumPolicy.getString("dfPrevious")
            masterDfpChunkExposure.readFits(dfPrevious)
            
            ipIsrIllum.illuminationCorrection(masterChunkExposure, masterDfpChunkExposure, masterIcpChunkMaskedImage, isrPolicy)    

    # QQQ: Do we need to preserve dynamic range by stretching 65K ADU by some factor??
    # pexLog.Trace("In %s:" % (stage,), 4, "Applying stretch factor: %s." % (stretchFactor,)
    # masterChunkMaskedImage *= stretchFactor

    # Divide the Chunk Exposure by the normalized Master Flat Field
    # Chunk Exposure.  Allow for additional scaling if desired.

    pexLog.Trace("In %s:" % (stage,), 4, "Dividing MaskedImage by Master Flat Field MaskedImage.")
    
    if flatFieldScale:
        pexLog.Trace("In %s:" % (stage,), 4, "Scaling Flat Field MaksedImage by %s." % (flatFieldScale,))
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
                              
