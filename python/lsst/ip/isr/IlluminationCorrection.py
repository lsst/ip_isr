"""
@brief: Implementation of the stage, Illumination Correction, for the
 Data Release and Nightly Instrument Signature Removal Pipeline

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
    pexLog.Trace("%s" % (stage,), 4, "Entering ISR Illumination Correction for DR Stage")

    # Parse the Policy File
    pexLog.Trace("%s" % (stage,), 4, "Parsing the ISR Policy File." )
    try:
        illumPolicy = isrPolicy.getPolicy("illuminationPolicy") 
        chunkType = isrPolicy.getString("chunkType")
        binSize = illumPolicy.getInt("binSize")
        kernel = illumPolicy.getString("kernel")
        kernelSize = illumPolicy.getInt("kernelSize") 
    except pexEx.LsstExceptionStack, e:
        pexLog.Trace("%s" % (stage,), 4, "Cannot parse the ISR Policy File: %s" % e )
        raise pexExcept.NotFound, "Can not obtain Illumination Correction Policy Parameters."
   
    masterChunkMaskedImage = masterChunkExposure.getMaskedImage()
    masterChunkMetadata = masterChunkMaskedImage.getImage().getMetaData()
    masterSFChunkMaskedImage = masterSFChunkExposure.getMaskedImage()
    masterSFChunkMetadata = masterSFChunkMaskedImage.getImage().getMetaData()

    # Check that the Master Illumination Chunk Exposure and Chunk Exposure are
    # the same size.

    pexLog.Trace("%s" % (stage,), 4, "Verifying Master and Chunk Exposures are the same size." )
    numCols = masterChunkMaskedImage.getCols()
    numRows = masterChunkMaskedImage.getRows() 

    mnumCols = masterSFChunkMaskedImage.getCols()
    mnumRows = masterSFChunkMaskedImage.getRows() 

    if numCols != mnumCols or numRows != mnumRows:
        raise pexExcept.LengthError, "In %s: Chunk Exposure and Master Flat Field Chunk Exposure are not the same size." % (stage,)
    else:
        pexLog.Trace("%s" % (stage,), 4, "Success: Master and Chunk Exposures are the same size." )

    # Check that the Master Ilumination Chunk Exposure and Chunk Exposure are
    # derived from the same pixels (eg. both are from the same amp,
    # CCD, or raft).

    pexLog.Trace("%s" % (stage,), 4, "Verifying Master and Chunk Exposures are derived from the same pixels." )
    if chunkType == "AMP":
        
        ampidField = masterChunkMetadata.findUnique("AMPLIST")
        mampidField = masterSFChunkMetadata.findUnique("AMPLIST")
        if ampField and mampField:
            ampid = ampidField.getValueString()
            print "AMPID chunk: ", ampid
            mampid = mampidField.getValueString()
            print "AMPID master: ", mampid
        else:
            raise pexExcept.NotFound, "In %s: Could not get AMPID from the Metadata." %(stage,)

        if ampid != mampid:
            raise pexEx.RangeError, "In %s: Chunk Exposure and Master Bias Chunk Exposure are not derived from the same pixels." % (stage,)
        else:
            pexLog.Trace("%s" % (stage,), 4, "Success: Master and Chunk Exposures are derived from the same pixels.")

    elif chunkType == "CCD":
        
        ccdidField = masterChunkMetadata.findUnique("CCDID")
        mccdidField = masterSFChunkMetadata.findUnique("CCDID")
        if ccdField and mccdField:
            ccdid = ccdidField.getValueString()
            mccdid = mccdidField.getValueString()
        else:
            raise pexExcept.NotFound, "In %s: Could not get CCDID from the Metadata." %(stage,)

        if ccdid != mccdid:
            raise pexEx.RangeError, "In %s: Chunk Exposure and Master Bias Chunk Exposure are not derived from the same pixels." % (stage,)
        else:
            pexLog.Trace("%s" % (stage,), 4, "Success: Master and Chunk Exposures are derived from the same pixels.")

    elif chunkType == "RAFT":
        
        raftidField = masterChunkMetadata.findUnique("RAFTID")
        mraftidField = masterSFChunkMetadata.findUnique("RAFTID")
        if raftField and mraftField:
            raftid = raftidField.getValueString()
            mraftid = mraftidField.getValueString()
        else:
            raise pexExcept.NotFound, "In %s: Could not get RAFTID from the Metadata." %(stage,)

        if raftid != mraftid:
            raise pexEx.RangeError, "In %s: Chunk Exposure and Master Bias Chunk Exposure are not derived from the same pixels." % (stage,)
        else:
            pexLog.Trace("%s" % (stage,), 4, "Success: Master and Chunk Exposures are derived from the same pixels." )
        
    else:
        raise pexExcept.NotFound, "In %s: Chunk Type Not Implemented. Use 'AMP', 'CCD', or 'RAFT'." % (stage,)
   
    # Get the relevant metadata 

    pexLog.Trace("%s" % (stage,), 4, "Obtaining additional parameters from the metadata." )
    fileNameField = masterSFChunkMetadata.findUnique("FILENAME")
    if fileNameField:
        fileName = fileNameField.getValue()
    else:
        raise pexExcept.NotFound,"In %s: Could not get FILENAME from the Metadata." %(stage,) 

    filterField = chunkMetadata.findUnique("FILTER")
    mfilterField = masterSFChunkMetadata.findUnique("FILTER")
    if filterField and mfilterField:
        filter = filterField.getValueString()
        print "FILTER chunk: ", filter
        mfilter = mfilterField.getValueString()
        print "FILTER master: ", mfilter
    else:
        raise pexExcept.NotFound,"In %s: Could not get FILTER from the Metadata." %(stage,)

    # Make sure the chunk and master have the same filter designation
    pexLog.Trace("%s" % (stage,), 4, "Verifying Master and Chunk Exposures are the same FILTER." )
    if filter != mfilter:
        raise pexExcept.NotFound,"In %s, Chunk Exposure and Master Flat Field Chunk Exposure are not from the same FILTER." % (stage,) 
    else:
        pexLog.Trace("%s" % (stage,), 4, "Success: Master and Chunk Exposures are from the same filter." )
        
    # Apply the Illumination Correction
    # assuming that masterSFChunkExposure has been normalized, Fs/Fd, and smoothed with some kernel
    pexLog.Trace("%s" % (stage,), 4, "Applying the Illumination Correction for DR." )
    masterChunkMaskedImage *= masterSFChunkMaskedImage
    
    # Record final stage provenance to the Image Metadata
    pexLog.Trace("%s" % (stage,), 4, "Recording final provenance information." )
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
    pexLog.Trace("%s" % (stage,), 4, "Entering ISR Illumination Coreection for Nightly Processing Stage")

    # Parse the Policy File
    pexLog.Trace("%s" % (stage,), 4, "Parsing the ISR Policy File." )
    try:
        illumPolicy = isrPolicy.getPolicy("illuminationPolicy") 
        chunkType = isrPolicy.getString("chunkType")
        binSize = illumPolicy.getInt("binSize")
        kernelType = illumPolicy.getString("kernelType")
        kernelWidth = illumPolicy.getInt("kernelWidth")
        kernelCol = illumPolicy.getInt("kernelCol")
        kernelRow = illumPolicy.getInt("kernelRow")
        edgeMaskBit = illumPolicy.getInt("edgeMaskBit")
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

    pexLog.Trace("%s" % (stage,), 4, "Verifying Master Exposures are the same size." )
    numCols = masterChunkMaskedImage.getCols()
    numRows = masterChunkMaskedImage.getRows() 

    mnumCols = masterDfpChunkMaskedImage.getCols()
    mnumRows = masterDfpChunkMaskedImage.getRows()

    inumCols = masterIcpChunkMaskedImage.getCols()
    inumRows = masterIcpChunkMaskedImage.getRows() 
    

    if numCols != mnumCols or numRows != mnumRows or numCols != inumRows or numRows != inumRows:
        raise pexExcept.LengthError, "In %s: Master Illumination Chunk Exposure and Master Flat Field Chunk Exposures are not the same size." % (stage,)
    else:
        pexLog.Trace("%s" % (stage,), 4, "Success: Master and Chunk Exposures are the same size." )

    # Check that the Master Ilumination Chunk Exposure and Master Dome
    # Chunk Exposures are derived from the same pixels (eg. all are
    # from the same amp, CCD, or raft).

    pexLog.Trace("%s" % (stage,), 4, "Verifying Master and Chunk Exposures are derived from the same pixels." )
    if chunkType == "AMP":
        
        ampidField = masterChunkMetadata.findUnique("AMPLIST")
        mampidField = masterDfpChunkMetadata.findUnique("AMPLIST")
        iampidField = masterIcpChunkMetadata.findUnique("AMPLIST")
        if ampField and mampField and iampField:
            ampid = ampidField.getValueString()
            mampid = mampidField.getValueString()
            iampid = iampidField.getValueString()
        else:
            raise pexExcept.NotFound, "In %s: Could not get AMPID from the Metadata." % (stage,)

        if ampid != mampid or ampid != iampid:
            raise pexEx.RangeError, "In %s: Master Chunk Exposures are not derived from the same pixels." % (stage,)
        else:
            pexLog.Trace("%s" % (stage,), 4, "Success: Master Chunk Exposures are derived from the same pixels.")
            
    elif chunkType == "ccd":
        
        ccdidField = masterChunkMetadata.findUnique("CCDID")
        mccdidField = masterDfpChunkMetadata.findUnique("CCDID")
        iccdidField = masterIcpChunkMetadata.findUnique("CCDID")
        if ccdField and mccdField and iccdField:
            ccdid = ccdidField.getValueString()
            mccdid = mccdidField.getValueString()
            iccdid = iccdidField.getValueString()
        else:
            raise pexExcept.NotFound, "In %s: Could not get CCDID from the Metadata." % (stage,)

        if ccdid != mccdid or ccdid != iccdid:
            raise pexEx.RangeError, "In %s: Master Chunk Exposures are not derived from the same pixels." % (stage,)
        else:
            pexLog.Trace("%s" % (stage,), 4, "Success: Master Chunk Exposures are derived from the same pixels.")
        
    elif chunkType == "raft":
        
        raftidField = masterChunkMetadata.findUnique("RAFTID")
        mraftidField = masterDfpChunkMetadata.findUnique("RAFTID")
        iraftidField = masterIcpChunkMetadata.findUnique("RAFTID")
        if raftField and mraftidField and iraftidField:
            raftid = raftidField.getValueString()
            mraftid = mraftidField.getValueString()
            iraftid = iraftidField.getValueString()
        else:
            raise pexExcept.NotFound, "In %s: Could not get RAFTID from the Metadata." % (stage,)

        if raftid != mraftid or rafttid != iraftid:
            raise pexEx.RangeError, "In %s: Master Chunk Exposures are not derived from the same pixels." % (stage,)
        else:
            pexLog.Trace("%s" % (stage,), 4, "Success: Master Chunk Exposures are derived from the same pixels.")

    else:
        raise pexExcept.NotFound, "In %s: Chunk Type Not Implemented. Use 'AMP', 'CCD', or 'RAFT'." % (stage,)
    
    # Get the relevant metadata 
    pexLog.Trace("%s" % (stage,), 4, "Obtaining additional parameters from the metadata." )
    fileNameField = masterDfpChunkMetadata.findUnique("FILENAME")
    mfileNameField = masterDfpChunkMetadata.findUnique("FILENAME")
    ifileNameField = masterIcpChunkMetadata.findUnique("FILENAME")
    if fileNameField and mfileNameField and ifileNameField:
        fileName = fileNameField.getValueString()
        mfilename = mfileNameField.getValueString()
        ifileName = ifileNameField.getValueString()
    else:
        raise pexExcept.NotFound,"In %s: Could not get FILENAME from the Metadata." %(stage,) 

    filterField = chunkMetadata.findUnique("FILTER")
    mfilterField = masterDfpChunkMetadata.findUnique("FILTER")
    ifilterField = masterIcpChunkMetadata.findUnique("FILTER")
    if filterField and mfilterField and ifilterField:
        filter = filterField.getValueString()
        mfilter = mfilterField.getValueString()
        ifilter = ifilterField.getValueString()
    else:
        raise pexExcept.NotFound,"In %s: Could not get FILTER from the Metadata." %(stage,)

    # Make sure the masters have the same filter designation
    pexLog.Trace("%s" % (stage,), 4, "Verifying Master Chunk Exposures are the same FILTER." )
    if filter != mfilter != ifilter:
        raise pexExcept.NotFound,"In %s, Master Flat Field Chunk Exposures are not from the same FILTER." % (stage,) 
    else:
        pexLog.Trace("%s" % (stage,), 4, "Success: Master Chunk Exposures are from the same filter." )
    
    # Have all of the Master Chunk Exposures been normalized?  The
    # 'Flat Field Correction' stage checks the current night's
    # masterChunkExposure
    pexLog.Trace("%s" % (stage,), 4, "Checking Master Flat Field Exposure for normalization." )
    misrNormalize = masterDfpChunkMetadata.findUnique("ISR_NC")
    if misrNormalize:
        pexLog.Trace("In %s:" % (stage,), 4, "Master Flat Field Chunk Exposure has been normalized by the ISR.")
    else:
        
        # Normalize the Master Dome (or twilight) Flat Field Exposure
        # from a previous night
        pexLog.Trace("In %s:" % (stage,), 4, "Normalizing the Master Flat Field Exposure.")
        mu = afwMath.make_Statistics(masterDfpChunkMaskedImage.getImage(), afwMath.MEAN).getValue(afwMath.MEAN)
        masterDfpChunkMaskedImage /= mu
        
        pexLog.Trace("%s" % (stage,), 4, "Recording normalization provenance information." )
        masterDfpChunkMetadata.addProperty(dafBase.Dataproperty("NORM_MU", mu))
        #masterDfpChunkMetadata.addProperty(dafBase.DataProperty("NORM_SD", sigma))      
        masterDfpChunkMetadata.addProperty(dafBase.DataProperty("NORM_END", "Completed Successfully")) 
        masterDfpChunkMaskedImage.setMetadata(masterDfpChunkMetadata)

    pexLog.Trace("%s" % (stage,), 4, "Checking Master Illumination Exposure for normalization." )
    iisrNormalize = masterIcpChunkMetadata.findUnique("ISR_NC")
    if iisrNormalize:
        pexLog.Trace("In %s:" % (stage,), 4, "Master Illumination Chunk Exposure has been normalized by the ISR.")
    else:
        
        # Normalize the Master Illumination Correction Exposure from a
        # previous night
        pexLog.Trace("In %s:" % (stage,), 4, "Normalizing the Master Illumination Exposure.")
        mu = afwMath.make_Statistics(masterIcpChunkMaskedImage.getImage(), afwMath.MEAN).getValue(afwMath.MEAN)
        masterIcpChunkMaskedImage /= mu

        pexLog.Trace("%s" % (stage,), 4, "Recording normalization provenance information." )
        masterIcpChunkMetadata.addProperty(dafBase.Dataproperty("NORM_MU", mu))
        #masterIcpChunkMetadata.addProperty(dafBase.DataProperty("NORM_SD", sigma))      
        masterIcpChunkMetadata.addProperty(dafBase.DataProperty("NORM_END", "Completed Successfully")) 
        masterIcpChunkMaskedImage.setMetadata(masterIcpChunkMetadata)

    # Divide the previous night's Dome Flat with the current night's Dome Flat.

    pexLog.Trace("%s" % (stage,), 4, "Dividing previous night Flat by current night Flat." )
    masterDpfChunkMaskedImage /= masterChunkMaskedImage

    # Will add other kernelTypes later...a simple gaussian is good for
    # now.

    if kernelType == "GAUSSIAN": 
    
        # Smooth the new masterDpfChunkMaskedImage with a kernel NOTE:
        # this is going to be slow.  Bin the maskedImage and this will go
        # faster...

        pexLog.Trace("%s" % (stage,), 4, "Obtaining kernel parameters for %s." % (kernelType,))
        smoothedMasterMaskedImage = afwImage.MaskedImageD()
        gaussianFunc = afwMath.GaussianFunction2D()
        gaussianFunc(kernelWidth, kernelWidth) 
        kernel = afwMath.AnalyticKernel()
        kernel(kernelCol, kernelRow, gaussianFunc) 
        
        # Do the convolution 
        pexLog.Trace("%s" % (stage,), 4, "Convolving the Master Flat with the kernel.")
        smoothedMasterMaskedImage(masterDpfChunkMaskedImage.getDimensions()) 
        afwMath.convolve(smoothedMasterMaskedImage, masterDpfChunkMaskedImage, kernel, edgeMaskBit, true) 

    else:
        raise pexExcept.NotImplemented, "In %s: Invalid kernel type." % (stage,)

    # Construct the final illumination correction
    #RETURN THIS TO THE CLIPBOARD...DON'T WRITE IT OUT

    pexLog.Trace("%s" % (stage,), 4, "Constucting the finel Illumination Correction.")
    smoothedMasterMaskedImage *= masterIcpChunkMaskedImage
    smoothedMetadata = smoothedMasterChunkMaskedImage.getImage().getMetaData() 

    pexLog.Trace("%s" % (stage,), 4, "Recording smoothing provenance information.")
    smoothedMetadata.addProperty(dafBase.DataProperty("ILNP_BS", binSize))
    smoothedMetadata.addProperty(dafBase.DataProperty("ILNP_KS", kernelSize))
    smoothedMetadata.addProperty(dafBase.DataProperty("ILNP_KS", kernelWidth))
    smoothedMetadata.addProperty(dafBase.DataProperty("ILNP_KT", kernelType))
    smoothedMasterChunkMaskedImage.setMetadata(smoothedMetadata)
    
    # Apply the Illumination Correction to the Master Flat Field MaskedImage
    pexLog.Trace("%s" % (stage,), 4, "Applying the Illumination Correction.")
    masterChunkMaskedImage *= smoothedMasterMaskedImage

    # Record final stage provenance to the Image Metadata
    pexLog.Trace("%s" % (stage,), 4, "Recording final provenance information." )
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
