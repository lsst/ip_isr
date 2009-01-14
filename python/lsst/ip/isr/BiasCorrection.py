"""

@brief Implementation of the stage, Bias Correction, for the Nightly
       or Data Release Instrument Signature Removal Pipelines

@author Nicole M. Silvestri
        University of Washington
        nms@astro.washington.edu
        created: Mon Nov 24, 2008  

@file
"""
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.daf.base as dafBase
import lsst.pex.exceptions as pexEx
import lsst.pex.logging as pexLog
import lsst.pex.policy as pexPolicy
import lsst.ip.isr as ipIsr

def biasCorrection(chunkExposure, masterExposure, isrPolicy):

    """Bias Correction

    @brief The appropriate Master Bias Exposure is subtracted from the
    Chunk Exposure to correct for structure in the bias (offsets from
    the overscan level).

    @return chunkExposure corrected for bias

    @throw LengthError if chunk and master Exposures are different sizes
    @throw RangeError if chunk and master Exposures are derived from different pixels 
    @throw NotFound if any Policy or metadata value can not be obtained
   
    TO DO (as of 11/24/08):
    - add any additional SDQA statistics requested by SDQA team
    """

    stage = "lsst.ip.isr.biasCorrection"   
    pexLog.Trace("%s" % (stage,), 4, "Entering ISR Bias Correction stage." )

    # Parse the Policy File
    pexLog.Trace("%s" % (stage,), 4, "Parsing the ISR Policy File." )
    try:
        biasPolicy = isrPolicy.getPolicy("biasPolicy")
        chunkType = isrPolicy.getString("chunkType")
        biasScale = biasPolicy.getDouble("biasScale")
        chunkField = biasPolicy.getString("chunkField")
        fileNameField = biasPolicy.getString("fileName")
        sigClip = biasPolicy.getBool("sigClip")
        if sigClip == "true":
            sigClipVal = biasPolicy.getDouble("sigClipVal")
            # add additional sigClipping policy info here - not yet implemented
    except pexEx.LsstExceptionStack, e:
        pexLog.Trace("%s" % (stage,), 4, "Can not parse the ISR Policy File." % (e,)) 
        raise pexExcept.NotFound, "in %s: Can not obtain policy parameters from the ISR Policy File." % (stage,)
   
    chunkMaskedImage = chunkExposure.getMaskedImage()
    chunkMetadata = chunkMaskedImage.getImage().getMetaData()
    masterMaskedImage = masterExposure.getMaskedImage()
    masterMetadata = chunkMaskedImage.getImage().getMetaData()

    # Verify that the Master Bias Exposure and Chunk Exposure are the
    # same size.

    pexLog.Trace("%s" % (stage,), 4, "Verifying Master and Chunk Exposures are the same size." )
    numCols = chunkMaskedImage.getCols()
    numRows = chunkMaskedImage.getRows() 
    pexLog.Trace("%s" % (stage,), 4, "Chunk Exposure NumCols, NumRows: %s, %s" % numCols, numRows )
    
    mnumCols = masterMaskedImage.getCols()
    mnumRows = masterMaskedImage.getRows() 
    pexLog.Trace("%s" % (stage,), 4, "Master Exposure NumCols, NumRows: %s, %s" % mnumCols, mnumRows )

    if numCols != mnumCols or numRows != mnumRows:
        raise pexExcept.LengthError, "In %s: Chunk Exposure and Master Bias Exposure are not the same size." % (stage,)
    else:
        pexLog.Trace("%s" % (stage,), 4, "Success: Master and Chunk Exposures are the same size." )

    # Check that the Master Bias Exposure and Chunk Exposure are
    # derived from the same pixels (eg. both are from the same 'chunk').

    pexLog.Trace("%s" % (stage,), 4, "Verifying Master and Chunk Exposures are derived from the same pixels." )
    
    chunkField = chunkMetadata.findUnique(chunkField)
    mchunkField = masterMetadata.findUnique(chunkField);
    if chunkField and mchunkField:
        chunkId = chunkField.getValueString()
        mchunkId = mchunkField.getValueString()
    else:
        raise pexExcept.NotFound, "In %s: Could not get %s info from the Metadata." % (stage, chunkType)

    if chunkId != mchinkId:
        raise pexExcept.RangeError, "In %s: Chunk and Master Bias Exposure are not derived from the same pixels." % (stage,)
    else:
        pexLog.Trace("%s" % (stage,), 4, "Success: Master and Chunk Exposures are derived from the same pixels." )

    # Get additional metadata 
    pexLog.Trace("%s" % (stage,), 4, "Obtaining additional parameters from the metadata." )
    fileNameKey = masterMetadata.findUnique(fileNameField)
    if fileNameKey:
        fileName = fileNameKey.getValueString()
        pexLog.Trace("%s" % (stage,), 4, "Master Bias Filename: %s" % (fileName,))
    else:
        raise pexExcept.NotFound,"In %s: Could not get FILENAME from the Metadata." %(stage,) 

#    meanBiasField = masterMetadata.findUnique("MEAN");
#    if meanBiasField:
#        meanBias = meanBiasField.getValueDouble()
#        print "Mean of Master Bias: ", meanBias
#    else:
#        raise pexExcept.NotFound,"In %s: Could not get MEAN from the master Metadata." %(stage,) 

    # subtract the Master Bias MaskedImage from the Chunk MaskedImage.
    # Allow for aditional scaling if desired.

    pexLog.Trace("%s" % (stage,), 4, "Subtracting Master Bias Exposure from Chunk Exposure." )
    if biasScale:
        pexLog.Trace("%s" % (stage,), 4, "Additional scale factor applied to Master Bias: %s" %(biasScale,))
        masterMaskedImage *= biasScale
        chunkMaskedImage -= masterMaskedImage
    else:
        chunkMaskedImage -= masterMaskedImage
        
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
                              
