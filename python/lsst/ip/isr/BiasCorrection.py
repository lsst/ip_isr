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
import lsst.pex.exceptions as pexExcept
import lsst.pex.logging as pexLog
import lsst.pex.policy as pexPolicy
import lsst.ip.isr as ipIsr
from VerifyMasterFile import VerifyMasterFile
import time

# global variables
STAGE_SIGNATURE = 'ISR_BIAS'

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

    if not VerifyMasterFile(chunkExposure, masterExposure, stage, isrPolicy.getPolicy("biasPolicy")):
        raise pexExcept.LsstException, "%s: Can not verify Master file for Bias correction" % (stage,)
        
    # Parse the Policy File
    pexLog.Trace("%s" % (stage,), 4, "Parsing the ISR Policy File." )
    try:
        chunkType     = isrPolicy.getString("chunkType")
        biasPolicy    = isrPolicy.getPolicy("biasPolicy")
        
        biasScale     = biasPolicy.getDouble("biasScale")
        fileNameField = biasPolicy.getString("fileName")
        sigClip       = biasPolicy.getBool("sigClip")
        
        if sigClip == "true":
            sigClipVal = biasPolicy.getDouble("sigClipVal")
            # add additional sigClipping policy info here - not yet implemented
    except:
        pexLog.Trace("%s" % (stage,), 4, "Can not parse the ISR Policy File.")
        raise pexExcept.LsstException, "%s: Can not obtain policy parameters from the ISR Policy File." % (stage,)

    # Get the MaskedImage and Metadata
    chunkMaskedImage = chunkExposure.getMaskedImage()
    chunkMetadata    = chunkExposure.getMetadata()
    if chunkMetadata == None:
        chunkMetadata = dafBase.PropertySet()
    if chunkMetadata.exists(STAGE_SIGNATURE):
        pexLog.Trace("%s" % (stage,), 4, "BiasCorrection has already been run")
        return

    masterMetadata    = masterExposure.getMetadata()
    masterMaskedImage = masterExposure.getMaskedImage()
    # Get additional metadata 
    pexLog.Trace("%s" % (stage,), 4, "Obtaining additional parameters from the metadata." )
    fileName = masterMetadata.getString(fileNameField, 'None')

    # Subtract the Master Bias MaskedImage from the Chunk MaskedImage.
    pexLog.Trace("%s" % (stage,), 4, "Subtracting Master Bias Exposure from Chunk Exposure." )
    if biasScale:
        # Not sure if we will ever use this feature
        pexLog.Trace("%s" % (stage,), 4, "Additional scale factor applied to Master Bias: %s" % (biasScale,))
        masterMaskedImage *= biasScale
        chunkMaskedImage  -= masterMaskedImage
    else:
        chunkMaskedImage  -= masterMaskedImage
        
    # Record final stage provenance to the Image Metadata
    pexLog.Trace("%s" % (stage,), 4, "Recording final provenance information." )
    chunkMetadata.setString("BIAS_IM", fileName)
    chunkMetadata.setString(STAGE_SIGNATURE, 'Completed Successfully: %s' % (time.asctime()))

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
                              
