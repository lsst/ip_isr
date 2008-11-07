// -*- LSST-C++ -*- // fixed format comment for emacs
/**
  * \file
  *
  * \ingroup isr
  *
  * \brief Implementation of the templated subStage, Illumination Correction, of
  * the Instrument Signature Removal stage for the Data Release LSST Image
  * Processing Pipeline.
  *
  * \author Nicole M. Silvestri, University of Washington
  *
  * Contact: nms@astro.washington.edu
  *
  * \version
  *
  * LSST Legalese here...
  */
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

#include "boost/shared_ptr.hpp"
#include "boost/cstdint.hpp"
#include "boost/format.hpp"

#include <lsst/afw/image/Exposure.h>
#include <lsst/afw/math/Function.h>
#include <lsst/afw/image/Mask.h>
#include <lsst/afw/image/MaskedImage.h>
#include <lsst/daf/base/DataProperty.h>
#include <lsst/pex/exceptions/Exception.h>
#include <lsst/pex/logging/Trace.h>
#include <lsst/pex/policy/Policy.h>

#include "lsst/ip/isr/isr.h"

/** \brief Create a Master Illumination Correction Chunk Exposure from the
  * Master Dome and Night Sky Flat Field Chunk Exposures to account for the
  * differences in dome vs. sky scattered-light illumination.  This corrects for
  * structure in the system response for large-scale surveys.  This method is to
  * be employed in the Data Release Pipeline (not for nightly processing)
  * because it requires the night sky flat produced from a median-combine of the
  * science images from the evening.  This method follows the proposed
  * illumination correction described by A. Rest et al. (10/19/2007).
  *
  * \return masterChunkExposure corrected for scattered light
  *
  * \throw Runtime if this sub-stage has been run previously
  * \throw NotFound if any Policy or metadata value can not be obtained
  * \throw LengthError if chunk and master Exposure's are different sizes
  * \throw RangeError if chunk and master Exposure's are derived from different pixels
  * 
  * The assumption is that the Master Night Sky Flat Field Chunk Exposure has
  * had the following performed on it before it arrives at this stage: 
  *
  * (1) Each individual image was calibrated with the Master Dome (or Twilight)
  *     Flat Field Chunk Exposure.
  * (2) All stars have been masked out of each individual image
  * (3) All images have been normalized to teh same average sky value.
  * These images are then combined to produce the Master Night Sky Flat Field
  * Chunk Exposure used in this sub-stage.
  *
  * The illumination correction (I) is described as follows:
  * I = smoothK((sum(Fs)/sum(Fd))^(-1))
  * where Fs = Master Night Sky Flat Field Chunk Exposure
  *       Fd = Master Dome (or Twilight) Flat Field Chunk Exposure
  *       smoothK = smoothing kernel  
  *
  * The final illumination corrected Master Dome (or Twilight) Flat Field Chunk Exposure (Fi)   
  * is described as follows:
  * Fi = Fd * I 
  * where Fd and I are normalized to 1.0.
  *
  * QQQ: filter dependence of ilumination corrections??
  * 
  * TO DO (as of  11/05/2008): 
  * - need to add code for a smoothing kernel
  * - needs to be cleaned of bugs after scons & instantiation
  * - need to add this use case version to the ISR's EA model
  */

typedef double vectorType;
typedef double funcType;

template<typename ImageT, typename MaskT>
void lsst::ip::isr::illuminationCorrectionDR(
    lsst::afw::image::Exposure<ImageT, MaskT> &masterChunkExposure, // the Master Dome (or Twilight) Flat Field Chunk Exposure
    lsst::afw::image::Exposure<ImageT, MaskT> &masterSfChunkExposure, // the Master Night Sky FF Chunk Exposure
    lsst::pex::policy::Policy &isrPolicy,  // the main ISR Policy File containing the Illumination Policy info.
    lsst::pex::policy::Policy &datasetPolicy // policy file with info. specific to the dataset being processed
    ) {

    std::string subStage = "Illumination Correction Data Release";

    // Get the Master Chunk MaskedImage and Image Metadata from the Master Dome FF Chunk Exposure
    lsst::afw::image::MaskedImage<ImageT, MaskT> masterChunkMaskedImage = masterChunkExposure.getMaskedImage();
    lsst::daf::base::DataProperty::PtrType masterChunkMetadata = masterChunkMaskedImage.getImage()->getMetadata();

     // Get the Master Chunk MaskedImage and Image Metadata from the Master Night Sky FF Chunk Exposure
    lsst::afw::image::MaskedImage<ImageT, MaskT> masterSfChunkMaskedImage = masterSfChunkExposure.getMaskedImage();
    lsst::daf::base::DataProperty::PtrType masterSfChunkMetadata = masterSfChunkMaskedImage.getImage()->getMetadata();

    // Check that the Master Dome FF Chunk Exposure and the Master Night Sky FF
    // Chunk Exposure are the same size.

    const int numCols = static_cast<int>(masterChunkMaskedImage.getCols());
    const int numRows = static_cast<int>(masterChunkMaskedImage.getRows()); 

    const int mnumCols = static_cast<int>(masterSfChunkMaskedImage.getCols());
    const int mnumRows = static_cast<int>(masterSfChunkMaskedImage.getRows()); 

    if (numCols != mnumCols || mnumRows != numRows) {
        throw lsst::pex::exceptions::LengthError(std::string("In ") + __func__ + std::string(": Master Dome and Night Sky Flat Field Chunk Exposures are not the same size."));
    }
    
    // Check that the Master Chunk Exposures are derived from the same pixels
    // (eg. both are from the same amp, CCD, or raft).
  
    lsst::pex::policy::Policy::Ptr illumPolicy = isrPolicy.getPolicy("illuminationPolicy"); 
    std::string chunkType = illumPolicy->getString("chunkType");
    if (chunkType == "amp") {
        
        lsst::daf::base::DataProperty::PtrType ampidField = masterChunkMetadata->findUnique("AMPID");
        int ampid;
        if (ampidField) {
            ampid = boost::any_cast<const int>(ampidField->getValue());
        } else {
            throw lsst::pex::exceptions::NotFound(std::string("In ") + __func__ + std::string(": Could not get AMPID from the Master Dome Flat Field Chunk Metadata."));
        }
     
        lsst::daf::base::DataProperty::PtrType mampidField = masterSfChunkMetadata->findUnique("AMPID");
        int mampid;
        if (mampidField) {
            mampid = boost::any_cast<const int>(mampidField->getValue());
        } else {
            throw lsst::pex::exceptions::NotFound(std::string("In ") + __func__ + std::string(": Could not get AMPID from the Master Night Sky Flat Field Chunk Metadata."));
        }
     
        if (ampid != mampid) {
            throw lsst::pex::exceptions::RangeError(std::string("In ") + __func__ + std::string("Master Dome and Night Sky Chunk Exposure are not derived from the same pixels."));
        }
    // CHECK IT IF ITS A CCD
    } else if (chunkType == "ccd") {
        lsst::daf::base::DataProperty::PtrType ccdidField = masterChunkMetadata->findUnique("CCDID");
        int ccdid;
        if (ccdidField) {
            ccdid = boost::any_cast<const int>(ccdidField->getValue());
        } else {
            throw lsst::pex::exceptions::NotFound(std::string("In ") + __func__ + std::string(": Could not get CCDID from the Master Dome Flat Filed Chunk Metadata."));
        }
     
        lsst::daf::base::DataProperty::PtrType mccdidField = masterSfChunkMetadata->findUnique("CCDID");
        int mccdid;
        if (mccdidField) {
            mccdid = boost::any_cast<const int>(mccdidField->getValue());
        } else {
            throw lsst::pex::exceptions::NotFound(std::string("In ") + __func__ + std::string(": Could not get CCDID from the Master Night Sky Flat Field Chunk Metadata."));
        }
     
        if (ccdid != mccdid) {
            throw lsst::pex::exceptions::RangeError(std::string("In ") + __func__ + std::string(": Master Dome and Night Sky Flat Field Chunk Exposure are not derived from the same pixels."));
        }
    // CHECK IT IF ITS A RAFT 
    } else {
        // raft level check
        // not yet implemented
    }

    // Parse the ISR Policy file for bias sub-stage information
    const int binSize = illumPolicy->getInt("binSize");
    std::string kernel = illumPolicy->getString("kernel");
    const int kernelSize = illumPolicy->getInt("kernelSize");
   
    // Normalize the Master Night Sky Flat Field Chunk Exposure - divide by the mean
    
    lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> masterChunkRowAcc(masterSfChunkMaskedImage);
    long int n = 0;
    double sum = 0;
       
    for (int chunkRow = 0; chunkRow < mnumRows; chunkRow++, masterChunkRowAcc.nextRow()) {
        lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> masterChunkColAcc = masterChunkRowAcc;
        for (int chunkCol = 0; chunkCol < mnumCols; chunkCol++, masterChunkColAcc.nextCol()) {       
            n++;
            sum += static_cast<double>(*masterChunkColAcc.image);
        } // for column loop
    } // for row loop     
    
    // the mean
    double mu = sum/static_cast<double>(n);	

    double sumSq = 0;
    for (int chunkRow = 0; chunkRow < mnumRows; chunkRow++, masterChunkRowAcc.nextRow()) {
        lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> masterChunkColAcc = masterChunkRowAcc;
        for (int chunkCol = 0; chunkCol < mnumCols; chunkCol++, masterChunkColAcc.nextCol()) {       
            double val = static_cast<double>(*masterChunkColAcc.image) - mu;
            sumSq += val * val;
        } // for column loop
    } // for row loop     

    //the standard deviation
    double sigma = std::sqrt(sumSq/static_cast<double>(n));

    // the normalized Master Night Sky Flat Field Chunk Exposure
    masterSfChunkMaskedImage /= mu;    

    masterSfChunkMetadata->addProperty(lsst::daf::base::DataProperty("ISR_NORMCOR"));
    lsst::daf::base::DataProperty::PtrType normCorProp = masterSfChunkMetadata->findUnique("ISR_NORMCOR");
    std::string exitTrue = "Completed Successfully";
    normCorProp->setValue(boost::any_cast<std::string>(exitTrue));
    masterSfChunkMetadata->addProperty(lsst::daf::base::DataProperty("ISR_NORMCOR_MEAN"));
    lsst::daf::base::DataProperty::PtrType normMeanProp = masterSfChunkMetadata->findUnique("ISR_NORMCOR_MEAN");
    normMeanProp->setValue(boost::any_cast<double>(mu));
    masterSfChunkMetadata->addProperty(lsst::daf::base::DataProperty("ISR_NORMCOR_STDEV"));
    lsst::daf::base::DataProperty::PtrType normSigmaProp = masterSfChunkMetadata->findUnique("ISR_NORMCOR_STDEV");
    normSigmaProp->setValue(boost::any_cast<double>(sigma));
    masterSfChunkMaskedImage.setMetadata(masterSfChunkMetadata);

    // multiplicative inverse...
    lsst::afw::image::MaskedImage<ImageT, MaskT> masterTempChunkMaskedImage;
    masterSfChunkMaskedImage /= masterChunkMaskedImage;
    //masterTempChunkMaskedImage = 1/masterSfChunkMaskedImage;
    
    // Smooth the temporary MaskedImage with a kernel
    // NEED TO WRITE/FINISH THIS!

    lsst::afw::image::MaskedImage<ImageT, MaskT> masterIcChunkMaskedImage;
    

    // Construct the final illumination corrected Master Dome (or Twilight) Flat
    // Field Chunk Exposure 

    masterChunkMaskedImage *= masterIcChunkMaskedImage;

    //RETURN THIS TOO...DON'T WRITE IT OUT
    // Write the Illumination Correction to Fits Storage
   //  std::string illumMiName = illumPolicy->getString("illumMiName");
//     masterIcChunkMaskedImage.writeFits(illumMiName);

    // Record the final sub-stage provenance to the Image Metadata
    masterChunkMetadata->addProperty(lsst::daf::base::DataProperty("ISR_ILLUMCOR_BIN_SIZE"));
    lsst::daf::base::DataProperty::PtrType binSizeProp = masterChunkMetadata->findUnique("ISR_ILLUMCOR_BIN_SIZE");
    binSizeProp->setValue(boost::any_cast<const int>(binSize));

    masterChunkMetadata->addProperty(lsst::daf::base::DataProperty("ISR_ILLUMCOR_KERNEL_SIZE"));
    lsst::daf::base::DataProperty::PtrType kernelSizeProp = masterChunkMetadata->findUnique("ISR_ILLUMCOR_KERNEL_SIZE");
    kernelSizeProp->setValue(boost::any_cast<const int>(kernelSize));
    
    masterChunkMetadata->addProperty(lsst::daf::base::DataProperty("ISR_ILLUMCOR_KERNEL_TYPE"));
    lsst::daf::base::DataProperty::PtrType kernelTypeProp = masterChunkMetadata->findUnique("ISR_ILLUMCOR_KERNEL_TYPE");
    kernelTypeProp->setValue(boost::any_cast<std::string>(kernel));

    masterChunkMetadata->addProperty(lsst::daf::base::DataProperty("ISR_ILLUMCOR"));
    lsst::daf::base::DataProperty::PtrType illumCorProp = masterChunkMetadata->findUnique("ISR_ILLUMCOR");
    illumCorProp->setValue(boost::any_cast<std::string>(exitTrue));
    masterChunkMaskedImage.setMetadata(masterChunkMetadata);

    // Calculate additional SDQA metrics. 


    // Issue a logging message indicating that the sub-stage executed without issue
    lsst::pex::logging::TTrace<7>("ISR sub-stage, %s, completed successfully.", subStage);

}

/************************************************************************/
/* Explicit instantiations */

template
void lsst::ip::isr::illuminationCorrectionDR(
    lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> &masterChunkExposure,
    lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> &masterSfChunkExposure,
    lsst::pex::policy::Policy &isrPolicy,
    lsst::pex::policy::Policy &datasetPolicy
    );

template
void lsst::ip::isr::illuminationCorrectionDR(
    lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> &masterChunkExposure,
    lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> &masterSfChunkExposure,
    lsst::pex::policy::Policy &isrPolicy,
    lsst::pex::policy::Policy &datasetPolicy
    );

/************************************************************************/
