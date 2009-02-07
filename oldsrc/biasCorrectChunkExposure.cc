// -*- LSST-C++ -*- // fixed format comment for emacs
/**
  * \file
  *
  * \ingroup isr
  *
  * \brief Implementation of the templated stage, Bias Correct Chunk Exposure,
  * of the Instrument Signature Removal stage for the nightly and Data Release
  * LSST Image Processing Pipeline.
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

#include "boost/cstdint.hpp"
#include "boost/format.hpp"
#include "boost/shared_ptr.hpp"

#include <lsst/afw/image/Exposure.h>
#include <lsst/afw/image/Mask.h>
#include <lsst/afw/image/MaskedImage.h>
#include <lsst/afw/image/PixelAccessors.h>
#include <lsst/daf/base/DataProperty.h>
#include <lsst/pex/exceptions/Exception.h>
#include <lsst/pex/logging/Trace.h>
#include <lsst/pex/policy/Policy.h>

#include "lsst/ip/isr/isr.h"

/** \brief The appropriate Master Bias Chunk Exposure is retrieved from the
  * Clipboard and is subtracted from the Chunk Exposure to correct for structure
  * in the bias (offsets from the overscan level).
  * 
  * \return chunkExposure corrected for bias level
  *
  * \throw Runtime if sub-stage has been run previously
  * \throw NotFound if any policy or metadata can not be obtained
  * \throw LengthError if chunk and master Exposures are different sizes
  * \throw RangeError if chunk and master are derived from different pixels 
  * 
  * TODO (as of Wed 11/06/08):
  * - implement sigma-clipping
  * - implement verification that images are derived from same raft if chunk=raft? 
  * - add more metadata
  */

std::string biasStage = "lsst.ip.isr.biasCorrectChunkExposure";

template <typename ImageT, typename MaskT>
void lsst::ip::isr::biasCorrectChunkExposure(
    lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure,
    lsst::afw::image::Exposure<ImageT, MaskT> const &masterChunkExposure,
    lsst::pex::policy::Policy &isrPolicy,
    lsst::pex::policy::Policy &datasetPolicy
    ) {

    lsst::pex::logging::TTrace<3>("Entering ISR stage: %s", biasStage);

    // Get the Chunk MaskedImage and Image Metadata from the Chunk Exposure 
    lsst::afw::image::MaskedImage<ImageT, MaskT> chunkMaskedImage = chunkExposure.getMaskedImage();
    lsst::daf::base::DataProperty::PtrType chunkMetadata = chunkMaskedImage.getImage()->getMetaData();

    // Get the Master Bias Chunk MaskedImage and Image Metadata from the Master
    // Bias Chunk Exposure
    lsst::afw::image::MaskedImage<ImageT, MaskT> masterChunkMaskedImage = masterChunkExposure.getMaskedImage();
    lsst::daf::base::DataProperty::PtrType masterChunkMetadata = masterChunkMaskedImage.getImage()->getMetadata();

    // Check that this ISR sub-stage has not been run previously on this Chunk
    // Exposure.  If it has, terminate the stage.  Looking for ISR's bias stage flag, "BIAS_END".
    lsst::daf::base::DataProperty::PtrType isrBiasField = chunkMetadata->findUnique("BIAS_END");
    if (isrBiasField) {
        throw lsst::pex::exceptions::Runtime(std::string("Bias Subtraction previously performed - terminating stage."));
    }

    // Check that the Master Bias Chunk Exposure and Chunk Exposure are the same
    // size.

    const int numCols = static_cast<int>(chunkMaskedImage.getCols());
    const int numRows = static_cast<int>(chunkMaskedImage.getRows()); 

    const int mnumCols = static_cast<int>(masterChunkMaskedImage.getCols());
    const int mnumRows = static_cast<int>(masterChunkMaskedImage.getRows()); 

    if (numCols != mnumCols || numRows != mnumRows) {
        throw lsst::pex::exceptions::LengthError(std::string("In ") + __func__ + std::string(": Chunk Exposure and Master Bias Chunk Exposure are not the same size."));
    }


    // Check that the Master Bias Chunk Exposure and Chunk Exposure are derived
    // from the same pixels (eg. both are from the same amp, CCD, or raft).
  
    lsst::pex::policy::Policy::Ptr biasPolicy = isrPolicy.getPolicy("biasPolicy"); 
    std::string chunkType = biasPolicy->getString("chunkType");
    if (chunkType == "amp") {
        
        lsst::daf::base::DataProperty::PtrType ampidField = chunkMetadata->findUnique("AMPID");
        int ampid;
        if (ampidField) {
            ampid = boost::any_cast<const int>(ampidField->getValue());
        } else {
            throw lsst::pex::exceptions::NotFound(std::string("In ") + __func__ + std::string(": Could not get AMPID from the Chunk Metadata."));
        }
     
        lsst::daf::base::DataProperty::PtrType mampidField = masterChunkMetadata->findUnique("AMPID");
        int mampid;
        if (mampidField) {
            mampid = boost::any_cast<const int>(mampidField->getValue());
        } else {
            throw lsst::pex::exceptions::NotFound(std::string("In ") + __func__ + std::string(": Could not get AMPID from the Master Bias Chunk Metadata."));
        }
     
        if (ampid != mampid) {
            throw lsst::pex::exceptions::RangeError(std::string("In ") + __func__ + std::string(": Chunk Exposure and Master Bias Chunk Exposure are not derived from the same pixels."));
        }
    // CHECK IT IF ITS A CCD
    } else if (chunkType == "ccd") {
        lsst::daf::base::DataProperty::PtrType ccdidField = chunkMetadata->findUnique("CCDID");
        int ccdid;
        if (ccdidField) {
            ccdid = boost::any_cast<const int>(ccdidField->getValue());
        } else {
            throw lsst::pex::exceptions::NotFound(std::string("In ") + __func__ + std::string(": Could not get CCDID from the Chunk Metadata."));
        }
     
        lsst::daf::base::DataProperty::PtrType mccdidField = masterChunkMetadata->findUnique("CCDID");
        int mccdid;
        if (mccdidField) {
            mccdid = boost::any_cast<const int>(mccdidField->getValue());
        } else {
            throw lsst::pex::exceptions::NotFound(std::string("In ") + __func__ + std::string(": Could not get CCDID from the Master Bias Chunk Metadata."));
        }
     
        if (ccdid != mccdid) {
            throw lsst::pex::exceptions::RangeError(std::string("In ") + __func__ + std::string(": Chunk Exposure and Master Bias Chunk Exposure are not derived from the same pixels."));
        }
    // CHECK IT IF ITS A RAFT 
    } else {
        // raft level check
        // not yet implemented
    }

    // Parse the ISR Policy file for bias sub-stage information
   
    double biasScale = biasPolicy->getDouble("biasScale");
    bool sigClip = biasPolicy->getBool("sigClip");
    if (sigClip = true) {
        double sigClipVal = biasPolicy->getDouble("sigClipVal");
        // add sigClipping here - not yet implemented
    }
    
    // Get the relevant metadata 

    lsst::daf::base::DataProperty::PtrType fileNameField = masterChunkMetadata->findUnique("FILENAME");
    std::string fileName;
    if (fileNameField) {
        fileName = boost::any_cast<std::string>(fileNameField->getValue());
    }

    lsst::daf::base::DataProperty::PtrType meanBiasField = masterChunkMetadata->findUnique("MEAN");
    doubel meanBias;
    if (meanBiasField) {
        meanBias = boost::any_cast<double>(meanBiasField->getValue());
    }
        

    // Subtract the Master Bias Chunk Exposure from the Chunk Exposure.

    if (biasScale) {
        masterChunkMaskedImage *= biasScale;
        chunkMaskedImage -= masterChunkMaskedImage;
    } else {
        chunkMaskedImage -= masterChunkMaskedImage;
    }

    // Record the sub-stage provenance to the Image Metadata


    lsst::daf::base::DataProperty::PtrType isrMasterCalFlag(new lsst::daf::base::DataProperty("BIAS_MC", fileName));
    chunkMetadata->addProperty(isrMasterCalFlag);
    // ISR version number std::string version; 
    // date run std::string date;
   
    lsst::daf::base::DataProperty::PtrType isrMeanFlag(new lsst::daf::base::DataProperty("BIAS_MU", meanBias));
    chunkMetadata->addProperty(isrMeanFlag);
   
    lsst::daf::base::DataProperty::PtrType isrEndFlag(new lsst::daf::base::DataProperty("BIAS_END", std::sring("Completed Successfully")));
    chunkMetadata->addProperty(isrEndFlag);

    chunkMaskedImage.setMetadata(chunkMetadata);

    // Calculate additional the SDQA metrics here ?

    // Issue a logging message indicating that the sub-stage executed without issue

    lsst::pex::logging::TTrace<3>("ISR stage: %s completed successfully.", biasStage);
    lsst::pex::logging::TTrace<3>("Leaving ISR stage: %s", biasStage);
	
}

/************************************************************************/
/* Explicit instantiations */

template
void lsst::ip::isr::biasCorrectChunkExposure(
    lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> &chunkExposure,
    lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> const &masterChunkExposure,
    lsst::pex::policy::Policy &isrPolicy,
    lsst::pex::policy::Policy &datasetPolicy
    );

template
void lsst::ip::isr::biasCorrectChunkExposure(
    lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> &chunkExposure,
    lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> const &masterChunkExposure,
    lsst::pex::policy::Policy &isrPolicy,
    lsst::pex::policy::Policy &datasetPolicy
    );
/************************************************************************/
