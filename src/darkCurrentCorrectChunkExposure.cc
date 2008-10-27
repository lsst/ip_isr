// -*- LSST-C++ -*- // fixed format comment for emacs
/**
  * \file
  *
  * \ingroup imageproc
  *
  * \brief Implementation of the templated subStage, Dark Current Correct Chunk
  * Exposure, of the Instrument Signature Removal stage for the nightly LSST
  * Image Processing Pipeline.
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
#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>

#include <boost/cstdint.hpp>
#include <boost/format.hpp>
#include <boost/shared_ptr.hpp>

#include <lsst/afw/image/Exposure.h>
#include <lsst/afw/image/Mask.h>
#include <lsst/afw/image/MaskedImage.h>
#include <lsst/afw/image/PixelAccessors.h>
#include <lsst/afw/math/Function.h>
#include <lsst/daf/base/DataProperty.h>
#include <lsst/pex/exceptions/Exception.h>
#include <lsst/pex/logging/Trace.h>
#include <lsst/pex/policy/Policy.h>

#include "lsst/ip/isr/isr.h"

/** \brief The appropriate Master Dark Curent Chunk Exposure is retrieved from
  * the Clipboard, scaled, and subtracted from the Chunk Exposure to correct for
  * the thermal noise contribution of the electronics.
  *
  * \return chunkExposure corrected for the dark current contribution
  *
  * \throw Runtime if this sub-stage has been run previously
  * \throw NotFound if any policy or metadata information can not be obtained
  * \throw LengthError if chunk and master Exposures are different sizes
  * \throw RangeError if chunk and master Exposures are derived from different 
  *        pixels
  * 
  * TODO (as of Wed 10/22/08):
  * - implement raft-level check of chunk and master Exposures  
  */

template <typename ImageT, typename MaskT>
lsst::afw::image::Exposure<ImageT, MaskT> darkCurrentCorrectChunkExposure(
	lsst::afw::image::Exposure<ImageT, MaskT> const &chunkExposure,
	lsst::afw::image::Exposure<ImageT, MaskT> const &masterChunkExposure,
	lsst::pex::policy::Policy &isrPolicy,
	lsst::pex::policy::Policy &datasetPolicy
 	) {

    // Get the Chunk MaskedImage and Image Metadata from the Chunk Exposure 

    lsst::afw::image::MaskedImage<ImageT, MaskT> chunkMaskedImage = chunkExposure.getMaskedImage();
    lsst::daf::base::DataProperty::PtrType chunkMetadata = chunkMaskedImage.getImage->getMetadata();

   // Get the Master Dark Current Chunk MaskedImage and Image Metadata from the
   // Master Dark Current  Chunk Exposure

    lsst::afw::image::MaskedImage<ImageT, MaskT> masterChunkMaskedImage = masterChunkExposure.getMaskedImage();
    lsst::daf::base::DataProperty::PtrType masterChunkMetadata = masterChunkMaskedImage.getImage->getMetadata();

    // Check that this ISR sub-stage has not been run previously on this Chunk
    // Exposure.  If it has, terminate the stage.

    lsst::daf::base::DataProperty::PtrType isrDarkField = chunkMetadata->findUnique("DARK");
    if (isrDarkField) {
        lsst::pex::logging::TTrace<3>(std::string("In ") + __func__ + std::string(": Exposure has already been Dark Current Corrected.  Terminating ISR sub-stage for this Chunk Exposure."));
        throw lsst::pex::exceptions::Runtime(std::string("Dark Current Subtraction previously performed."));
    }

    // Check that the Master Dark Current Chunk Exposure and Chunk Exposure are
    // the same size.

    const int numCols = static_cast<int>(chunkExposure.getCols());
    const int numRows = static_cast<int>(chunkExposure.getRows()); 

    const int mnumCols = static_cast<int>(masterChunkExposure.getCols());
    const int mnumRows = static_cast<int>(masterChunkExposure.getRows()); 

     if (numCols != mnumCols || numRows != mnumRows) {
        throw lsst::pex::exceptions::LengthError(std::string("In ") + __func__ + std::string(": Chunk Exposure and Master Bias Chunk Exposure are not the same size."));
    }

    // Check that the Master Dark Current Chunk Exposure and Chunk Exposure are
    // derived from the same pixels.

     lsst::pex::policy::Policy darkPolicy = isrPolicy.getPolicy("darkPolicy");
     std::string chunkType = darkPolicy.getString("chunkType");

     if (chunkType = "amp") {
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
             throw lsst::pex::exceptions::NotFound(std::string("In ") + __func__ + std::string(": Could not get AMPID from the Master Dark Current Chunk Metadata."));
         }
   
         if (ampid != mampid) {
             throw lsst::pex::exceptions::RangeError(std::string("In ") + __func__ + std::string(": Chunk Exposure and Master Dark Current Chunk Exposure are not derived from the same pixels."));
         }
     } else if (chunkType = "ccd") {
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
             throw lsst::pex::exceptions::NotFound(std::string("In ") + __func__ + std::string(": Could not get CCDID from the Master Dark Current Chunk Metadata."));
         }
   
         if (ccdid != mccdid) {
             throw lsst::pex::exceptions::RangeError(std::string("In ") + __func__ + std::string(": Chunk Exposure and Master Dark Current Chunk Exposure are not derived from the same pixels."));
         }
     } else {
         // raft level check
         // not yet implemented
     }

    // Get the rest of the necessary information from the Image Metadata

    lsst::daf::base::DataProperty::PtrType exptimeField = chunkMetadata->findUnique("EXPTIME");
    float exptime;
    if (exptimeField) {
        exptime = boost::any_cast<const float>(exptimeField->getValue());
    } else {
        throw lsst::pex::exceptions::NotFound(std::string("In ") + __func__ + std::string(": Could not get EXPTIME from Chunk Metadata."));
    }

    lsst::daf::base::DataProperty::PtrType mexptimeField = masterChunkMetadata->findUnique("EXPTIME");
    float mexptime;
    if (mexptimeField) {
        mexptime = boost::any_cast<const float>(mexptimeField->getValue());
    } else {
        throw lsst::pex::exceptions::NotFound(std::string("In ") + __func__ + std::string(": Could not get EXPTIME from Master Dark Current Chunk Metadata."));
    }

    // Parse the ISR policy file for dark current correction information
   
    double darkScale = darkPolicy.getDouble("darkScale");

    // Scale the master dark current exposure by the Chunk Exposure's exposure
    // time if the exposure time for the Master Chunk Exposure is different
    if (exptime != mexptime) {
        double scale = exptime/mexptime; 
        masterChunkExposure *= scale;
    }

    // Subtract the Master Dark Chunk Chunk Exposure from the Chunk Exposure.
    // Hopefully RHL has fixed the Image class so that it properly computes the
    // varaince...

    // additional scaling?
    if (darkScale) {
        chunkExposure -= (masterChunkExposure * darkScale);
    } else {
        chunkExposure -= masterChunkExposure;
    }

     // Record the final sub-stage provenance to the Image Metadata
     chunkMetadata->addProperty(lsst::daf::base::DataProperty("ISR_DARKCOR", "Complete"));
     chunkMaskedImage.setMetadata(chunkMetadata);

     // Calculate additional SDQA metrics here. 

     // Issue a logging message indicating that the sub-stage executed without issue
     lsst::pex::logging::TTrace<7>(std::string("ISR sub-stage") + __func__ + std::string("completed successfully."));
	
}

/************************************************************************/
/* Explicit instantiations */

// template
// lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> darkCurrentCorrectChunkExposure(
//     lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> const &chunkExposure,
//     lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> const &masterChunkExposure,
//     lsst::pex::policy::Policy &isrPolicy,
//     lsst::pex::policy::Policy &datasetPolicy
//     );

// template
// lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> darkCurrentCorrectChunkExposure(
//     lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> const &chunkExposure,
//     lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> const &masterChunkExposure,
//     lsst::pex::policy::Policy &isrPolicy,
//     lsst::pex::policy::Policy &datasetPolicy
//     );
/************************************************************************/
