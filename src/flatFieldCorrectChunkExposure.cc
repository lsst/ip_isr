// -*- LSST-C++ -*- // fixed format comment for emacs
/**
  * \file
  *
  * \ingroup isr
  *
  * \brief Implementation of the templated subStage, Flat Field Correct Chunk
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
#include <cctype>
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

/** \brief Divide the Chunk Exposure by the (normalized) Master Flat Field Chunk
  * Exposure(s) to correct for pixel-to-pixel variations (eg. optics, vignetting,
  * thickness variations, gain, etc.). The Master Flat Field Chunk Exposure can
  * be one of potentially three different types of flats (dome, twilight, night
  * sky) with further sub-divisions into LSST filters (ugrizy) or bandpasses.
  *
  * Calls "Illumination Correction" sub-stage (DR or nightly) to perform an
  * illumination correction.
  *
  * Dome Flats: correct for the pixel-to-pixel variations in the response og the
  * CCD.  These will be the 'Stubb's' tunable laser flats.
  *
  * Twilight Flats: correct for the large-scale illumination of the Chunk
  * Exposure (compensates for any brightness gradients in the dome flats).  These
  * will be more rare as the time to take them in asronomical twilight may not be
  * enough to get these in all filters slated for observing for an evening.
  *
  * Night Sky Flats: correct for large-scale illumination effects.  These will be
  * derived from the Science Chunk Exposures.
  * 
  * NOTE: The bias subtraction sub-stage of the ISR must be run BEFORE this sub-stage.
  *
  * \return chunkExposure flat field corrected
  * 
  * \throw Runtime if this sub-stage has been run previously (for the particular flat)
  * \throw LengthError if chunk and master Exposure's are different sizes
  * \throw RangeError if chunk and master Exposure's are derived from different pixels
  * \throw NotFound if any requested Policy or metadata information can not be obtained
  * 
  * TODO (as of Wed 10/22/08):
  * - perform raft-level check for chunk and master Exposures
  * - handle stretch and scale factors better, if needed
  *
  * QQQ: do we need to sig-clip here?
  */

template <typename ImageT, typename MaskT>
lsst::afw::image::Exposure<ImageT, MaskT> flatFieldCorrectChunkExposure(
    lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure,
    lsst::afw::image::Exposure<ImageT, MaskT> &masterChunkExposure, // Master Dome (or Twilight) Flat Field Chunk Exposure
    lsst::pex::policy::Policy &isrPolicy,
    lsst::pex::policy::Policy &datasetPolicy
    ) {

    // Get the Chunk MaskedImage and Image Metadata from the Chunk Exposure 
    lsst::afw::image::MaskedImage<ImageT, MaskT> chunkMaskedImage = chunkExposure.getMaskedImage();   
    lsst::daf::base::DataProperty::PtrType chunkMetadata = chunkMaskedImage.getImage()->getMetadata();

    // Get the Master Flat Field Chunk MaskedImage and Image Metadata from the
    // Master Flat Field Chunk Exposure
    lsst::afw::image::MaskedImage<ImageT, MaskT> masterChunkMaskedImage = masterChunkExposure.getMaskedImage();
    lsst::daf::base::DataProperty::PtrType masterChunkMetadata = masterChunkMaskedImage.getImage()->getMetadata();

    std::string subStage = "Flat Field Correct Chunk Exposure";

    // Check that this ISR sub-stage has not been run previously on this Chunk
    // Exposure.  If it has, terminate the stage.
    lsst::daf::base::DataProperty::PtrType isrFlatField = chunkMetadata->findUnique("ISR_FLATCOR");
    if (isrFlatField) {
        lsst::pex::logging::TTrace<3>("In %s: Exposure has already been corrected.  Terminating ISR sub-stage for this Chunk Exposure.", subStage);
        throw lsst::pex::exceptions::Runtime(std::string("Flat Field correction previously performed."));
    }

    // Check that the Master Flat Field Chunk Exposure and Chunk Exposure are
    // the same size.

    const int numCols = static_cast<int>(chunkMaskedImage.getCols());
    const int numRows = static_cast<int>(chunkMaskedImage.getRows()); 

    const int mnumCols = static_cast<int>(masterChunkMaskedImage.getCols());
    const int mnumRows = static_cast<int>(masterChunkMaskedImage.getRows()); 

    if (numCols != mnumCols || numRows != mnumRows) {
        throw lsst::pex::exceptions::LengthError(std::string("In ") + __func__ + std::string(": Chunk Exposure and Master Flat Field Chunk Exposure are not the same size."));
    }

    // Check that the Master Flat Field Chunk Exposure and Chunk Exposure are
    // derived from the same pixels.

    lsst::pex::policy::Policy::Ptr flatPolicy = isrPolicy.getPolicy("flatPolicy");
    std::string chunkType = flatPolicy->getString("chunkType");
    if (chunkType == "amp") {
        lsst::daf::base::DataProperty::PtrType ampidField = chunkMetadata->findUnique("AMPID");
        unsigned int ampid;
        if (ampidField) {
            ampid = boost::any_cast<const int>(ampidField->getValue());
        } else {
            throw lsst::pex::exceptions::NotFound(std::string("In ") + __func__ + std::string(": Could not get AMPID from the Chunk Metadata."));
        }
   
        lsst::daf::base::DataProperty::PtrType mampidField = masterChunkMetadata->findUnique("AMPID");
        unsigned int mampid;
        if (mampidField) {
            mampid = boost::any_cast<const int>(mampidField->getValue());
        } else {
            throw lsst::pex::exceptions::NotFound(std::string("In ") + __func__ + std::string(": Could not get AMPID from the Master Flat Field Chunk Metadata."));
        }
   
        if (ampid != mampid) {
            throw lsst::pex::exceptions::RangeError(std::string("In ") + __func__ + std::string(": Chunk Exposure and Master Flat Field Chunk Exposure are not derived from the same pixels."));
        }
    // CHECK IT IF ITS A CCD
    } else if (chunkType == "ccd") {
        lsst::daf::base::DataProperty::PtrType ccdidField = chunkMetadata->findUnique("CCDID");
        unsigned int ccdid;
        if (ccdidField) {
            ccdid = boost::any_cast<const int>(ccdidField->getValue());
        } else {
            throw lsst::pex::exceptions::NotFound(std::string("In ") + __func__ + std::string(": Could not get CCDID from the Chunk Metadata."));
        }
   
        lsst::daf::base::DataProperty::PtrType mccdidField = masterChunkMetadata->findUnique("CCDID");
        unsigned int mccdid;
        if (mccdidField) {
            mccdid = boost::any_cast<const int>(mccdidField->getValue());
        } else {
            throw lsst::pex::exceptions::NotFound(std::string("In ") + __func__ + std::string(": Could not get CCDID from the Master Flat Field Chunk Metadata."));
        }
   
        if (ccdid != mccdid) {
            throw lsst::pex::exceptions::RangeError(std::string("In ") + __func__ + std::string(": Chunk Exposure and Master Flat Field Chunk Exposure are not derived from the same pixels."));
        }
    } else {
        // check for raft-level compliance
        // not yet implemented
    }

   // Check that the Master Chunk Exposure and Chunk Chunk Exposure are taken in
   // the same filter

    lsst::daf::base::DataProperty::PtrType filterField = chunkMetadata->findUnique("FILTER");
    int filter;
    if (filterField) {
        //  Determine if the filter field value is a number (1-6?)or a string (ugrizY?) 
        filter = boost::any_cast<const int>(filterField->getValue());
    } else {
        throw lsst::pex::exceptions::NotFound(std::string("In ") + __func__ + std::string(": Could not get FILTER from the Chunk Metadata."));
    }    
    if (isalpha(filter)) 
       lsst::pex::logging::TTrace<3>("In %s: Filter Name: %s", subStage, filter); 
//    } else if {
//        filter equal to LSST numerical designations for filters ...do something else

    lsst::daf::base::DataProperty::PtrType mfilterField = masterChunkMetadata->findUnique("FILTER");
    int mfilter;
    if (mfilterField) {
        //  Determine if the filter field value is a number (1-6?)or a string (ugrizY?) 
        mfilter = boost::any_cast<const int>(mfilterField->getValue());
    } else {
        throw lsst::pex::exceptions::NotFound(std::string("In ") + __func__ + std::string(": Could not get FILTER from the Master Flat Field Chunk Metadata."));
    }    
    if (isalpha(mfilter)) 
        lsst::pex::logging::TTrace<3>("In %s: Filter Name: %s", subStage, mfilter); 
            
    if (mfilterField) {
    
        // assuming that the filter will be a number for lsst (1-6) not a string
        // (ugrizY)??

        mfilter = boost::any_cast<const int>(mfilterField->getValue());
    } else {
        throw lsst::pex::exceptions::NotFound(std::string("In ") + __func__ + std::string(": Could not get AMPID from the Master Flat Field Chunk Metadata."));
    }
   
    if (filter != mfilter) {
        throw lsst::pex::exceptions::DomainError(std::string("In ") + __func__ + std::string(": Chunk Exposure and Master Flat Field Chunk Exposure are not from the same FILTER."));
    }

    // Has the Master Flat Field Chunk Exposure been normalized?

    // CFHT data lists all image processing flags as 'IMRED_processingStep'
    // eg. 'IMRED_NF' = elixir normalized the master flat field.  Will need to
    // ask for processing flags in the policy for all datasets (the
    // datasetSpecificPolicy).

    //std::string normalizeKey = datasetPolicy.getString("normalizeKey");
    //lsst::daf::base::DataProperty::PtrType isrNormalize = chunkMetadata->findUnique(normalizeKey);
    lsst::daf::base::DataProperty::PtrType isrNormalize = chunkMetadata->findUnique("ISR_NORMCOR");
    if (isrNormalize) {
        lsst::pex::logging::TTrace<3>("In %s: Master Flat Field Chunk Exposure has been normalized.", subStage);
    } else {

        // Normalize the Master Flat Field Chunk Exposure by dividing the Master
        // Flat Field Chunk Exposure by the mean value of the entire Master Flat
        // Field Chunk Exposure. Borrowing R. Owen's code in Kaiser-coadd for
        // computing the mean (mu) and standard deviation (sigma).
	
        // Compute the number of elements (n), mean (mu) and standard deviation
        // (sigma)

        const int mnumCols = static_cast<int>(masterChunkMaskedImage.getCols());
        const int mnumRows = static_cast<int>(masterChunkMaskedImage.getRows()); 

        lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> masterChunkRowAcc(masterChunkMaskedImage);
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

	// the normalized Master Flat Field Exposure
        masterChunkMaskedImage /= mu;           

        masterChunkMetadata->addProperty(lsst::daf::base::DataProperty("ISR_NORMCOR"));
        lsst::daf::base::DataProperty::PtrType normCorProp = masterChunkMetadata->findUnique("ISR_NORMCOR");
        std::string exitTrue = "Completed Successfully";
        normCorProp->setValue(boost::any_cast<std::string>(exitTrue));
        masterChunkMetadata->addProperty(lsst::daf::base::DataProperty("ISR_NORMCOR_MEAN"));
        lsst::daf::base::DataProperty::PtrType normMeanProp = masterChunkMetadata->findUnique("ISR_NORMCOR_MEAN");
        normMeanProp->setValue(boost::any_cast<double>(mu));
        masterChunkMetadata->addProperty(lsst::daf::base::DataProperty("ISR_NORMCOR_STDEV"));
        lsst::daf::base::DataProperty::PtrType normSigmaProp = masterChunkMetadata->findUnique("ISR_NORMCOR_STDEV");
        normSigmaProp->setValue(boost::any_cast<double>(sigma));
        masterChunkMaskedImage.setMetadata(masterChunkMetadata);
    }
        

    // Has an Illumination Correction been previously applied to the Chunk Exposure?  
  
    lsst::daf::base::DataProperty::PtrType isrIllumination = chunkMetadata->findUnique("ISR_ILLUMCOR");
    if (isrIllumination) {
        lsst::pex::logging::TTrace<3>("In %s: Master Flat Field Chunk Exposure has been corrected for scattered light.", subStage);
    } else {

        // If not, lets correct the Master Dome (or Twilight) Flat Field Chunk
        // Exposure for scattered light.

        // This correction is different depending on which processing pipeline
        // the ISR stage is being run in.  So first, determine if the current
        // ISR pipeline being run for Data Release or Nightly Processing?

        lsst::pex::policy::Policy::Ptr illumPolicy = isrPolicy.getPolicy("illumPolicy");
        std::string run = isrPolicy.getString("run");
        if (run == "DR"){
            lsst::afw::image::Exposure<ImageT, MaskT> masterSfChunkExposure; // Master Night Sky Flat Field Chunk Exposure
            std::string sfCurrent = illumPolicy->getString("sfCurrent");
            masterSfChunkExposure.readFits(sfCurrent);
            lsst::ip::isr::illuminationCorrectionDR<ImageT, MaskT>(masterChunkExposure, masterSfChunkExposure, isrPolicy, datasetPolicy);
        } 
        if (run == "nightly"){
            lsst::afw::image::MaskedImage<ImageT, MaskT> masterIcpChunkMaskedImage; // Master Night Sky Flat Field Chunk Exposure from a previous night
            std::string icPrevious = illumPolicy->getString("icPrevious");
            masterIcpChunkMaskedImage.readFits(icPrevious);
            lsst::afw::image::Exposure<ImageT, MaskT> masterDfpChunkExposure; // Master Dome (or Twilight) Flat Field Chunk Exposure from a previous night
            std::string dfPrevious = illumPolicy->getString("dfPrevious");
            masterDfpChunkExposure.readFits(dfPrevious);
            lsst::ip::isr::illuminationCorrection<ImageT, MaskT>(masterChunkExposure, masterDfpChunkExposure, masterIcpChunkMaskedImage, isrPolicy, datasetPolicy);
        } 
    }

    // Parse the main ISR Policy file for Flat Field sub-stage parameters.
    double flatFieldScale = flatPolicy->getDouble("flatFieldScale");
    // do we need to preserve dynamic range by stretching 65K ADU by some factor??
    double stretchFactor = flatPolicy->getDouble("stretchFactor");
    masterChunkMaskedImage *= stretchFactor;
    bool sigClip = flatPolicy->getBool("sigClip");
    double sigClipVal = flatPolicy->getDouble("sigClipVal");

    // Divide the Chunk Exposure by the normalized Master Flat Field Chunk
    // Exposure.  Hopefully RHL has fixed the Image class so that it properly
    // computes the varaince...

    if (flatFieldScale) {
        masterChunkMaskedImage *= flatFieldScale;
        chunkMaskedImage /= masterChunkMaskedImage;
    } else {
        chunkMaskedImage /= masterChunkMaskedImage;
    }

    // Record the final sub-stage provenance to the Image Metadata
    chunkMetadata->addProperty(lsst::daf::base::DataProperty("ISR_FLATCOR"));
    lsst::daf::base::DataProperty::PtrType flatCorProp = chunkMetadata->findUnique("ISR_FLATCOR");
    std::string exitTrue = "Completed Successfully";
    flatCorProp->setValue(boost::any_cast<std::string>(exitTrue));
    chunkMaskedImage.setMetadata(chunkMetadata);

    // Calculate additional SDQA metrics. 


    // Issue a logging message indicating that the sub-stage executed without issue
    lsst::pex::logging::TTrace<7>("ISR sub-stage, %s, completed successfully.", subStage);
	
}

/************************************************************************/
/* Explicit instantiations */

template
lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> flatFieldCorrectChunkExposure(
    lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> &chunkExposure,
    lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> &masterChunkExposure,
    lsst::pex::policy::Policy &isrPolicy,
    lsst::pex::policy::Policy &datasetPolicy
    );

template
lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> flatFieldCorrectChunkExposure(
    lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> &chunkExposure,
    lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> &masterChunkExposure,
    lsst::pex::policy::Policy &isrPolicy,
    lsst::pex::policy::Policy &datasetPolicy
    );

/************************************************************************/
