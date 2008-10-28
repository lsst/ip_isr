// -*- LSST-C++ -*- // fixed format comment for emacs
/**
  * \file
  *
  * \ingroup isr
  *
  * \brief Implementation of the templated subStage, Saturation Correction for
  * Chunk Exposure, of the Instrument Signature Removal stage of the
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
#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>

#include <boost/cstdint.hpp>
#include <boost/format.hpp>
#include <vw/Math/Functions.h> 
#include <vw/Math/Vector.h> 

#include <lsst/afw/image/Exposure.h>
#include <lsst/afw/image/Mask.h>
#include <lsst/afw/image/MaskedImage.h>
#include <lsst/afw/image/PixelAccessors.h>
#include <lsst/daf/base/DataProperty.h>
#include <lsst/detection/Footprint.h>
#include <lsst/pex/exceptions/Exception.h>
#include <lsst/pex/logging/Trace.h>
#include <lsst/pex/policy/Policy.h>

#include "lsst/ip/isr/isr.h"
#include "lsst/ip/isr/interpolateOverMaskedPixels.h"

/** \brief Detect and mask pixels that are saturated (in the A/D converter or
  * have excessive non-linear response) in the Chunk Exposure.  
  *
  * Grow by additional pixels (as given in the ISR Policy) to mask charge
  * spillover.  Set appropriate bits in the Mask and interpolate over masked
  * pixels using the 'InterpolateOverMaskedPixels utility function.
  *
  * \return chunkExposure with saturated pixels masked and interpolated
  *
  * \throw Runtime if this sub-stage has been run previously on the image
  * \throw NotFound if any metadata parameter can not be obtained
  * 
  * TO DO (as of Wed 10/22/08):
  * - delineate between A/D saturated pixels and other
  * - Calculate SDQA metrics as requested by SDQA team
  * - do a better job at selecting saturated pixel and growing
  *   - don't grow into other saturated pixels or grow'n pixels, etc 
  */

template<typename ImageT, typename MaskT>
lsst::afw::image::Exposure<ImageT, MaskT> saturationCorrectionForChunkExposure(
    lsst::afw::image::Exposure<ImageT, MaskT> const &chunkExposure,    
    lsst::pex::policy::Policy &isrPolicy, 
    lsst::pex::policy::Policy &datasetPolicy
    //std::vector<float> &saturationLookUpTable
    ) { 

    // Start with some setup and information gathering...
    // Get the Chunk MaskedImage and Chunk Metadata from the Chunk Exposure

    lsst::afw::image::MaskedImage<ImageT, MaskT> chunkMaskedImage = chunkExposure.getMaskedImage();    
    typename lsst::afw::image::MaskedImage<ImageT, MaskT>::MaskedImagePtrT chunkMaskedImagePtr; 
    lsst::daf::base::DataProperty::PtrType chunkMetadata = chunkMaskedImage.getImage->getMetadata();

    // Make sure we haven't run this sub-stage previously
     lsst::daf::base::DataProperty::PtrType isrSaturation = chunkMetadata->findUnique("ISR_SATCOR");  // saturation stage processing flag
    int isrSat;
    if (isrSaturation) {
        isrSat = boost::any_cast<const int>(isrSaturation->getValue());
    } else {
    throw lsst::pex::exceptions::Runtime(std::string("In ") + __func__ + std::string(": Could not get sub-stage processing flag from the image metadata"));
    }
    
    if (isrSat != 0) {
        throw lsst::pex::exceptions::Runtime(std::string("Saturation Correction has already been applied to this Chunk Exposure"));
    }

    // Parse the ISR policy file for the saturation sub-stage information
    lsst::pex::policy::Policy::Ptr saturationPolicy = isrPolicy.getPolicy("saturationPolicy");
    bool satTable = saturationPolicy->getBool("satTable");
    bool useDefSat = saturationPolicy->getBool("useDefSat");
    std::string satTableName = saturationPolicy->getString("satTableName"); 
    const int ccdNum = datasetPolicy.getInt("ccdName");
    const int satGrow = saturationPolicy->getInt("grow");
    const int threshold = saturationPolicy->getInt("threshold"); // cushion around sat limit

    // Get the saturation limit for the Chunk Exposure
    double satLimit;   
    lsst::daf::base::DataProperty::PtrType saturationField = chunkMetadata->findUnique("SATURATE");
    if (saturationField) {
        // First, try to get the saturation limit for the chunk form the metadata
        satLimit = boost::any_cast<const double>(saturationField->getValue());

    } else if (satTable = true){
        // next, try to get the saturation limit for the chunk from a lookup table
        for (unsigned int tableIter = satTableName.begin(); tableIter < satTableName.end(); tableIter++){
            if (tableIter == ccdNum){
                satLimit = satTableName[tableIter];
            }   
        }
    } else if (useDefSat = true){
        // use a generic saturation limit from the policy file
        satLimit = datasetPolicy.getDouble("satLimit");       
    } else {
        // can't find the saturation limit for the chunk anywhere, I give up!
        throw lsst::pex::exceptions::NotFound(std::string("Can not get Saturation limit for Chunk Exposure."));
    }  

    // Setup the bad pixel mask plane
    MaskT const satMaskBit = chunkMaskedImage.getMask()->getPlaneBitMask("SAT");
   
    const int numCols = static_cast<int>(chunkMaskedImage.getCols());
    const int numRows = static_cast<int>(chunkMaskedImage.getRows()); 
   
    // Save the saturated pixels in a vector for later use as a footprint. May
    // want a class with more functionality (like CRPixel in detection::CR.h to
    // hold saturated pixel info??

    std::vector<float> satPix; // newly detected saturated pixels
    typedef typename std::vector<float>::iterator satPixIter;  // if need to iterate
    std::vector<float> satPixOut; // final list of grow'd saturated pixels

    // Find all of the saturated pixels
    lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> chunkRowAcc(chunkMaskedImage);

    unsigned int numSatPix = 0;
    for (int chunkRow = 0; chunkRow < numRows; chunkRow++, chunkRowAcc.nextRow()) {
        lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> chunkColAcc = chunkRowAcc;
        for (int chunkCol = 0; chunkCol < numCols; chunkCol++, chunkColAcc.nextCol()) {
            if (*chunkColAcc.image >= satLimit){
                
                // store the saturated pixel
                satPix.push_back(*chunkColAcc.image);
                // Sum the saturated pixels for SDQA
                numSatPix += 1;
            }
        }
    }     
  
    // grow around all of the saturated pixels on the Chunk
    for (std::vector<lsst::detection::Footprint::PtrType>::iterator pixIter = satPix.begin(); pixIter != satPix.end(); ++pixIter) { 

        // Need to create a bounding box so I can turn the saturated pixels into
        // footprints and grow around them
 
        vw::BBox2i const & bbox = (*pixIter)->getBBox(); 
        vw::Vector2i const minVec(bbox.min().x() - satGrow, bbox.min().y() - satGrow); 
        vw::Vector2i const maxVec(bbox.max().x() + satGrow, bbox.max().y() + satGrow); 
        vw::BBox2i const fpBBox(minVec, maxVec); 
              
        // lets turn it into a sub image 
        try { 
            chunkMaskedImagePtr = chunkMaskedImage.getSubImage(fpBBox);                  
        } catch (lsst::pex::exceptions::ExceptionStack &e) { 
            continue; 
        } 
             
        // Create a new footprint with the grow'd box
        lsst::detection::Footprint::PtrType fpGrow(new lsst::detection::Footprint(fpBBox)); 
        satPixOut.push_back(fpGrow);        
    } 

    // Mask all of those pixels
    lsst::detection::setMaskFromFootprintList<MaskT>(chunkMaskedImage.getMask(),satPixOut, satMaskBit);

    // Interpolate over all masked pixels. 
    lsst::ip::isr::interpolateOverMaskedPixels<ImageT, MaskT>(chunkExposure, isrPolicy);

    // Record the sub-stage provenance to the Image Metadata
    chunkMetadata->addProperty(lsst::daf::base::DataProperty('ISR_SATCOR', 'Complete'));
    chunkMaskedImage.setMetadata(chunkMetadata);

    // Calculate additional SDQA Metrics??

    //Issue a logging message if the sub-stage executes without issue to this point!
    lsst::pex::logging::TTrace<7>(std::string("ISR sub-stage") +__func__ + std::string("completed successfully."));
         
}

/************************************************************************/
/* Explicit instantiations */

template
lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> saturationCorrectionForChunkExposure(
    lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> const &chunkExposure,    
    lsst::pex::policy::Policy &isrPolicy, 
    lsst::pex::policy::Policy &datasetPolicy
    //std::vector<float> &saturationLookUpTable
    );

template
lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> saturationCorrectionForChunkExposure(
    lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> const &chunkExposure,    
    lsst::pex::policy::Policy &isrPolicy, 
    lsst::pex::policy::Policy &datasetPolicy
    //std::vector<float> &saturationLookUpTable
    );

/************************************************************************/
