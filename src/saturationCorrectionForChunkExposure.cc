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
#include <sstream>
#include <vector>
#include <cmath>

#include "boost/cstdint.hpp"
#include "boost/format.hpp"
#include "vw/Math/Functions.h" 
#include "vw/Math/Vector.h" 

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

/** \brief Detect and mask pixels that are saturated in the A/D converter or
  * have excessive non-linear response in the Chunk Exposure.  Uses the
  * lsst::detection::DetectionSet to create footprints of saturated pixels above
  * a threshold as given in the ISR Policy.
  *
  * Grow by additional pixels (as given in the ISR Policy) to mask charge
  * spillover.  Set appropriate bits in the Mask and interpolate over masked
  * pixels using the 'InterpolateOverMaskedPixels' utility function.
  *
  * \return chunkExposure with saturated pixels masked and interpolated
  *
  * \throw Runtime if this sub-stage has been run previously on the image
  * \throw NotFound if any metadata parameter can not be obtained
  * 
  * TO DO (as of Wed 10/29/08):
  * - delineate between A/D saturated pixels and other?
  * - Calculate SDQA metrics as requested by SDQA team
  */

template<typename ImageT, typename MaskT>
lsst::afw::image::Exposure<ImageT, MaskT> saturationCorrectionForChunkExposure(
    lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure,    
    lsst::pex::policy::Policy &isrPolicy, 
    lsst::pex::policy::Policy &datasetPolicy
    //std::vector<float> &saturationLookUpTable
    ) { 

    // Start with some setup and information gathering...
    // Get the Chunk MaskedImage and Chunk Metadata from the Chunk Exposure

    lsst::afw::image::MaskedImage<ImageT, MaskT> chunkMaskedImage = chunkExposure.getMaskedImage();    
    //typename lsst::afw::image::MaskedImage<ImageT, MaskT>::MaskedImagePtrT chunkMaskedImagePtr; 
    lsst::daf::base::DataProperty::PtrType chunkMetadata = chunkMaskedImage.getImage()->getMetadata();
    std::string subStage = "Saturation Correction for Chunk Exposure";

    // Make sure we haven't run this sub-stage previously

     lsst::daf::base::DataProperty::PtrType isrSaturation = chunkMetadata->findUnique("ISR_SATCOR");  // saturation stage processing flag
    if (isrSaturation) {
        throw lsst::pex::exceptions::Runtime(std::string("Saturation Correction has already been applied to this Chunk Exposure."));
    } 
    
    // Parse the ISR policy file for the saturation sub-stage information
    // First get the Saturation Sub-Stage Policy from the ISR Policy

    lsst::pex::policy::Policy::Ptr saturationPolicy = isrPolicy.getPolicy("saturationPolicy");
    bool satTable = saturationPolicy->getBool("satTable");
    bool useDefSat = saturationPolicy->getBool("useDefSat");
    std::vector<double> satTableName = saturationPolicy->getDoubleArray("satTableName");
    // std::string satTableName = saturationPolicy->getString("satTableName"); 
    const int ccdNum = datasetPolicy.getInt("ccdName");
    const int satGrow = saturationPolicy->getInt("grow");
    double threshold = saturationPolicy->getDouble("threshold");
    // const int nSatPixMin = saturationPolicy->getInt("nSatPixMin");

    // Get the saturation limit for the Chunk Exposure

    double satLimit;   
    lsst::daf::base::DataProperty::PtrType saturationField = chunkMetadata->findUnique("SATURATE");
    if (saturationField) {
        // First, try to get the saturation limit for the chunk form the metadata
        satLimit = boost::any_cast<const double>(saturationField->getValue());

    } else if (satTable = true){
       
        // next, try to get the saturation limit for the chunk from a lookup table
       
        satLimit = satTableName[ccdNum];
                  
    } else if (useDefSat = true){

        // use a generic saturation limit from the policy file
        satLimit = datasetPolicy.getDouble("satLimit");       
    } else {
        // can't find the saturation limit for the chunk anywhere, I give up!
        throw lsst::pex::exceptions::NotFound(std::string("Can not get Saturation limit for Chunk Exposure."));
    }  

    // Get the bad pixel mask and setup the "SAT" mask plane
    typename lsst::afw::image::Mask<MaskT>::MaskPtrT chunkMaskPtr;
    chunkMaskPtr = chunkMaskedImage.getMask();
    MaskT const satMaskBit = chunkMaskedImage.getMask()->getPlaneBitMask("SAT");
    
    typedef std::vector<lsst::detection::Footprint::PtrType> FootprintList;
    typedef typename FootprintList::iterator FootprintIter;

    FootprintList newSatFps;   // newly detected saturated pixel footprints
    FootprintIter satFpIter;   // if need to iterate
    FootprintList grownSatFps; // final list of grown saturated pixel footprints

    // Save the saturated pixels as a vector of footprints.  

    // QQQ: Is there anyway of returning the number of pixels in each footprint?

    lsst::detection::DetectionSet<ImageT, MaskT> detectionSet(chunkMaskedImage, lsst::detection::Threshold(threshold, lsst::detection::Threshold::VALUE)
        );       

    newSatFps = detectionSet.getFootprints();
  
    // Grow around all of the saturated pixel footprints.  

    // QQQ: Can we distinguish between pixels saturated in the A/D converter and
    // those just saurated on chip?  If so, we don't want to grow around the A/D
    // saturated pixels (in unbinned data).

    int numSatFootprints = 0;
    int numSatPix = 0;
    for (FootprintIter satFpIter = newSatFps.begin(); satFpIter < newSatFps.end(); satFpIter++) { 

        // Need to create a bounding box to turn the saturated footprints into
        // new grown footprints
 
        vw::BBox2i const & bbox = (*satFpIter)->getBBox(); 
        vw::Vector2i const minVec(bbox.min().x() - satGrow, bbox.min().y() - satGrow); 
        vw::Vector2i const maxVec(bbox.max().x() + satGrow, bbox.max().y() + satGrow); 
        vw::BBox2i const fpBBox(minVec, maxVec); 
              
        // lets turn each into a subImage and get the cols/rows and number of
        // pixels in each footprint so we can sum them
        typename lsst::afw::image::MaskedImage<ImageT, MaskT>::MaskedImagePtrT fpChunkMaskedImagePtr;
        fpChunkMaskedImagePtr = chunkMaskedImage.getSubImage(fpBBox); 
        const int numCols = static_cast<int>(fpChunkMaskedImagePtr->getCols());
        const int numRows = static_cast<int>(fpChunkMaskedImagePtr->getRows());
        numSatPix += (numCols * numRows);
               
        // Create a new footprint with the grown bbox and save the new
        // footprints in another vector.

        lsst::detection::Footprint::PtrType fpGrow(new lsst::detection::Footprint(fpBBox)); 
        grownSatFps.push_back(fpGrow);
        numSatFootprints += 1;
    } 

    // Mask all of those saturated pixel footprints.  Using "SAT" bitmask for
    // all pixels in the footprint.  

    // QQQ: Do we want to distinguish between grown pixels and those that were
    // actually saturated?  What bitmask would that be ("GROW")? Detection set
    // will mask pixels but its not yet implemented.
   
    lsst::detection::setMaskFromFootprintList<MaskT>(chunkMaskPtr, grownSatFps, satMaskBit);
    
    // Interpolate over all of the masked saturated pixels. 
    lsst::ip::isr::interpolateOverMaskedPixels<ImageT, MaskT>(chunkExposure, isrPolicy);

    // Record the sub-stage provenance to the Image Metadata
    chunkMetadata->addProperty(lsst::daf::base::DataProperty("ISR_SATCOR")); 
    lsst::daf::base::DataProperty::PtrType satCorProp = chunkMetadata->findUnique("ISR_SATCOR");
    std::string exTrue = "Completed Successfully"; 
    satCorProp->setValue(boost::any_cast<std::string>(exTrue));
    chunkMetadata->addProperty(lsst::daf::base::DataProperty("ISR_SATCOR_END"));
    chunkMaskedImage.setMetadata(chunkMetadata);

    // Calculate additional SDQA Metrics??
       //return the following fro SDQA:
       // numSatFootprints
       // numSatPix

    // Issue a logging message if the sub-stage executes without issue to this
    // point! Yay!!

    lsst::pex::logging::TTrace<7>("ISR sub-stage, %s, completed successfully.", subStage);
         
}

/************************************************************************/
/* Explicit instantiations */

template
lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> saturationCorrectionForChunkExposure(
    lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> &chunkExposure,    
    lsst::pex::policy::Policy &isrPolicy, 
    lsst::pex::policy::Policy &datasetPolicy
    //std::vector<float> &saturationLookUpTable
    );

template
lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> saturationCorrectionForChunkExposure(
    lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> &chunkExposure,    
    lsst::pex::policy::Policy &isrPolicy, 
    lsst::pex::policy::Policy &datasetPolicy
    //std::vector<double> &saturationLookUpTable
    );

/************************************************************************/
