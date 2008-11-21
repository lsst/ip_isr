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
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

#include "boost/cstdint.hpp"
#include "boost/format.hpp"
#include "boost/shared_ptr.hpp"
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
  * TO DO (as of Tue 11/11/08):
  * - delineate between A/D saturated pixels and other?
  * - Calculate additional SDQA metrics as requested by SDQA team
  * - use threshold or satLimit (from LookupTable?)??
  *
  */

typedef double vectorType;
std::string satStage = "lsst.ip.isr.saturationCorrectionForChunkExposure";

template<typename ImageT, typename MaskT>
void lsst::ip::isr::saturationCorrectionForChunkExposure(
    lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure,    
    lsst::pex::policy::Policy &isrPolicy, 
    lsst::pex::policy::Policy &datasetPolicy
    //std::vector<vectorType> &saturationLookUpTable
    ) { 

    // Start with some setup and information gathering...
    // Get the Chunk MaskedImage and Chunk Metadata from the Chunk Exposure

    lsst::pex::logging::TTrace<3>("Entering ISR stage: %s", satStage);

    lsst::afw::image::MaskedImage<ImageT, MaskT> chunkMaskedImage = chunkExposure.getMaskedImage();    

    typename lsst::afw::image::MaskedImage<ImageT, MaskT>::MaskedImagePtrT chunkMaskedImagePtr; 

    lsst::daf::base::DataProperty::PtrType chunkMetadata = chunkMaskedImage.getImage()->getMetaData();

    // Make sure we haven't run this sub-stage previously

     lsst::daf::base::DataProperty::PtrType isrSatFlag = chunkMetadata->findUnique("SATU_END");  // ISR saturation stage processing flag
    if (isrSatFlag) {
        throw lsst::pex::exceptions::Runtime(std::string("Saturation Correction previously performed - terminating stage."));
    } 
    
    // Parse the ISR policy file for the saturation sub-stage information
    // First get the Saturation Stage Policy from the ISR Policy
 
    lsst::pex::policy::Policy::Ptr saturationPolicy; //QQQ: shouldn't this be 'PtrType' in Policy.h??
    bool satTable;
    bool useDefSat;
    std::string satTableName;
    int satGrow;
    double threshold;
    // const int nSatPixMin;
    try {
        std::cout << "Attempting to parse the Policy File." << std::endl;
        saturationPolicy = isrPolicy.getPolicy("saturationPolicy");
        satTable = saturationPolicy->getBool("satTable");
        useDefSat = saturationPolicy->getBool("useDefSat");
        satTableName = saturationPolicy->getString("satTableName");
        std::string satTableName = saturationPolicy->getString("satTableName"); 
        satGrow = saturationPolicy->getInt("grow");
        threshold = saturationPolicy->getDouble("threshold");
        // nSatPixMin = saturationPolicy->getInt("nSatPixMin");
    } catch (std::exception e) {
        throw lsst::pex::exceptions::Runtime(std::string("Can not parse the main ISR policy file."));
    }

    std::string ccdName;
    try {
        ccdName = datasetPolicy.getString("ccdName");
    } catch (std::exception e) {
         throw lsst::pex::exceptions::Runtime(std::string("Can not parse the dataset policy file."));
    }

    // Get the saturation limit for the Chunk Exposure...QQQ: may not need
    // this if we use a threshold??

    int satLimit;  
    std::vector<vectorType> satTableList;
     // First, try to get the saturation limit for the chunk from the
     // chunkMetadata
    lsst::daf::base::DataProperty::PtrType saturationField = chunkMetadata->findUnique("SATURATE");
    if (saturationField) {
       
        satLimit = boost::any_cast<const int>(saturationField->getValue());

    } else if (satTable = true){
       
        // Next, try to get the saturation limit for the chunk from a lookup
        // table.  

        // QQQ: what is the form of this lookup table (chunk#
        // satLimit)??
       
//         std::ifstream in(satTableName);
//    const int numMiCols = static_cast<int>(chunkMaskedImage.getCols());
//    const int numMiRows = static_cast<int>(chunkMaskedImage.getRows());  
//         int num = numMiCols * numMiRows;
//         double ccdnum;
//         while(in >> ccdnum)
//            satTableList.push_back(ccdnum);
//         satLimit = satTableName[ccdNum];
                  
    } else if (useDefSat = true){

        // use a generic saturation limit from the policy file
        satLimit = datasetPolicy.getInt("satLimit");       
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

    lsst::detection::DetectionSet<ImageT, MaskT> detectionSet(chunkMaskedImage, lsst::detection::Threshold(threshold, lsst::detection::Threshold::VALUE));       

    newSatFps = detectionSet.getFootprints();

     // Grow around all of the saturated pixel footprints.  

    lsst::detection::DetectionSet<ImageT, MaskT> detectionSetGrown(newSatFps, grow);

    grownSatFps = detectionSetGrown.getFootprints();

    // QQQ: Can we distinguish between pixels saturated in the A/D converter and
    // those just saurated on chip?  If so, we don't want to grow around the A/D
    // saturated pixels (in unbinned data).

//     int numSatFootprints = 0;
//     int numSatPix = 0;
//     for (FootprintIter satFpIter = newSatFps.begin(); satFpIter < newSatFps.end(); satFpIter++) { 

//         // Need to create a bounding box to turn the saturated footprints into
//         // new grown footprints
 
//         vw::BBox2i const & bbox = (*satFpIter)->getBBox(); 
//         vw::Vector2i const minVec(bbox.min().x() - satGrow, bbox.min().y() - satGrow); 
//         vw::Vector2i const maxVec(bbox.max().x() + satGrow, bbox.max().y() + satGrow); 
//         vw::BBox2i const fpBBox(minVec, maxVec); 
              
//         // lets turn each into a subImage and get the cols/rows and number of
//         // pixels in each footprint so we can sum them
//         typename lsst::afw::image::MaskedImage<ImageT, MaskT>::MaskedImagePtrT fpChunkMaskedImagePtr;
        
//         // 'getSubImage will throw an exception if the requested subImage is
//         // outside of the image.  This happens when we grow a footprint that is
//         // close to the edge of the chunk.  QQQ: Need to get the EDGE bit and do
//         // something here to deal with this case better.  Catch these for now,
//         // log them and continue.

//         try {
//         fpChunkMaskedImagePtr = chunkMaskedImage.getSubImage(fpBBox); 
//         const int numCols = static_cast<int>(fpChunkMaskedImagePtr->getCols());
//         const int numRows = static_cast<int>(fpChunkMaskedImagePtr->getRows());
//         // use Footprint::setNpix?
//         numSatPix += (numCols * numRows);
//         } catch (lsst::pex::exceptions::ExceptionStack &e) {
//             lsst::pex::logging::TTrace<3>("In ISR stage %s, Requested footprint BBox, %d, is not contained within the original Image.", satStage, numSatFootprints + 1);
//             continue;
//         }
//         // Create a new footprint with the grown bbox and save the new
//         // footprints in another vector.

//         lsst::detection::Footprint::PtrType fpGrow(new lsst::detection::Footprint(fpBBox)); 
//         grownSatFps.push_back(fpGrow);
//         numSatFootprints += 1;
//     } 

    // Mask all of those saturated pixel footprints.  Using "SAT" bitmask for
    // all pixels in the footprint.  

    // QQQ: Do we want to distinguish between grown pixels and those that were
    // actually saturated?  What bitmask would that be ("GROW")? Detection set
    // will mask pixels but its not yet implemented.
   
    lsst::detection::setMaskFromFootprintList<MaskT>(chunkMaskPtr, grownSatFps, satMaskBit);
    
    // Interpolate over all of the masked saturated pixels. 
    lsst::ip::isr::interpolateOverMaskedPixels(chunkExposure, isrPolicy);

    // Record the sub-stage provenance to the Image Metadata

    lsst::daf::base::DataProperty::PtrType isrThresholdFlag(new lsst::daf::base::DataProperty("SATU_TH", threshold));
    chunkMetadata->addProperty(isrThresholdFlag);
    lsst::daf::base::DataProperty::PtrType isrPixelFlag(new lsst::daf::base::DataProperty("SATU_PIX", numSatPix));
    chunkMetadata->addProperty(isrPixelFlag);
    lsst::daf::base::DataProperty::PtrType isrFootprintFlag(new lsst::daf::base::DataProperty("SATU_FP", numSatFootprints));
    chunkMetadata->addProperty(isrFootprintFlag);
    lsst::daf::base::DataProperty::PtrType isrEndFlag(new lsst::daf::base::DataProperty("SATU_END", std::string("Completed Successfully"))); 
    chunkMetadata->addProperty(isrEndFlag); 
    chunkMaskedImage.setMetadata(chunkMetadata);

    // Calculate additional SDQA Metrics??
       //return the following for SDQA:
       // numSatFootprints
       // numSatPix

    // Issue a logging message if the sub-stage executes without issue to this
    // point! Yay!!

    lsst::pex::logging::TTrace<3>("ISR stage: %s completed successfully.", satStage);
    lsst::pex::logging::TTrace<3>("Leaving ISR stage: %s", satStage);
         
}

/************************************************************************/
/* Explicit instantiations */

template
void lsst::ip::isr::saturationCorrectionForChunkExposure(
    lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> &chunkExposure,    
    lsst::pex::policy::Policy &isrPolicy, 
    lsst::pex::policy::Policy &datasetPolicy
    //std::vector<vectorType> &saturationLookUpTable
    );

template
void lsst::ip::isr::saturationCorrectionForChunkExposure(
    lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> &chunkExposure,    
    lsst::pex::policy::Policy &isrPolicy, 
    lsst::pex::policy::Policy &datasetPolicy
    //std::vector<vectorType> &saturationLookUpTable
    );

/************************************************************************/
