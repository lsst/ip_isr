// -*- LSST-C++ -*- // fixed format comment for emacs
/**
  * \file
  *
  * \ingroup isr
  *
  * \brief Implementation of the templated subStage, Saturation Correction for
  * Chunk Exposure, of the Instrument Signature Removal stage for the nightly
  * LSST Image Processing Pipeline.
  *
  * \author Nicole M. Silvestri, University of Washington
  */
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>

#include <boost/shared_ptr.hpp>
#include <boost/cstdint.hpp>
#include <boost/format.hpp>

#include <lsst/afw/image/Exposure.h>
#include <lsst/afw/math/Function.h>
#include <lsst/afw/image/Mask.h>
#include <lsst/afw/image/MaskedImage.h>
#include <lsst/afw/image/PixelAccessors.h>
#include <lsst/daf/base/DataProperty.h>
#include <lsst/pex/exceptions/Exception.h>
#include <lsst/pex/logging/Trace.h>
#include <lsst/pex/policy/Policy.h>

#include "lsst/ip/isr/isr.h"
#include "lsst/ip/isr/interpolateOverMaskedPixels.h"

/** \brief Detect and mask pixels that are saturated (in the A/D converter or
  * have excessive non-linear response) in the Chunk Exposure.  
  *
  * Grow by additional pixels (number of pixels TBA) to mask
  * charge spillover.  Set appropriate bits in the Mask and Chunk Exposure and
  * update the metadata for the Chunk Exposure. 
  */

template<typename ImageT, typename MaskT>
lsst::afw::image::Exposure<ImageT, MaskT> saturationCorrectionForChunkExposure(
    lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure,    
    lsst::daf::base::DataProperty::PtrType &chunkMetaData, 
    lsst::pex::policy::Policy &policy
    // lsst::daf::base::DataProperty::PtrType &saturationLookUpTable
    ) { 

    double satLimit;
    // Get saturation correction sub-stage information from the Policy 
    // When a pixel is saturated, all pixels within this distance of the masked
    // saturated pixel are also masked */
    double grow = policy.getDouble("staurationCorrectionForChunkExposure.grow");
    // Get the saturation limit from the Chunk MetaData - what is the proper
    // MetaData keyword for this?
    lsst::daf::base::DataProperty::PtrType saturationField = chunkMetaData->findUnique("SATURATE");
    if (saturationField) {
        double satLimit = boost::any_cast<const double>(saturationField->getValue());
        return satLimit;
    } else {
        throw lsst::pex::exceptions::Runtime(std::string("in ") + __func__ + std::string(": Could not get SATURATE from the chunkMetaData"));

        // default saturation limit - to be used if the amplifier saturation
        // limit is not present in the Chunk MetaData - will assume this is
        // policy file for now but may be a lookup table.
        satLimit = policy.getDouble("staurationCorrectionForChunkExposure.satLimit");
    }

    /* Get saturatoin mask bit for the chunkExposure */
    int saturationMaskBit = chunkExposure.getMask()->getMaskPlane("SATURATE");

    const int numCols = static_cast<int>(chunkExposure.getCols());
    const int numRows = static_cast<int>(chunkExposure.getRows()); 

    lsst::afw::image::MaskedImage<ImageT, MaskT> chunkMi = chunkExposure.getMaskedImage();
    lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> chunkRowAcc(chunkMi);

    /* for each pixel in the chunkExposure */
   
    for (int chunkRow = 0; chunkRow < numRows; chunkRow++, chunkRowAcc.nextRow()) {
        lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> chunkColAcc = chunkRowAcc;
        for (int chunkCol = 0; chunkCol < numCols; chunkCol++, chunkColAcc.nextCol()) {
            if (chunkCol >= satLimit){
                
                *chunkColAcc.image = 0;
                *chunkColAcc.variance = 0;
                *chunkColAcc.mask = saturationMaskBit;
            }

            //   if (
                // need to grow around saturated pixel - not the pixels
                // saturated in the A/D conv.  Then need to mask the grow'n
                // pixels.
            //       ){

                //use vw::BBox2i or detection's grow/footprint method...
            // }

            //then need to interpolate over these masked pixels. Call
            //"InterpolateOverMaskedPixels" utility function.
            
        } // for column loop
    } // for row loop     
           
}

/************************************************************************/
/* Explicit instantiations */
