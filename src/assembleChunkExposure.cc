// -*- LSST-C++ -*- // fixed format comment for emacs
/**
  * \file
  *
  * \ingroup isr
  *
  * \brief Implementation of the templated subStage, Assemble Chunk
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

#include <boost/shared_ptr.hpp>
#include <boost/cstdint.hpp>
#include <boost/format.hpp>

#include <lsst/afw/image/Exposure.h>
#include <lsst/afw/image/Mask.h>
#include <lsst/afw/image/MaskedImage.h>
#include <lsst/afw/math/Function.h>
#include <lsst/daf/base/DataProperty.h>
#include <lsst/pex/exceptions/Exception.h>
#include <lsst/pex/logging/Trace.h>
#include <lsst/pex/policy/Policy.h>

#include <lsst/ip/isr/isr.h>
#include <lsst/ip/isr/interpolateOverMaskedPixels.h>

/** \brief Create a Chunk Exposure using the crosstalk-corrected Chunk Image
  * from the telescope, the Chunk Variance Image (synthesized from the pixel
  * values, gain and rdnoise), and the Chunk Mask (created using a copy of
  * initial Bad Pixel Mask from the telescope).  The Chunk Exposure is created
  * with no WCS.
  */

template<typename ImageT, typename MaskT>
lsst::afw::image::Exposure<ImageT, MaskT> assembleChunkExposure(
    lsst::afw::image::Image<ImageT> &chunkImage,
    lsst::daf::base::DataProperty::PtrType &chunkMetaData,
    lsst::afw::image::Mask<MaskT> &badPixelMask, 
    std::string const interpMethod
    ) { 

    double gain;
    double rdNoise;
    const int numCols = static_cast<int>(chunkImage.getCols());
    const int numRows = static_cast<int>(chunkImage.getRows());    

    // make copy of bad pixel mask to use as the ChunkMask.
    std::map<std::string, int> maskPlaneDict = badPixelMask->getMaskPlaneDict();  
    lsst::afw::image::Mask<MaskT> chunkMask(numCols, numRows, maskPlaneDict);
    *chunkMask = *badPixelMask;  

    // if we are not given a badPixelMask, will have to grab them from the Image itself??
    // int badMaskBit = chunkImage.getMask()->getMaskPlane("BAD");
    // badPixelMask = (badMaskBit < 0) ? 0 : (1 << badMaskBit);

    // Calculate the variance per pixel (in ADU**2). Assuming the units of the
    // pixels are ADUs, gain are electrons/ADU, and rdNoise are electrons.
    // Note: this is different from the variance calculation in the MaskedImage
    // class.

    //get the gain from the metadata
    lsst::daf::base::DataProperty::PtrType gainField = chunkMetaData->findUnique("GAIN");
    if (gainField) {
        gain = boost::any_cast<const double>(gainField->getValue());
        return gain;
    }
    throw lsst::pex::exceptions::Runtime(std::string("in ") + __func__ + std::string(": Could not get GAIN from the chunkMetaData"));

    // get the read noise from the metadata
    lsst::daf::base::DataProperty::PtrType rdNoiseField = chunkMetaData->findUnique("RDNOISE");
    if (rdNoiseField) {
        rdNoise = boost::any_cast<const double>(rdNoiseField->getValue());
        return rdNoise;
    }
    throw lsst::pex::exceptions::Runtime(std::string("in ") + __func__ + std::string(": Could not get RDNOISE from the chunkMetaData"));

    // Create an empty Chunk Variance Image that is the same size as the
    // chunkImage and populate with new pixel values.
    lsst::afw::image::Image<ImageT> chunkVarianceImage(numCols, numRows);   
    *chunkVarianceImage = *chunkImage/gain + 1/pow(gain,2) * pow(rdNoise,2);
  
    // Create the MaskedImage from the Chunk Image, Chunk VarianceImage, and the Chunk Mask
    lsst::afw::image::MaskedImage<ImageT, MaskT> chunkMaskedImage(*chunkImage, *chunkVarianceImage, *chunkMask);
    
    // Interpolate over masked pixels in the new MaskedImage.
    lsst::afw::image::MaskedImage<ImageT, MaskT> finalChunkMaskedImage = lsst::ip::isr::interpolateOverMaskedPixels<ImageT, MaskT>(chunkMaskedImage, chunkMetaData, interpMethod); 

    // Create the Chunk Exopsure
    lsst::afw::image::Exposure<ImageT, MaskT> chunkExposure(finalChunkMaskedImage);
    return chunkExposure;
}

/***************************************************************************/
/* Explicit instantiations */

