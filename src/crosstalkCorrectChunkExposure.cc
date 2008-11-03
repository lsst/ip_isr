// -*- LSST-C++ -*- // fixed format comment for emacs
/**
  * \file
  *
  * \ingroup isr
  *
  * \brief Implementation of the templated subStage, Crosstalk Correct Chunk
  * Exposure, of the Instrument Signature Removal stage forthe nightly LSST
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

/** \brief Correct the Chunk Exposure for crosstalk contributions from other
  * Chunks.
  *
  * \return chunkExposure with crosstalk correction
  *
  * \throw Runtime if this sub-stage has been run previously on the image
  * \throw NotFound if any metadata parameter can not be obtained
  * 
  * TO DO (as of ):
  * - Calculate SDQA metrics as requested by SDQA team
  */

template<typename ImageT, typename MaskT>
lsst::afw::image::Exposure<ImageT, MaskT> crosstalkCorrectChunkExposure(
    lsst::afw::image::Exposure<ImageT, MaskT> const &chunkExposure,    
    lsst::pex::policy::Policy &isrPolicy, 
    lsst::pex::policy::Policy &datasetPolicy
    //std::vector<float> &crosstalkLookUpTable
    ) { 


}

/************************************************************************/
/* Explicit instantiations */

template
lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> crosstalkCorrectChunkExposure(
    lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> const &chunkExposure,    
    lsst::pex::policy::Policy &isrPolicy, 
    lsst::pex::policy::Policy &datasetPolicy
    //std::vector<float> &crosstalkLookUpTable
    );

template
lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> crosstalkCorrectChunkExposure(
    lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> const &chunkExposure,    
    lsst::pex::policy::Policy &isrPolicy, 
    lsst::pex::policy::Policy &datasetPolicy
    //std::vector<double> &crosstalkLookUpTable
    );

/************************************************************************/
