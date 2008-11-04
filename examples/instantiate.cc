// -*- LSST-C++ -*- // fixed format comment for emacs
/**
  * \file
  *
  * \ingroup isr
  *
  * \brief Simple test code for explicit instantiations
  *
  * \author Nicole M. Silvestri, University of Washington
  *
  * Contact: nms@astro.washington.edu
  * 
  */
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>

#include <boost/cstdint.hpp>
#include <boost/format.hpp>
#include <boost/operators.hpp>
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

template<typename ImageT, typename MaskT>
lsst::afw::image::Exposure<ImageT, MaskT> saturationCorrectionForChunkExposure(
    lsst::afw::image::Exposure<ImageT, MaskT> const &chunkExposure,    
    lsst::pex::policy::Policy &isrPolicy, 
    lsst::pex::policy::Policy &datasetPolicy
    //std::vector<float> &saturationLookUpTable
    ) { 

    lsst::afw::image::MaskedImage<ImageT, MaskT> chunkMaskedImage = chunkExposure.getMaskedImage();
    lsst::pex::policy::Policy::Ptr saturationPolicy = isrPolicy.getPolicy("saturationPolicy");
   

}

// explicit instantiation

template
lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> saturationCorrectionForChunkExposure(
    lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> const &chunkExposure,    
    lsst::pex::policy::Policy &isrPolicy, 
    lsst::pex::policy::Policy &datasetPolicy
    //std::vector<float> &saturationLookUpTable
    );
