// -*- LSST-C++ -*- // fixed format comment for emacs
/**
  * \file
  *
  * \ingroup isr
  *
  * \brief Implementation of the templated subStage, Mask and Correct Additional
  * Artifacts, of the Instrument Signature Removal stage for the nightly LSST
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
#include <sstream>
#include <vector>
#include <cmath>

#include "boost/shared_ptr.hpp"
#include "boost/cstdint.hpp"
#include "boost/format.hpp"

#include <lsst/afw/image/Exposure.h>
#include <lsst/afw/math/Function.h>
#include <lsst/afw/image/Mask.h>
#include <lsst/afw/image/MaskedImage.h>
#include <lsst/daf/base/DataProperty.h>
#include <lsst/pex/exceptions/Exception.h>
#include <lsst/pex/logging/Trace.h>
#include <lsst/pex/policy/Policy.h>

#include "lsst/ip/isr/isr.h"

/** \brief 
  */


 template<typename ImageT, typename MaskT>
 lsst::afw::image::Exposure<ImageT, MaskT> maskAndCorrectAdditionalArtifacts(
     lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure,
     lsst::pex::policy::Policy &isrPolicy,
     lsst::pex::policy::Policy &datasetPolicy
     ){

}

/************************************************************************/
/* Explicit instantiations */

 template
 lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> maskAndCorrectAdditionalArtifacts(
     lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> &chunkExposure,
     lsst::pex::policy::Policy &isrPolicy,
     lsst::pex::policy::Policy &datasetPolicy
     );

 template
 lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> maskAndCorrectAdditionalArtifacts(
     lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> &chunkExposure,
     lsst::pex::policy::Policy &isrPolicy,
     lsst::pex::policy::Policy &datasetPolicy
     );

/************************************************************************/
