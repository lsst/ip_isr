// -*- LSST-C++ -*- // fixed format comment for emacs
/**
  * \file
  *
  * \ingroup isr
  *
  * \brief Implementation of the templated subStage, Illumination
  *  Correction, of the Instrument Signature Removal stage for the nightly 
  *  LSST Image Processing Pipeline.
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
  *
  * \return masterChunkExposure corrected for scattered light
  *
  * \throw Runtime if this sub-stage has been run previously
  * \throw NotFound if any Policy or metadata value can not be obtained
  * \throw InvalidParameter if functional form for the lineaization fit is invalid
  */
template<typename ImageT, typename MaskT>
lsst::afw::image::Exposure<ImageT, MaskT> illuminationCorrection(
    lsst::afw::image::Exposure<ImageT, MaskT> &masterChunkExposure, // the Master Flat Filed Exposure
    lsst::afw::image::Exposure<ImageT, MaskT> const &masterIcChunkExposure, //the Master Illumination Correction Exposure
    lsst::pex::policy::Policy &isrPolicy,
    lsst::pex::policy::Policy &datasetPolicy
    ) {


}

/************************************************************************/
/* Explicit instantiations */

template
lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> illuminationCorrection(
    lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> &masterChunkExposure,
    lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> const &masterIcChunkExposure,
    lsst::pex::policy::Policy &isrPolicy,
    lsst::pex::policy::Policy &datasetPolicy
    );

template
lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> illuminationCorrection(
    lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> &masterChunkExposure,
    lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> const &masterIcChunkExposure,
    lsst::pex::policy::Policy &isrPolicy,
    lsst::pex::policy::Policy &datasetPolicy
    );

/************************************************************************/
