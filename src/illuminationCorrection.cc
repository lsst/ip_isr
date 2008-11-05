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

/** \brief Correct Master Flat Field Chunk Exposures for the differences in dome
  * vs. sky illumination.  This corrects for large scale structure in the system
  * response for large-scale surveys.  This method is for Data Release only
  * because it requires the median night sky flat produced from the science
  * images whcih obviously can not be done during nightly processing.
  *
  * \return masterChunkExposure corrected for scattered light
  *
  * \throw Runtime if this sub-stage has been run previously
  * \throw NotFound if any Policy or metadata value can not be obtained
  * \throw InvalidParameter if functional form for the lineaization fit is invalid
  * 
  * QQQ: filter dependence of ilumination corrections.
  * 
  * TO DO (as of  11//2008): 
  */

typedef double vectorType;
typedef double funcType;

template<typename ImageT, typename MaskT>
lsst::afw::image::Exposure<ImageT, MaskT> illuminationCorrection(
    lsst::afw::image::Exposure<ImageT, MaskT> &masterDfChunkExposure, // the Master Dome Flat Field Chunk Exposure
    lsst::afw::image::Exposure<ImageT, MaskT> &masterSfChunkExposure, // the Master Sky Flat Field Chunk Exposure
    lsst::pex::policy::Policy &isrPolicy,  // the main ISR Policy File containing the Illumination Policy info.
    lsst::pex::policy::Policy &datasetPolicy // policy file with info. specific to the dataset being processed
    ) {


 // estimate the illumination correction (assuming that theglobal changes in illumination are 
 // negligible (or very small).  Use the Master Dome Flat Field Chunk Exposure from the current 
 // night (Df(t2), the Master Dome Flat Field Chunk Exposure from the provious night (or one close 
 // in time; Df(t1)) and the illumination correction from the previous night (or one close in 
 // time; I(t2)) as follows:
 // I(t2) = smoothK(Df(t1)/Df(t2))*I(t2) where smooth K is the smoothing kernel  

}

/************************************************************************/
/* Explicit instantiations */

template
lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> illuminationCorrection(
    lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> &masterChunkExposure,
    lsst::pex::policy::Policy &isrPolicy,
    lsst::pex::policy::Policy &datasetPolicy
    );

// template
// lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> illuminationCorrection(
//     lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> &masterChunkExposure,
//     lsst::pex::policy::Policy &isrPolicy,
//     lsst::pex::policy::Policy &datasetPolicy
//     );

/************************************************************************/
