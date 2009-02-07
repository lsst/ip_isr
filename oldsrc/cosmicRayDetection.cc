// -*- LSST-C++ -*- // fixed format comment for emacs
/**
  * \file
  *
  * \ingroup isr
  *
  * \brief Implementation of the templated subStage, Cosmic Ray Detection, of
  * the Instrument Signature Removal stage forthe nightly LSST Image Processing
  * Pipeline.
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

#include <lsst/afw/image.h>
#include <lsst/daf/base/DataProperty.h>
#include <lsst/pex/exceptions/Exception.h>
#include <lsst/pex/logging/Trace.h>
#include <lsst/pex/policy/Policy.h>

#include "lsst/ip/isr/isr.h"

typedef double vectorType;
typedef double funcType;

/** \brief Correct the Chunk Exposure for cosmic ray hits.  This is a simple
  * first pass using a crude PSF.
  *
  * \return chunkExposure with cosmic rays removed.
  *
  * \throw Runtime if this sub-stage has been run previously on the image
  * \throw NotFound if any metadata parameter can not be obtained
  * 
  * TO DO (as of ):
  * - Calculate SDQA metrics as requested by SDQA team
  */

template<typename ImageT, typename MaskT>
void lsst::ip::isr::cosmicRayDetection(
    lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure,    
    lsst::pex::policy::Policy &isrPolicy, 
    lsst::pex::policy::Policy &datasetPolicy
    ) { 


}

/************************************************************************/
/* Explicit instantiations */

template
void lsst::ip::isr::cosmicRayDetection(
    lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> &chunkExposure,    
    lsst::pex::policy::Policy &isrPolicy, 
    lsst::pex::policy::Policy &datasetPolicy
    );

template
void lsst::ip::isr::cosmicRayDetection(
    lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> &chunkExposure,    
    lsst::pex::policy::Policy &isrPolicy, 
    lsst::pex::policy::Policy &datasetPolicy
    );

/************************************************************************/

