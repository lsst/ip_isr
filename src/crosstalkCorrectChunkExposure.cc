// -*- LSST-C++ -*- // fixed format comment for emacs
/**
  * \file
  *
  * \ingroup isr
  *
  * \brief Implementation of the templated subStage, Trim Chunk
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

#include <boost/shared_ptr.hpp>
#include <boost/cstdint.hpp>
#include <boost/format.hpp>

#include <lsst/afw/Exposure.h>
#include <lsst/afw/Function.h>
#include <lsst/afw/Mask.h>
#include <lsst/afw/MaskedImage.h>
#include <lsst/daf/data/DataProperty.h>
#include <lsst/daf/exceptions/Exception.h>
#include <lsst/daf/utils/Trace.h>
#include <lsst/pex/policy/Policy.h>

#include "lsst/ip/isr/isr.h"

/** \brief Remove overscan strip region and other non-illuminated edge pixels
  * from the Chunk Exposure.  Valid region to be retained is defined by the
  * DATASEC (defines the four corners of the valid region).
  */
