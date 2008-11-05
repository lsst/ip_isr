// -*- LSST-C++ -*- // fixed format comment for emacs
/**
  * \file
  *
  * \ingroup isr
  *
  * \brief Implementation of the templated utility function, Interpolate Over
  * Masked Pixels.
  *
  * \author Nicole M. Silvestri, University of Washington
  *
  * Contact: nms@astro.washington.edu
  *
  * \version
  *
  * LSST Legalese here...
  */

#ifndef LSST_IP_ISR_INTERPOLATEOVERMASKEDPIXELS_H
#define LSST_IP_ISR_INTERPOLATEOVERMASKEDPIXELS_H

#include <lsst/daf/base.h>
#include <lsst/daf/data/LsstBase.h>	
#include <lsst/afw/image.h>
#include <lsst/pex/policy/Policy.h>

/** \brief Interpolate over all masked pixels in a Chunk Exposure's image using
  * one of several interpolation methods.
  */

namespace lsst {
namespace ip {
namespace isr {

    template<typename ImageT, typename MaskT>
    lsst::afw::image::Exposure<ImageT, MaskT> interpolateOverMaskedPixels(
        lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure, // the chunk exposure with pixels to be interpolated over
	lsst::pex::policy::Policy &stagePolicy // the main policy file for the stage being executed
        );

}}} // namespace lsst::ip::isr

#endif // !defined(LSST_IP_ISR_INTERPOLATEOVERMASKEDPIXELS_H)
