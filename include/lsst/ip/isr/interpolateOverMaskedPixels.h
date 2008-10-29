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

#include <lsst/afw/image/Exposure.h>
#include <lsst/pex/policy/Policy.h>

/** \brief Interpolate over all masked pixels in a Chunk Exposure's image using
  * one of several interpolation methods.
  */

namespace lsst {
namespace ip {
namespace isr {

    template<typename ImageT, typename MaskT>
    lsst::afw::image::Exposure<ImageT, MaskT> interpolateOverMaskedPixels(
        lsst::afw::image::Exposure<ImageT, MaskT> const &chunkExposure,
	lsst::pex::policy::Policy &isrPolicy
        );

}}} // namespace lsst::ip::isr

#endif // !defined(LSST_IP_ISR_INTERPOLATEOVERMASKEDPIXELS_H)
