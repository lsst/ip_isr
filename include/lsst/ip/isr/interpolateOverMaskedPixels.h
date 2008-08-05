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

#include <string>

#include <boost/shared_ptr.hpp>

#include <lsst/afw/image/Exposure.h>
#include <lsst/afw/image/Image.h>
#include <lsst/afw/image/Mask.h>
#include <lsst/afw/math/Function.h>
#include <lsst/daf/base/DataProperty.h>
#include <lsst/pex/policy/Policy.h>

/** \brief Interpolate over all masked pixels in a Chunk Exposure's image using
  * one of several interpolation methods.
  */

namespace lsst {
namespace ip {
namespace isr {

    typedef boost::uint16_t maskPixelType;

    template<typename ImageT, typename MaskT>
    lsst::afw::image::MaskedImage<ImageT, MaskT> interpolateOverMaskedPixels(
        lsst::afw::image::MaskedImage<ImageT, MaskT> &chunkMaskedImage,
        lsst::daf::base::DataProperty::PtrType &chunkMetaData,
        std::string const interpMethod
        );

}}} // namespace lsst::ip::isr

#endif // !defined(LSST_IP_ISR_INTERPOLATEOVERMASKEDPIXELS_H)
