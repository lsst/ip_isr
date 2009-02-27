
// -*- LSST-C++ -*- // fixed format comment for emacs
/**
  * \file
  *
  * \ingroup isr
  *
  * \brief Utility function to fit a polynomial function to a image.
  *
  * \author Nicole M. Silvestri, University of Washington
  *
  * Contact: nms@astro.washington.edu
  *
  * \version
  *
  * LSST Legalese here...
  */

#include "lsst/ip/isr/isr.h"

//Given the image and the polynomial function, 
//\return the image fitted for the polynomial.

template<typename ImagePixelT>
void lsst::ip::isr::fitFunctionToImage(
    lsst::afw::image::MaskedImage<ImagePixelT> &maskedImage,
    lsst::afw::math::Function1<double> const &function
) {
    typedef lsst::afw::image::MaskedImage<ImagePixelT> MaskedImage;

    const int miHeight = maskedImage.getHeight();

    // Set the pixels row by row, to avoid repeated checks for end-of-row
    for (int y = 0; y < miHeight; ++y) {
        for (typename MaskedImage::x_iterator miPtr = maskedImage.row_begin(y), end = maskedImage.row_end(y); miPtr != end; ++miPtr) {
            miPtr.image() = static_cast<ImagePixelT>(function(static_cast<double>(miPtr.image())));
        }
    }
}

/************************************************************************/
/* Explicit instantiations */

template
void lsst::ip::isr::fitFunctionToImage(
    lsst::afw::image::MaskedImage<float> &maskedImage,
    lsst::afw::math::Function1<double> const &function
    );

template
void lsst::ip::isr::fitFunctionToImage(
    lsst::afw::image::MaskedImage<double> &maskedImage,
    lsst::afw::math::Function1<double> const &function
    );

/************************************************************************/
