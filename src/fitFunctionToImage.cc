
// -*- LSST-C++ -*- // fixed format comment for emacs
/**
  * \file
  *
  * \ingroup isr
  *
  * \brief Utility function to fit a polynomial function to a maskedImage.
  *
  * \author Nicole M. Silvestri, University of Washington
  *
  * Contact: nms@astro.washington.edu
  *
  * \version
  *
  * LSST Legalese here...
  */

#include <lsst/afw/image.h>
#include <lsst/afw/math.h>
#include "lsst/ip/isr/isr.h"

//Given the maskedImage and the polynomial function, 
//\return the maskedImage fitted for the polynomial.

typedef double funcType;

template<typename ImageT, typename MaskT> 
void lsst::ip::isr::fitFunctionToImage(
    lsst::afw::image::MaskedImage<ImageT, MaskT> const &maskedImage,
    lsst::afw::math::PolynomialFunction1<funcType> const &polyFunction
    ) {

    const int miCols = static_cast<int>(maskedImage.getCols());
    const int miRows = static_cast<int>(maskedImage.getRows());

    lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> miRowAcc(maskedImage);
                
    for (int miRow = 0; miRow < miRows; miRow++, miRowAcc.nextRow()) {
        lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> miColAcc = miRowAcc;
        for (int miCol = 0; miCol < miCols; miCol++, miColAcc.nextCol()) {
            *miColAcc.image = static_cast<ImageT>(polyFunction(*miColAcc.image));
            *miColAcc.variance = static_cast<ImageT>(polyFunction(*miColAcc.image) * polyFunction(*miColAcc.image));
        }
    }
}

/************************************************************************/
/* Explicit instantiations */

template
void lsst::ip::isr::fitFunctionToImage(
    lsst::afw::image::MaskedImage<float, lsst::afw::image::maskPixelType> const &maskedImage,
    lsst::afw::math::PolynomialFunction1<funcType> const &polyFunction
    );

template
void lsst::ip::isr::fitFunctionToImage(
    lsst::afw::image::MaskedImage<double, lsst::afw::image::maskPixelType> const &maskedImage,
    lsst::afw::math::PolynomialFunction1<funcType> const &polyFunction
    );

/************************************************************************/
