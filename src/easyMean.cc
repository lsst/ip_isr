// -*- LSST-C++ -*- // fixed format comment for emacs
/**
  * \file
  *
  * \ingroup isr
  *
  * \brief Utility function to calculate the mean and standard deviation of an
  * image region.
  *
  * \author Nicole M. Silvestri, University of Washington
  *
  * Contact: nms@astro.washington.edu
  *
  * \version
  *
  * LSST Legalese here...
  */

#include <cmath>

#include <lsst/afw/image.h>
#include <lsst/afw/math.h>
#include "lsst/ip/isr/isr.h"

// Given the maskedImage, compute the number of elements (n), mean (mu) and
// standard deviation (sigma) borrowing some code from Russell Owen's
// medianBinApprox.
// \return mu, and sigma

typedef double funcType;

template <typename ImageT, typename MaskT>
double lsst::ip::isr::easyMean(
	lsst::afw::image::MaskedImage<ImageT, MaskT> &maskedImage
    ) {

    lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> miRowAcc(maskedImage);

    const int miCols = static_cast<int>(maskedImage.getCols());
    const int miRows = static_cast<int>(maskedImage.getRows());
    long int n = 0;
    double sum = 0;
       
    for (int miRow = 0; miRow < miRows; miRow++, miRowAcc.nextRow()) {
        lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> miColAcc = miRowAcc;
        for (int miCol = 0; miCol < miCols; miCol++, miColAcc.nextCol()) {
            n++;
            sum += static_cast<double>(*miColAcc.image);                 
        } 
    }     
    // the mean
    double mu = sum/static_cast<double>(n);
    
    return mu;

    double sumSq = 0;
    for (int miRow = 0; miRow < miRows; miRow++, miRowAcc.nextRow()) {
        lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> miColAcc = miRowAcc;
        for (int miCol = 0; miCol < miCols; miCol++, miColAcc.nextCol()) {
            n++;
            double value = static_cast<double>(*miColAcc.image) - mu; 
            sumSq += value * value; 
        } 
    }
    // the standard deviation
    double sigma = std::sqrt(sumSq/static_cast<double>(n));

}

/************************************************************************/
/* Explicit instantiations */

template
double lsst::ip::isr::easyMean(
    lsst::afw::image::MaskedImage<float, lsst::afw::image::maskPixelType> &maskedImage
    );

template
double lsst::ip::isr::easyMean(
    lsst::afw::image::MaskedImage<double, lsst::afw::image::maskPixelType> &maskedImage
    ); 

/************************************************************************/
