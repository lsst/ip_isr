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

#include <lsst/afw/image/image.h>
#include <lsst/afw/math.math.h>
#include "lsst/ip/isr/isr.h"

// Given the maskedImage, compute the number of elements (n), mean (mu) and
// standard deviation (sigma) borrowing some code from Russell Owen's
// medianBinApprox.
// \return n, mu, and sigma

typedef double funcType;

int lsst::ip::isr::easyMean(long int &n, double &mu, double &sigma, lsst::afw::image::MaskedImage<ImageT, MaskT> &maskedImage
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
