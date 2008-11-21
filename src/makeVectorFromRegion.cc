
// -*- LSST-C++ -*- // fixed format comment for emacs
/**
  * \file
  *
  * \ingroup isr
  *
  * \brief Utility function to create vectors from image regions.
  *
  * \author Nicole M. Silvestri, University of Washington
  *
  * Contact: nms@astro.washington.edu
  *
  * \version
  *
  * LSST Legalese here...
  */

#include <vector>

#include <lsst/afw/image/image.h>
#include "lsst/ip/isr/isr.h"

// Given two empty lists for holding the pixel values of the image and the
// variance, and two empty lists to hold the column and row numbers, and the
// maskedImage 
// \return the four populated lists.

typedef double vectorType;

template<typename vectorType> 
void lsst::ip::isr::makeVectorFromRegion(
    std::vector<vectorType> &list1, 
    std::vector<vectorType> &list2, 
    std::vector<vectorType> &list3, 
    std::vector<vectorType> &list4, 
    lsst::afw::image::MaskedImage<ImageT, MaskT> &maskedImage
    ) {

    lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> miRowAcc(maskedImage);
       
    const int miCols = static_cast<int>(maskedImage.getCols());
    const int miRows = static_cast<int>(maskedImage.getRows());

    for (int miRow = 0; miRow < miRows; miRow++, miRowAcc.nextRow()) {
        lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> miColAcc = miRowAcc;
        for (int miCol = 0; miCol < miCols; miCol++, miColAcc.nextCol()) {
            list1.push_back(*miColAcc.image);
            list2.push_back(*miColAcc.variance);
            list3.push_back(chunkCol);
            list4.push_back(chunkRow);
        }
    }   
}
