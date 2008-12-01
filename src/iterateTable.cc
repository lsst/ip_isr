
// -*- LSST-C++ -*- // fixed format comment for emacs
/**
  * \file
  *
  * \ingroup isr
  *
  * \brief Utility function to iterate through a table and apply the new value
  *        to the image.
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

// Given the maskedImage and the lookupTable, 
//\return the maskedImage adjusted to the table values

typedef double vectorType;

template<typename ImageT, typename MaskT> 
void lsst::ip::isr::iterateTable(
    lsst::afw::image::MaskedImage<ImageT, MaskT> &maskedImage,
    std::vector<vectorType> &lookupTable
    ) {

//    std::vector<vectorType>::iterator tableIter = lookupTable.begin();
//    std::vector<vectorType>::iterator tableIterEnd = lookupTable.end();

    int maxInd = lookupTable.size() - 1;
    const int miCols = static_cast<int>(maskedImage.getCols());
    const int miRows = static_cast<int>(maskedImage.getRows());

    lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> miRowAcc(maskedImage);
                
    for (int miRow = 0; miRow < miRows; miRow++, miRowAcc.nextRow()) {
        lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> miColAcc = miRowAcc;
        for (int miCol = 0; miCol < miCols; miCol++, miColAcc.nextCol()) {
            int ind = static_cast<int>(*miColAcc.image);
            if ((ind < 0) || (ind > maxInd)) {
            	 throw lsst::pex::exceptions::RangeError(std::string("Lookup table iterator out of range."));        
            } else {
                *miColAcc.image += lookupTable[ind]; 
	    } 
        }
    }
}

/************************************************************************/
/* Explicit instantiations */

template
void lsst::ip::isr::iterateTable(
    lsst::afw::image::MaskedImage<float, lsst::afw::image::maskPixelType> &maskedImage,
    std::vector<vectorType> &lookupTable
    );

template
void lsst::ip::isr::iterateTable(
    lsst::afw::image::MaskedImage<double, lsst::afw::image::maskPixelType> &maskedImage,
    std::vector<vectorType> &lookupTable
    );

/************************************************************************/
