
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

#include <lsst/afw/image/image.h>
#include <lsst/afw/math.math.h>
#include "lsst/ip/isr/isr.h"

// Given the maskedImage and the lookupTable, 
//\return the maskedImage adjusted to the table values

typedef double vectorType;

template<typename ImageT, typename MaskT> 
void lsst::ip::isr::iterateTable(
    lsst::afw::image::MaskedImage<ImageT, MaskT> &maskedImage,
    std::vector<vectorType> &lookupTable
    ) {

    std::vector<vectorType>::iterator tableIter = lookupTable.begin();
    std::vector<vectorType>::iterator tableIterEnd = lookupTable.end();
    const int miCols = static_cast<int>(maskedImage.getCols());
    const int miRows = static_cast<int>(maskedImage.getRows());

    lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> miRowAcc(maskedImage);
                
    for (int miRow = 0; miRow < miRows; miRow++, miRowAcc.nextRow()) {
        lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> miColAcc = miRowAcc;
        for (int miCol = 0; miCol < miCols; miCol++, miColAcc.nextCol()) {
            if (tableIter != tableIterEnd) {
                    
                double tableVal = lookupTable[tableIter];     
                // Note: if adding a constant with no associated uncertainty,
                // the variance does not change
                *chunkColAcc.image += static_cast<ImageT>(tableVal);
                *tableIter++;
            }
        }
    }
}
