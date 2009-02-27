
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

#include "lsst/pex/exceptions.h"
#include "lsst/ip/isr/isr.h"

/**
* \brief Given the maskedImage and the lookupTable, 
*
* \return the maskedImage adjusted to the table values
*
* \throw lsst::pex::exceptions::RangeErrorException if lookup table index out of range
*/
template<typename ImagePixelT>
void lsst::ip::isr::iterateTable(
    lsst::afw::image::MaskedImage<ImagePixelT> &maskedImage,
    std::vector<double> const &lookupTable
) {
    typedef lsst::afw::image::MaskedImage<ImagePixelT> MaskedImage;

    int maxInd = lookupTable.size() - 1;

    const int miHeight = maskedImage.getHeight();

    // Set the pixels row by row, to avoid repeated checks for end-of-row
    for (int y = 0; y < miHeight; ++y) {
        for (typename MaskedImage::x_iterator miPtr = maskedImage.row_begin(y), end = maskedImage.row_end(y); miPtr != end; ++miPtr) {
            int ind = static_cast<int>(miPtr.image());
            if ((ind < 0) || (ind > maxInd)) {
                 throw LSST_EXCEPT(lsst::pex::exceptions::RangeErrorException, "Lookup table index out of range.");
            } else {
                miPtr.image() += lookupTable[ind]; 
            } 
        }
    }
}

/************************************************************************/
/* Explicit instantiations */

template
void lsst::ip::isr::iterateTable(
    lsst::afw::image::MaskedImage<float> &maskedImage,
    std::vector<double> const &lookupTable
    );

template
void lsst::ip::isr::iterateTable(
    lsst::afw::image::MaskedImage<double> &maskedImage,
    std::vector<double> const &lookupTable
    );

/************************************************************************/
