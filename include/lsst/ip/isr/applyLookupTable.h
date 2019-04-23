// -*- LSST-C++ -*- 

/*
 * LSST Data Management System
 * Copyright 2016 LSST Corporation.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

#ifndef LSST_IP_ISR_APPLY_LOOKUP_TABLE
#define LSST_IP_ISR_APPLY_LOOKUP_TABLE

#include "ndarray.h"

#include "lsst/afw/image.h"

namespace lsst {
namespace ip {
namespace isr {

    /**
    Add the values in a lookup table to an image, e.g. for non-linearity correction

    The algorithm is as follows:
        numOutOfRange = 0
        For each i,j of the image:
            lookupInd = int(indOffset + image[i,j])
            if lookupInd not in range [0, table.size() - 1]:
                set lookupInd to nearest edge and increment numOutOfRange
            image[i,j] += table[lookupInd]
        return numOutOfRange

    @param[in,out] image  image to which to add the values; modified in place
    @param[in] table  lookup table
    @param[in] indOffset  scalar added to image value before truncating to lookup column

    @return the number of pixels whose values were out of range
    */
    template<typename PixelT>
    int applyLookupTable(
        afw::image::Image<PixelT> &image,
        ndarray::Array<PixelT, 1, 1> const &table,
        PixelT indOffset
    );

}}} // lsst::ip::isr

#endif
