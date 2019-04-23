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

#include <cstdint>

#include "lsst/pex/exceptions.h"
#include "lsst/ip/isr/applyLookupTable.h"

namespace lsst {
namespace ip {
namespace isr {

template<typename PixelT>
int applyLookupTable(
    afw::image::Image<PixelT> &image,
    ndarray::Array<PixelT, 1, 1> const &table,
    PixelT indOffset
) {
    if (table.size() == 0u) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthError,
            "Lookup table has zero size."
        );
    }
    int numOutOfRange = 0;
    int const maxLookupCol = table.size() - 1;
    for (int col = 0, imHeight = image.getHeight(); col < imHeight; ++col) {
        for (auto imPtr = image.row_begin(col), end = image.row_end(col); imPtr != end; ++imPtr) {
            int lookupCol = indOffset + *imPtr;
            if (lookupCol < 0) {
                lookupCol = 0;
                ++numOutOfRange;
            } else if (lookupCol > maxLookupCol) {
                lookupCol = maxLookupCol;
                ++numOutOfRange;
            }
            *imPtr += table[lookupCol];
        }
    }
    return numOutOfRange;
}

#define INSTANTIATE(T) \
    template int applyLookupTable<T>(afw::image::Image<T> &, ndarray::Array<T, 1, 1> const &, T);

INSTANTIATE(float);
INSTANTIATE(double);

}}} // lsst::ip::isr
