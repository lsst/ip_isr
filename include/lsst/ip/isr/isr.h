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

/*
 * Implementation of the templated Instrument Signature Removal
 * stage of the nightly LSST Image Processing Pipeline.
 */

#ifndef LSST_IP_ISR_ISR_H
#define LSST_IP_ISR_ISR_H

#include <memory>
#include <string>
#include <vector>
#include <cmath>

#include "lsst/afw/math.h"
#include "lsst/afw/math/Statistics.h"
#include "lsst/afw/image.h"
#include "lsst/pex/exceptions/Exception.h"

namespace lsst {
namespace ip {
namespace isr {

    /**
     * Remove all non-astronomical counts from the Chunk Exposure's pixels.
     */
    template <typename ImageT, typename MaskT=lsst::afw::image::MaskPixel>
    class CountMaskedPixels {
    public:
        typedef typename lsst::afw::image::MaskedImage<ImageT>::x_iterator x_iterator;
        CountMaskedPixels() :
            _count(0) {} ;
        virtual ~CountMaskedPixels() {};

        // Clear the accumulator
        void reset() { _count = 0; }

        // Count pixels
        void apply(lsst::afw::image::MaskedImage<ImageT> const& image,
                   MaskT bitmask) {
            reset();
            for (int y = 0; y != image.getHeight(); ++y) {
                for (x_iterator ptr = image.row_begin(y); ptr != image.row_end(y); ++ptr) {
                    if ( ((*ptr).mask() & bitmask) == bitmask ) {
                        _count += 1;
                    }
                }
            }
        }

        // Return the total counts
        int getCount() const { return _count; }

    private:
        int _count;
    };

    /// Mask NANs in an image
    ///
    /// NANs in the image or variance that are not already masked by
    /// the 'allow' value are masked with the 'maskVal'.
    ///
    /// @return Number of pixels masked
    template <typename PixelT>
    size_t maskNans(
        afw::image::MaskedImage<PixelT> const& mi, ///< Input image
        afw::image::MaskPixel maskVal,  ///< Bit mask value to give a NaN
        afw::image::MaskPixel allow=0 ///< Retain NANs with this bit mask (0 to mask all NANs)
        );


    template<typename ImagePixelT>
    std::vector<double> fitOverscanImage(
        lsst::afw::image::MaskedImage<ImagePixelT> const& overscan,
        std::vector<std::string> badPixelMask,
        bool isTransposed
        );

}}} // namespace lsst::ip::isr

#endif // !defined(LSST_IP_ISR_ISR_H)
