// -*- lsst-c++ -*-

/*
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
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

#include <cmath>

#include "lsst/geom.h"
#include "lsst/afw/math.h"
#include "lsst/afw/math/Statistics.h"
#include "lsst/ip/isr/isr.h"

namespace lsst { namespace ip { namespace isr {

template <typename PixelT>
size_t maskNans(afw::image::MaskedImage<PixelT> const& mi, afw::image::MaskPixel maskVal,
                afw::image::MaskPixel allow)
{
    typedef typename afw::image::MaskedImage<PixelT>::x_iterator x_iterator;
    size_t nPix = 0;
    for (int y = 0; y != mi.getHeight(); ++y) {
        for (x_iterator ptr = mi.row_begin(y), end = mi.row_end(y); ptr != end; ++ptr) {
            if (!(ptr.mask() & allow) && (!std::isfinite(ptr.image()) ||
                                          !std::isfinite(ptr.variance()))) {
                nPix += 1;
                ptr.mask() |= maskVal;
            }
        }
    }
    return nPix;
}

template<typename ImagePixelT>
std::vector<double> fitOverscanImage(
    afw::image::Image<ImagePixelT> const& overscan,
    bool isTransposed
) {
    typedef afw::image::Image<ImagePixelT> Image;

    /**
    This is transposed here to match the existing numpy-array ordering.
    This effectively transposes the image for us.
    **/
    const int height = overscan.getHeight();
    const int width  = overscan.getWidth();

    int length = height;
    if (isTransposed) {
        length = width;
    }

    std::vector<double> values(length);

    const int x0 = overscan.getX0();
    const int y0 = overscan.getY0();
    for (int x = 0; x < length; ++x) {
        /**
        geom::Box2I bbox       = geom::Box2I( geom::Point2I(0, y),
                                              geom::Point2I(0, width) );
        The above was how this was defined before ticket #1556.  As I understand it
        the following is the new way to do this
        **/
        geom::Box2I bbox;
        if (isTransposed) {
            bbox = geom::Box2I(geom::Point2I(x0 + x,y0), geom::Extent2I(1,height));
        }
        else {
            bbox = geom::Box2I(geom::Point2I(x0,y0 + x), geom::Extent2I(width,1));
        }
        Image mi         = Image(overscan, bbox);
        afw::math::Statistics stats = afw::math::makeStatistics(mi, afw::math::MEDIAN);

        values[x]    = stats.getValue(afw::math::MEDIAN);
    }
    return values;
}

std::string between(std::string &s, char ldelim, char rdelim) {
    std::string::iterator b(s.begin());
    std::string::iterator e(s.end());
    std::string::iterator lp;
    std::string::iterator rp;

    std::string result;

    if((lp = std::find(b, e, ldelim)) != e)
        if((rp = std::find(++lp, e, rdelim)) != e)
            result = std::string(lp, rp);

    return result;
}

// Explicit instantiations

template
std::vector<double> fitOverscanImage<int>(
    afw::image::Image<int> const&, bool isTransposed);
template
std::vector<double> fitOverscanImage<float>(
    afw::image::Image<float> const&, bool isTransposed);
template
std::vector<double> fitOverscanImage<double>(
    afw::image::Image<double> const&, bool isTransposed);

template class CountMaskedPixels<float>;
template class CountMaskedPixels<double>;
template class CountMaskedPixels<int>;

// Function to mask nans in a masked image
template size_t maskNans<float>(afw::image::MaskedImage<float> const&, afw::image::MaskPixel,
                                afw::image::MaskPixel);
template size_t maskNans<double>(afw::image::MaskedImage<double> const&, afw::image::MaskPixel,
                                 afw::image::MaskPixel);
template size_t maskNans<int>(afw::image::MaskedImage<int> const&, afw::image::MaskPixel,
                              afw::image::MaskPixel);

}}} // namespace lsst::ip::isr
