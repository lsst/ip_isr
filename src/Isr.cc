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
 

#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/math.h"
#include "lsst/afw/math/Statistics.h"
#include "lsst/ip/isr.h"

namespace lsst { namespace ip { namespace isr {

// Functions
template <typename ImageT>
void LookupTableReplace<ImageT>::apply(afw::image::MaskedImage<ImageT> &image, float gain) const {
    double igain = 1.0 / gain;
    int nPixTooHigh = 0;
    int nPixTooLow = 0;
    for (int y = 0; y != image.getHeight(); ++y) {
        for (x_iterator ptr = image.row_begin(y), end = image.row_end(y); ptr != end; ++ptr) {
            int ind = static_cast<int>(ptr.image() + 0.5);  // Rounded pixel value
            if (ind < 0) {
                ind = 0;
                ++nPixTooLow;
            } else if (ind >= _max) {
                ind = _max - 1;
                ++nPixTooHigh;
            }
            PixelT p = PixelT(_table[ind],  (*ptr).mask(), _table[ind] * igain);
            *ptr = p;
        }
    }
    if ((nPixTooHigh > 0) || (nPixTooLow > 0)) {
        // log message
        pex::logging::TTrace<1>("lsst.ip.isr.LookupTableReplace.apply", 
            "Data truncated; %d pixels were < 0; %d pixels were >= %d", nPixTooLow, nPixTooHigh, _max);
    }
}

template <typename PixelT>
size_t maskNans(afw::image::MaskedImage<PixelT> const& mi, afw::image::MaskPixel maskVal,
                afw::image::MaskPixel allow)
{
    typedef typename afw::image::MaskedImage<PixelT>::x_iterator x_iterator;
    size_t nPix = 0;
    for (int y = 0; y != mi.getHeight(); ++y) {
        for (x_iterator ptr = mi.row_begin(y), end = mi.row_end(y); ptr != end; ++ptr) {
            if (!(ptr.mask() & allow) && (!utils::lsst_isfinite(ptr.image()) ||
                                          !utils::lsst_isfinite(ptr.variance()))) {
                nPix += 1;
                ptr.mask() |= maskVal;
            }
        }
    }
    return nPix;
}

template<typename ImagePixelT, typename FunctionT>
void fitOverscanImage(
    boost::shared_ptr< afw::math::Function1<FunctionT> > &overscanFunction,
    afw::image::MaskedImage<ImagePixelT> const& overscan,
    double ssize,
    int sigma
) {
    typedef afw::image::MaskedImage<ImagePixelT> MaskedImage;


    const int height = overscan.getHeight();
    const int width  = overscan.getWidth();
    std::vector<double> values(height);
    std::vector<double> errors(height);
    std::vector<double> positions(height);

    std::vector<double> parameters(overscanFunction->getNParameters(), 0.);
    std::vector<double> stepsize(overscanFunction->getNParameters(), ssize);
    
    for (int y = 0; y < height; ++y) {
        /**
        afw::geom::Box2I bbox       = afw::geom::Box2I( afw::geom::Point2I(0, y),
                                              afw::geom::Point2I(0, width) );
        The above was how this was defined before ticket #1556.  As I understand it
        the following is the new way to do this
        **/
        afw::geom::Box2I bbox = afw::geom::Box2I(afw::geom::Point2I(0,y), afw::geom::Point2I(width,y));
        MaskedImage mi         = MaskedImage(overscan, bbox, afw::image::PARENT);
        afw::math::Statistics stats = afw::math::makeStatistics(*(mi.getImage()), afw::math::MEAN | afw::math::STDEV);

        values[y]    = stats.getValue(afw::math::MEAN);
        errors[y]    = stats.getValue(afw::math::STDEV);
        positions[y] = y;
     
    }
    afw::math::FitResults fitResults = afw::math::minimize(
        *overscanFunction,
        parameters,
        stepsize,
        values,
        errors,
        positions,
        sigma
        );
    
    overscanFunction->setParameters(fitResults.parameterList);
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
void fitOverscanImage(
     boost::shared_ptr<afw::math::Function1<double> > &overscanFunction, 
    afw::image::MaskedImage<float> const& overscan,
    double ssize,
    int sigma);

template
void fitOverscanImage(
     boost::shared_ptr<afw::math::Function1<double> > &overscanFunction,
    afw::image::MaskedImage<double> const& overscan,
    double ssize,
    int sigma);

template class CountMaskedPixels<float>;
template class CountMaskedPixels<double>;

// Integer classes make no sense for multiplicative table
//   unless you change the image type
template class LookupTableMultiplicative<float>;
template class LookupTableMultiplicative<double>;

// Only integer images make sense for a replacement table
template class LookupTableReplace<int>;
// But we turn our images into floats immediately, so use it
template class LookupTableReplace<float>;
// Function to mask nans in a masked image
template size_t maskNans<float>(afw::image::MaskedImage<float> const&, afw::image::MaskPixel,
                                afw::image::MaskPixel);
template size_t maskNans<double>(afw::image::MaskedImage<double> const&, afw::image::MaskPixel,
                                 afw::image::MaskPixel); 

}}} // namespace lsst::ip::isr
