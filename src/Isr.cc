// -*- lsst-c++ -*-

#include <lsst/pex/logging/Trace.h>
#include <lsst/afw/math.h>
#include <lsst/afw/math/Statistics.h>
#include <lsst/ip/isr/Isr.h>

namespace pexLog = lsst::pex::logging;
namespace afwImage = lsst::afw::image;
namespace afwMath = lsst::afw::math;
namespace ipIsr = lsst::ip::isr;

// Functions
template <typename ImageT>
void ipIsr::LookupTableReplace<ImageT>::apply(afwImage::MaskedImage<ImageT> &image, float gain) const {
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
        pexLog::TTrace<1>("lsst.ip.isr.LookupTableReplace.apply", 
            "Data truncated; %d pixels were < 0; %d pixels were >= %d", nPixTooLow, nPixTooHigh, _max);
    }
}

template<typename ImagePixelT, typename FunctionT>
void ipIsr::fitOverscanImage(
    boost::shared_ptr<afwMath::Function1<FunctionT> > &overscanFunction,
    afwImage::MaskedImage<ImagePixelT> const& overscan,
    double ssize,
    int sigma
) {
    typedef afwImage::MaskedImage<ImagePixelT> MaskedImage;


    const int height = overscan.getHeight();
    const int width  = overscan.getWidth();
    std::vector<double> values(height);
    std::vector<double> errors(height);
    std::vector<double> positions(height);

    std::vector<double> parameters(overscanFunction->getNParameters(), 0.);
    std::vector<double> stepsize(overscanFunction->getNParameters(), ssize);
    
    for (int y = 0; y < height; ++y) {
        afwImage::BBox bbox       = afwImage::BBox( afwImage::PointI(0, y),
                                              afwImage::PointI(0, width) );
        MaskedImage mi         = MaskedImage(overscan, bbox);
        afwMath::Statistics stats = afwMath::makeStatistics(*(mi.getImage()), afwMath::MEAN | afwMath::STDEV);

        values[y]    = stats.getValue(afwMath::MEAN);
        errors[y]    = stats.getValue(afwMath::STDEV);
        positions[y] = y;
     
    }
    afwMath::FitResults fitResults = afwMath::minimize(
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

afwImage::BBox ipIsr::BBoxFromDatasec(std::string datasection) { 
    
    const char begin('[');
    const char end(']');
    const char delim1(':');
    const char delim2(',');
    
    std::string temp(between(datasection, delim2, end));
    std::size_t colonPos = temp.find(":");
    
    // NOTE: atoi() needs to be passed a c_str() to get the int out
    // ALSO: note the fits convention requires an adjustment by 1
    afwImage::PointI startPt(
        atoi(between(datasection, begin, delim1).c_str()) - 1,
        atoi(between(datasection, delim2, delim1).c_str()) - 1);
    //std::cout << "start = " << startPt.getX() << ", " << startPt.getY() << std::endl;
    
    
    afwImage::PointI endPt(
        atoi(between(datasection, delim1, delim2).c_str()) - 1,
        atoi(temp.substr(colonPos + 1).c_str()) - 1);
    //std::cout << "end = " << endPt.getX() << ", " << endPt.getY() << std::endl;

    afwImage::BBox bbox = afwImage::BBox(startPt, endPt);
    return bbox;
}



// Explicit instantiations

template
void ipIsr::fitOverscanImage(
     boost::shared_ptr<afwMath::Function1<double> > &overscanFunction, 
    afwImage::MaskedImage<float> const& overscan,
    double ssize,
    int sigma);

template
void ipIsr::fitOverscanImage(
     boost::shared_ptr<afwMath::Function1<double> > &overscanFunction,
    afwImage::MaskedImage<double> const& overscan,
    double ssize,
    int sigma);

template class ipIsr::CountMaskedPixels<float>;
template class ipIsr::CountMaskedPixels<double>;

// Integer classes make no sense for multiplicative table
//   unless you change the image type
template class ipIsr::LookupTableMultiplicative<float>;
template class ipIsr::LookupTableMultiplicative<double>;

// Only integer images make sense for a replacement table
template class ipIsr::LookupTableReplace<int>;
// But we turn our images into floats immediately, so use it
template class ipIsr::LookupTableReplace<float>;
