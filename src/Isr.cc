// -*- lsst-c++ -*-

#include <lsst/ip/isr/Isr.h>
#include <lsst/afw/math.h>
#include <lsst/afw/math/Statistics.h>

namespace isr   = lsst::ip::isr;
namespace image = lsst::afw::image;
namespace math  = lsst::afw::math;

// Functions
template<typename ImagePixelT, typename FunctionT>
void isr::fitOverscanImage(
    boost::shared_ptr<math::Function1<FunctionT> > &overscanFunction,
    image::MaskedImage<ImagePixelT> const& overscan,
    double ssize,
    int sigma
) {
    typedef image::MaskedImage<ImagePixelT> MaskedImage;


    const int height = overscan.getHeight();
    const int width  = overscan.getWidth();
    std::vector<double> values(height);
    std::vector<double> errors(height);
    std::vector<double> positions(height);

    std::vector<double> parameters(overscanFunction->getNParameters(), 0.);
    std::vector<double> stepsize(overscanFunction->getNParameters(), ssize);
    
    for (int y = 0; y < height; ++y) {
        image::BBox bbox       = image::BBox( image::PointI(0, y),
                                              image::PointI(0, width) );
        MaskedImage mi         = MaskedImage(overscan, bbox);
        math::Statistics stats = math::makeStatistics(*(mi.getImage()), math::MEAN | math::STDEV);

        values[y]    = stats.getValue(math::MEAN);
        errors[y]    = stats.getValue(math::STDEV);
        positions[y] = y;
     
    }
    lsst::afw::math::FitResults fitResults = math::minimize(
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

image::BBox isr::BBoxFromDatasec(std::string datasection) { 
    
    const char begin('[');
    const char end(']');
    const char delim1(':');
    const char delim2(',');
    
    std::string temp(between(datasection, delim2, end));
    std::size_t colonPos = temp.find(":");
    
    // NOTE: atoi() needs to be passed a c_str() to get the int out
    // ALSO: note the fits convention requires an adjustment by 1
    image::PointI startPt(
        atoi(between(datasection, begin, delim1).c_str()) - 1,
        atoi(between(datasection, delim2, delim1).c_str()) - 1);
    //std::cout << "start = " << startPt.getX() << ", " << startPt.getY() << std::endl;
    
    
    image::PointI endPt(
        atoi(between(datasection, delim1, delim2).c_str()) - 1,
        atoi(temp.substr(colonPos + 1).c_str()) - 1);
    //std::cout << "end = " << endPt.getX() << ", " << endPt.getY() << std::endl;

    image::BBox bbox = image::BBox(startPt, endPt);
    return bbox;
}



// Explicit instantiations

template
void isr::fitOverscanImage(
     boost::shared_ptr<math::Function1<double> > &overscanFunction, 
    image::MaskedImage<float> const& overscan,
    double ssize,
    int sigma);

template
void isr::fitOverscanImage(
     boost::shared_ptr<math::Function1<double> > &overscanFunction,
    image::MaskedImage<double> const& overscan,
    double ssize,
    int sigma);

// Integer classes make no sense for multiplicative table
//   unless you change the image type
template class isr::LookupTableMultiplicative<float>;
template class isr::LookupTableMultiplicative<double>;

// Only integer images make sense for a replacement table
template class isr::LookupTableReplace<int>;
