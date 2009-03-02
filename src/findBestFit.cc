
// -*- LSST-C++ -*- // fixed format comment for emacs
/**
  * \file
  *
  * \ingroup isr
  *
  * \brief Utility function to find the best fit function to an image region.
  *
  * \author Nicole M. Silvestri 
  *         University of Washington
  *         nms@astro.washington.edu
  *
  * \version
  *
  * LSST Legalese here...
  */

#include <string>
#include <vector>

#include "lsst/pex/exceptions.h"
#include "lsst/ip/isr/Isr.h"

/**\brief Iterate through a given MaskedImage to fill vectors for the pixel
  * values, variances, and col and row positions.  Pass these vectors with policy
  * information to a minimization routine to perform Chi^2.
  *
  * \return FitResults for best parameter fit.
  *
  * \throw lsst::pex::exceptions::RuntimeErrorException if equested function is not implemented.
  */ 

template<typename ImagePixelT>
lsst::afw::math::FitResults lsst::ip::isr::findBestFit(
    lsst::afw::image::MaskedImage<ImagePixelT> const &maskedImage,
    std::string const &funcForm,
    int funcOrder,
    double stepSize
) {
    typedef lsst::afw::image::MaskedImage<ImagePixelT> MaskedImage;

    // Find the best fit function to the MaskedImage region. Best fit determined via Chi^2 minimization.

    std::vector<double> parameterList;
    std::vector<double> stepSizeList;
    
    //std::fill(parameterList.begin(), parameterList.end(), numParams);
    //std::fill(stepSizeList.begin(), stepSizeList.end(), stepSize);

    // collapse the image region down to a vector.  Get vectors of measurements,
    // variances, and x,y positions

    std::vector<double> measurementList;
    std::vector<double> varianceList;
    std::vector<double> colPositionList;
    
    const int miHeight = maskedImage.getHeight();
    

    for (int y = 0; y < miHeight; ++y) {
        int miCol = 0;
        for (typename MaskedImage::x_iterator miPtr = maskedImage.row_begin(y), end = maskedImage.row_end(y); miPtr != end; ++miPtr, ++miCol) {
            measurementList.push_back(miPtr.image());
            varianceList.push_back(miPtr.variance());
            colPositionList.push_back(miCol);
        }
    }   
      
    // Ideally we need a function factory to which we feed different
    // functional forms and get out a Function1...but for now,
    // pick a few typical choices and code them here (poly for now
    // 10/22/08...ADD THE OTHERS LATER)
    
    if (funcForm == "POLYNOMIAL") {

        lsst::afw::math::PolynomialFunction1<double> polyFunc1(funcOrder);
        double sigma = 1.0; // initial guess

        const unsigned int nParams = polyFunc1.getNParameters();

        std::vector<double> parameterList = polyFunc1.getParameters();
        std::vector<double> stepSizeList(nParams);
        polyFunc1.setParameters(parameterList);
        for (unsigned int i = 0; i < nParams; ++i) {
            stepSizeList[i] = stepSize;
        }

        // find the best fit function
        lsst::afw::math::FitResults fitResults = lsst::afw::math::minimize(
            polyFunc1,
            parameterList,
            stepSizeList,
            measurementList, ///< pixel values
            varianceList,    ///< variance for each pixel value
            colPositionList, ///< x positions
            sigma
            );

        return fitResults;
        
    } else if (funcForm == "SPLINE") {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException, "Function 'SPLINE' Not Implemented Yet.");

    } else if (funcForm == "CHEVYCHEV") {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException, "Function 'CHEVYCHEV' Not Implemented Yet.");

    } else {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException, "Function Not Implemented. Use 'POLYNOMIAL', 'SPLINE', or 'CHEVYCHEV'.");
    }
}

/************************************************************************/
/* Explicit instantiations */

template
lsst::afw::math::FitResults lsst::ip::isr::findBestFit(
    lsst::afw::image::MaskedImage<float> const &maskedImage,
    std::string const &funcForm,
    int funcOrder,
    double stepSize
    );

template
lsst::afw::math::FitResults lsst::ip::isr::findBestFit(
    lsst::afw::image::MaskedImage<double> const &maskedImage,
    std::string const &funcForm,
    int funcOrder,
    double stepSize
    );
/************************************************************************/
