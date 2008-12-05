
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

#include <lsst/afw/image.h>
#include <lsst/afw/math.h>
#include "lsst/ip/isr/isr.h"

/**\brief Iterate through a given MaskedImage to fill vectors for the pixel
  * values, variances, and col and row positions.  Pass these vectors with policy
  * information to a minimization routine to perform Chi^2.
  *
  * \return FitResults for best parameter fit.
  *
  * \throw Runtime if equested function is not implemented.
  */ 

typedef double vectorType;
typedef double funcType;

template<typename ImageT, typename MaskT> 
lsst::afw::math::FitResults lsst::ip::isr::findBestFit(
    lsst::afw::image::MaskedImage<ImageT, MaskT> const &maskedImage,
    std::string const &funcForm,
    int funcOrder,
    double stepSize
    ) {

    // Find the best fit function to the MaskedImage region.  Best fit
    // determined via Chi^2 minimization.

    std::vector<vectorType> parameterList;
    std::vector<vectorType> stepSizeList;
    
    //std::fill(parameterList.begin(), parameterList.end(), numParams);
    //std::fill(stepSizeList.begin(), stepSizeList.end(), stepSize);

    // collapse the image region down to a vector.  Get vectors of measurements,
    // variances, and x,y positions

    std::vector<vectorType> measurementList;
    std::vector<vectorType> varianceList;
    std::vector<vectorType> colPositionList;
    std::vector<vectorType> rowPositionList;
    
    lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> miRowAcc(maskedImage);
       
    const int miCols = static_cast<int>(maskedImage.getCols());
    const int miRows = static_cast<int>(maskedImage.getRows());

    for (int miRow = 0; miRow < miRows; miRow++, miRowAcc.nextRow()) {
        lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> miColAcc = miRowAcc;
        for (int miCol = 0; miCol < miCols; miCol++, miColAcc.nextCol()) {
            measurementList.push_back(*miColAcc.image);
            varianceList.push_back(*miColAcc.variance);
            colPositionList.push_back(miCol);
            rowPositionList.push_back(miRow);
        }
    }   
      
    // Ideally we need a function factory to which we feed different
    // functional forms and get out a Function1...but for now,
    // pick a few typical choices and code them here (poly for now
    // 10/22/08...ADD THE OTHERS LATER)
    
    if (funcForm == "POLYNOMIAL") {

        lsst::afw::math::PolynomialFunction1<funcType> polyFunc1(funcOrder);
        double sigma = 1.0; // initial guess

	const unsigned int nParams = polyFunc1.getNParameters();

	std::vector<vectorType> parameterList = polyFunc1.getParameters();
        std::vector<vectorType> stepSizeList(nParams);
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
        throw lsst::pex::exceptions::Runtime(std::string("Function 'SPLINE' Not Implemented Yet."));

    } else if (funcForm == "CHEVYCHEV") {
        throw lsst::pex::exceptions::Runtime(std::string("Function 'CHEVYCHEV' Not Implemented Yet."));

    } else {
        throw lsst::pex::exceptions::Runtime(std::string("Function Not Implemented. Use 'POLYNOMIAL', 'SPLINE', or 'CHEVYCHEV'."));
    }
}

/************************************************************************/
/* Explicit instantiations */

template
lsst::afw::math::FitResults lsst::ip::isr::findBestFit(
    lsst::afw::image::MaskedImage<float, lsst::afw::image::maskPixelType> const &maskedImage,
    std::string const &funcForm,
    int funcOrder,
    double stepSize
    );

template
lsst::afw::math::FitResults lsst::ip::isr::findBestFit(
    lsst::afw::image::MaskedImage<double, lsst::afw::image::maskPixelType> const &maskedImage,
    std::string const &funcForm,
    int funcOrder,
    double stepSize
    );
/************************************************************************/
