
// -*- LSST-C++ -*- // fixed format comment for emacs
/**
  * \file
  *
  * \ingroup isr
  *
  * \brief Utility function to find the best fit function to an image region.
  *
  * \author Nicole M. Silvestri, University of Washington
  *
  * Contact: nms@astro.washington.edu
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

// \brief
//
// \return FitResults for best parameter fit

typedef double vectorType;
typedef double funcType;

template<typename ImageT, typename MaskT> 
lsst::afw::math::FitResults lsst::ip::isr::findBestFit(
    lsst::afw::image::MaskedImage<ImageT, MaskT> const &maskedImage,
    std::string &funcForm,
    int &funcOrder
    ) {

    // Find the best fit function to the overscan region and apply this
    // function to the Chunk Exposure.  Best fit determined via Chi^2 minimization.

    std::vector<vectorType> parameterList;
    std::vector<vectorType> stepSizeList;
    std::fill(parameterList.begin(), parameterList.end(), numParams);
    std::fill(stepSizeList.begin(), stepSizeList.end(), stepSize);

    // collapse the overscan region down to a vector.  Get vectors of
    // measurements, variances, and x,y positions

    std::vector<vectorType> overscanMeasurementList;
    std::vector<vectorType> overscanVarianceList;
    std::vector<vectorType> colPositionList;
    std::vector<vectorType> rowPositionList;
    
    lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> miRowAcc(maskedImage);
       
    const int miCols = static_cast<int>(maskedImage.getCols());
    const int miRows = static_cast<int>(maskedImage.getRows());

    for (int miRow = 0; miRow < miRows; miRow++, miRowAcc.nextRow()) {
        lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> miColAcc = miRowAcc;
        for (int miCol = 0; miCol < miCols; miCol++, miColAcc.nextCol()) {
            overscanMeasurementList.push_back(*miColAcc.image);
            overscanVarianceList.push_back(*miColAcc.variance);
            colPositionList.push_back(miCol);
            rowPositionList.push_back(miRow);
        }
    }   
      
    // ideally we need a function factory to which we feed different
    // functional forms and get out a Function1...but for now,
    // pick a few typical choices and code them here (poly for now
    // 10/22/08...ADD THE OTHERS LATER)
    
    if (funcForm == "polynomial") {

        lsst::afw::math::PolynomialFunction1<funcType> polyFunc1(funcOrder);
        double sigma = 1.0; // initial guess

        // find the best fit function
        lsst::afw::math::FitResults overscanFit  = lsst::afw::math::minimize(
            polyFunc1,
            parameterList,
            stepSizeList,
            overscanMeasurementList, ///< overscan values
            overscanVarianceList,    ///< variance for each value
            colPositionList,         ///< x position
            sigma
            );

        return overscanFit;
    } else {
        throw lsst::pex::exceptions::Runtime(std::string("Function Not Implemented"));
    }
}

/************************************************************************/
/* Explicit instantiations */


template
void lsst::ip::isr::findBestFit(
    lsst::afw::image::MaskedImage<float, lsst::afw::image::maskPixelType> &maskedImage,
    std::string &funcForm,
    int &funcOrder
    );

template
void lsst::ip::isr::findBestFit(
    lsst::afw::image::MaskedImage<double, lsst::afw::image::maskPixelType> &maskedImage,
    std::string &funcForm,
    int &funcOrder
    );

