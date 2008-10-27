// -*- LSST-C++ -*- // fixed format comment for emacs
/**
  * \file
  *
  * \ingroup isr
  *
  * \brief Implementation of the templated subStage, Overscan Correct
  * and Trim Chunk Exposure, of the Instrument Signature Removal stage of 
  * the LSST Image Processing Pipeline.
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
#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>

#include <boost/cstdint.hpp>
#include <boost/format.hpp>

#include <vw/Math/BBox.h>

#
#include <lsst/afw/image/Exposure.h>
#include <lsst/afw/image/Mask.h>
#include <lsst/afw/image/MaskedImage.h>
#include <lsst/afw/image/PixelAccessors.h>
#include <lsst/afw/math/Function.h>
#include <lsst/afw/math/FunctionLibrary.h>
#include <lsst/afw/math/minimize.h>
#include <lsst/daf/base/DataProperty.h>
#include <lsst/pex/exceptions/Exception.h>
#include <lsst/pex/logging/Trace.h>
#include <lsst/pex/policy/Policy.h>

#include "lsst/ip/isr/isr.h"

// Simple string iterator to help parse data sections.
// Finds the string between two delimiting characters
// Adapted from C-code snippet by Mike Wahler.

std::string between(const std::string& s, char ldelim, char rdelim)
{
    static std::string::const_iterator b(s.begin());
    static std::string::const_iterator e(s.end());
    static std::string::const_iterator lp;
    static std::string::const_iterator rp;
    
    std::string result;
    
    if((lp = std::find(b, e, ldelim)) != e)
        if((rp = std::find(++lp, e, rdelim)) != e)
            result = std::string(lp, rp);
    
    return result;
}

/** \brief Subtract a single constant (using all input pixels for the statistic)
  * or one-dimensional function from the Chunk Exposure which varies along the
  * overscan.  Use spline or polynomial fit along the overscan region (want a
  * slowly varying function).  For the 1-D representation, Pan-STARRS determines
  * the input values to the fit as representations on the coordinate along the
  * overscan, with the statistic derived from the pixels in the perpendicular
  * direction at each location.  Sigma-clipping on input data is a necessary
  * option.  
  *
  * The Image is also trimmed to remove the overscan region and any 
  * unexposed pixels.
  *
  * \return chunkExposure corrected for overscan and trimmed of all
  * non-illuminated pixels
  *
  * \throw Runtime if sage has already been run on the image
  * \throw Runtime if any Policy or metadata value can not be obtained
  * \throw InvalidParameter if functional form for the overscan fit is invalid 
  * 
  * TO DO (as of Wed 10/22/08):
  * - Implement sigma-clipping
  * - Trim ramp-up on overscan (do not include in fits)
  * - Implement additional constant modes and spline fit
  * - add a smoothing option to the functional fit? 
  * - add any additional SDQA statistics requested by SDQA team
  */

template<typename ImageT, typename MaskT>
lsst::afw::image::Exposure<ImageT, MaskT> overscanCorrectAndTrimChunkExposure(
    lsst::afw::image::Exposure<ImageT, MaskT> const &chunkExposure,
    lsst::pex::policy::Policy &isrPolicy,
    lsst::pex::policy::Policy &datasetPolicy
    ) {

    // Get the Chunk MaskedImage from the Chunk Exposure
    lsst::afw::image::MaskedImage<ImageT, MaskT> chunkMaskedImage = chunkExposure.getMaskedImage();
    // Get the Image Metadata from the Chunk Science Image
    lsst::daf::base::DataProperty::PtrType chunkMetadata = chunkMaskedImage.getImage->getMetadata();
    
    // Make sure we haven't run this sub-stage previously
    lsst::daf::base::DataProperty::PtrType isrOverscan = chunkMetadata->findUnique("ISR_OVERSCAN");  // overscan processing flag
    int isrOver;
    if (isrOverscan) {
        isrOver = isrOverscan->getValue();
    } else {
    throw lsst::pex::exceptions::Runtime(std::string("In ") + __func__ + std::string(": Could not get sub-stage processing flag from the image metadata"));	
    }
    
    if (isrOver != 0) {
        throw lsst::pex::exceptions::Runtime(std::string("Overscan Correction has already been applied to this Chunk Exposure"));
    }

    // Parse the ISR policy file for the overscan and trim sub-stage information
    lsst::pex::policy::Policy overscanPolicy = isrPolicy.getPolicy("overscanPolicy");
    std::string function;
    std::string constantMeth;
    std::string method = overscanPolicy.getString("method");
    std::string funcForm;
    const int funcOrder;
    double numParams;
    double stepSize;
    if (method = "function"){
        funcForm = overscanPolicy.getString("funcForm");
        funcOrder = overscanPolicy.getInt("funcOrder");
        stepSize = overscanPolicy.getDouble("stepSize");
        numParams = overscanPolicy.getDouble("numParams");
    } else {
        constantMeth = overscanPolicy.getString("constantMeth");
    }
    bool smooth = overscanPolicy.getBool("smooth");
    const int smoothSigma;
    const int smoothOrder;
    std::string smoothOpt;
    if (smooth = true) {
        smoothOpt = overscanPolicy.getString("smoothOpt");
        if (smoothOpt = "gaussian") {
            smoothSigma = overscanPolicy.getInt("smoothSigma");
        } else {
            smoothOrder = overscanPolicy.getInt("smoothOrder");
        }
    }
    bool sigClip = overscanPolicy.getBool("sigClip");
    const int sigClipVal;
    if (sigClip = true) {
        sigClipVal = overscanPolicy.getInt("sigClipVal");
    }

    // Get the necessary metadata for the sub-stage
    lsst::daf::base::DataProperty::PtrType ccdTypeProp = chunkMetadata->findUnique("CCDTYPE");  // ccd image type (bias, dark, etc.)
    std::string ccdType;
    if (ccdTypeProp) {
        ccdType = ccdTypeProp->getValue();
    } else {
        throw lsst::pex::exceptions::Runtime(std::string("In ") + __func__ + std::string(": Could not get CCDTYPE from the image metadata"));	
    }
    lsst::daf::base::DataProperty::PtrType datasecProp = chunkMetadata->findUnique("DATASEC");
    std::string datasec = datasecProp->getValue();
    lsst::daf::base::DataProperty::PtrType biassecProp = chunkMetadata->findUnique("BIASSEC");
    std::string biassec = biassecProp->getValue();
    lsst::daf::base::DataProperty::PtrType trimsecProp = chunkMetadata->findUnique("TRIMSEC");  
    std::string trimsec = trimsecProp->getValue();

    // need to parse the biassec and trimsec into column and row numbers...
    // typical pattern for the data sections is: [####:####,####:####]
    
    const char begin('[');
    const char end(']');
    const char delim1(':');
    const char delim2(',');

    std::string temp1(between(biassec, delim2, end));

    // note: atoi() needs to be passed a c_str() to get the int
    const int overscanColsStart = atoi(between(biassec, begin, delim1).c_str());
    const int overscanColsEnd = atoi(between(biassec, delim1, delim2).c_str());
    const int overscanRowsStart = atoi(between(biassec, delim2, delim1).c_str());
    const int overscanRowsEnd = atoi(between(temp1, delim1, temp1.end()).c_str());

    std::string temp2(between(trimsec, delim2, end));

    const int trimColsStart = atoi(between(trimsec, begin, delim1).c_str());
    const int trimColsEnd = atoi(between(trimsec, delim1, delim2).c_str());
    const int trimRowsStart = atoi(between(trimsec, delim2, delim1).c_str());
    const int trimRowsEnd= atoi(between(temp2, delim1, temp2.end()).c_str());

    // create a MaskedImage holding overscan region to be subtracted and trimmed
    // NEED TO TRIM RAMP UP HERE...
    const vw::BBox2i overscanBbox = vw::BBox2i(overscanColsStart, 
                                               overscanRowsStart, 
                                               overscanColsEnd, 
                                               overscanRowsEnd);

    lsst::afw::image::Exposure<ImageT, MaskT> overscanExposure = chunkExposure.getSubExposure(overscanBbox);
    lsst::afw::image::MaskedImage<ImageT, MaskT> overscanMaskedImage = overscanExposure.getMaskedImage();
    const int overscanCols = static_cast<int>(overscanExposure.getCols());
    const int overscanRows = static_cast<int>(overscanExposure.getRows());
    const int numCols = static_cast<int>(chunkExposure.getCols());
    const int numRows = static_cast<int>(chunkExposure.getRows());

    // trim the Chunk Exposure so we can now fit or apply the overscan
    // correction to only the pixels associated with the Chunk Exposures pixels

    const vw::BBox2i trimmedBBox = vw::BBox2i(trimColsStart,
                                              trimRowsStart,
                                              trimColsEnd,
                                              trimRowsEnd);

    lsst::afw::image::Exposure<ImageT, MaskT> trimmedChunkExposure = chunkExposure.getSubExposure(trimmedBBox);
    // replace the old datasec with the new one which is actually the trimsec
    chunkMetadata->deleteAll("DATASEC", false);
    chunkMetadata->addProperty(lsst::daf::base::DataProperty("DATASEC", trimsec));
    //get the new Chunk MaskedImage  
    lsst::afw::image::MaskedImage<ImageT, MaskT> trimmedChunkMaskedImage = trimmedChunkExposure.getMaskedImage();
    lsst::daf::base::DataProperty::PtrType trimmedChunkMetadata = trimmedChunkMaskedImage.getImage->getMatedata();
 
    const int chunkCols = trimmedChunkMaskedImage.getCols();
    const int chunkRows = trimmedChunkMaskedImage.getRows();

    if (method = "function"){

        // Find the best fit function to the overscan region and apply this
        // function to the Chunk Exposure.  Best fit determined via Chi^2 minimization.

        std::vector<double> parameterList;
        std::vector<double> stepSizeList; 
        std::fill(parameterList.begin(), parameterList.end(), numParams);
        std::fill(stepSizeList.begin(), stepSizeList.end(), stepSize); 
        
        // collapse the overscan region down to a vector.  Get vectors of
        // measurements, variances, and x,y positions
 
        std::vector<double> overscanMeasurementList;
        std::vector<double> overscanVarianceList;
        std::vector<double> colPositionList;
        std::vector<double> rowPositionList;
        double sigma = 1.0; // initial guess

        lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> overscanRowAcc(overscanMaskedImage);
       
        for (int chunkRow = 0; chunkRow < overscanRows; chunkRow++, overscanRowAcc.nextRow()) {
            lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> overscanColAcc = overscanRowAcc;
            for (int chunkCol = 0; chunkCol < overscanCols; chunkCol++, overscanColAcc.nextCol()) {
                overscanMeasurementList.push_back(*overscanColAcc.image);
                overscanVarianceList.push_back(*overscanColAcc.variance);
                colPositionList.push_back(chunkCol);
                rowPositionList.push_back(chunkRow);
            } // for column loop
        } // for row loop     

        // ideally we need a function factory to which we feed different
        // functional forms and get out a Function1 or Function2...but for now,
        // pick a few typical choices and code them here (poly for now
        // 10/22/08...ADD THE OTHERS LATER)

        if (funcForm = "polynomial") {

            //lsst::afw::math::PolynomialFunction1<double> function(funcOrder);
            lsst::afw::math::Function1<double> function;

            // find the best fit function
            lsst::afw::math::FitResults overscanFit  = lsst::afw::math::minimize( 
                function, 
                parameterList, 
                stepSizeList, 
                overscanMeasurementList, ///< overscan values 
                overscanVarianceList,    ///< variance for each value 
                colPositionList,   ///< x position  
                //rowPositionList,   ///< y position (if Function2 function)
                sigma 
                );  

            std::vector<double> parameters;
            for (unsigned int i = 0; i < overscanFit.parameterList.size(); ++i){
                parameters[i] = overscanFit.parameterList[i];      
            }
            unsigned int order = overscanFit.parameterList.size() - 1;
            lsst::afw::math::PolynomialFunction1<double> polyFunction(order);
            polyFunction.setParameters(parameters);

            lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> chunkRowAcc(trimmedChunkMaskedImage);
                
            for (int chunkRow = 0; chunkRow < numRows; chunkRow++, chunkRowAcc.nextRow()) {
                lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> chunkColAcc = chunkRowAcc;
                for (int chunkCol = 0; chunkCol < numCols; chunkCol++, chunkColAcc.nextCol()) {
                    *chunkColAcc.image = static_cast<ImageT>(polyfunc(*chunkColAcc.image));
                    *chunkColAcc.varaince = static_cast<ImageT>(polyfunc(*chunkColAcc.image) * polyfunc(*chunkColAcc.image));
                }
            }
            // Record the provenance
            trimmedChunkMetadata->addProperty(lsst::daf::base::DataProperty("ISR_OVERSCAN_FUNCTION", funcForm));
            trimmedChunkMetadata->addProperty(lsst::daf::base::DataProperty("ISR_OVERSCAN_ORDER", funcOrder));
            //trimmedChunkMetadata->addProperty(lsst::daf::base::DataProperty("ISR_OVERSCAN_STDEV", overscanFit.sigma));
            trimmedChunkMaskedImage.setMetadata(trimmedChunkMetadata);  

        } else if (funcForm = "spline") {
            // not yet implemented
        } else {
            throw lsst::pex::exceptions::InvalidParameter(std::string("Invalid functional form for overscan fit."));
        }
    } else {

        // Subtract a constant value.  For now, compute the mean for the
        // constant value to be subtracted. ADD THE OTHERS LATER...
        
        if (constantMeth = "mean") {

            // Compute the number of elements (n), mean (mu) and standard
            // deviation (sigma) borrowing some code from Russell Owen's
            // medianBinApprox in coadd/kaiser to do so

            lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> overscanRowAcc(overscanMaskedImage);
            long int n = 0;
            double sum = 0;
       
            for (int chunkRow = 0; chunkRow < overscanRows; chunkRow++, overscanRowAcc.nextRow()) {
                lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> overscanColAcc = overscanRowAcc;
                for (int chunkCol = 0; chunkCol < overscanCols; chunkCol++, overscanColAcc.nextCol()) {
                    n++;
                    sum += static_cast<double>(*overscanColAcc.image);                 
                } 
            }     
            // the mean
            double mu = sum/static_cast<double>(n);
            // subtract off the overscan mean form the Chunk Exposure
            trimmedChunkExposure -= mu;

            double sumSq = 0;
            for (int chunkRow = 0; chunkRow < overscanRows; chunkRow++, overscanRowAcc.nextRow()) {
                lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> overscanColAcc = overscanRowAcc;
                for (int chunkCol = 0; chunkCol < overscanCols; chunkCol++, overscanColAcc.nextCol()) {
                    n++;
                    double val = static_cast<double>(*overscanColAcc.image) - mu; 
                    sumSq += val * val; 
                } 
             }
            // the standard deviation
            double sigma = std::sqrt(sumSq/static_cast<double>(n));

            // Record the provenance
            trimmedChunkMetadata->addProperty(lsst::daf::base::DataProperty("ISR_OVERSCAN_MEAN", mu));
            trimmedChunkMetadata->addProperty(lsst::daf::base::DataProperty("ISR_OVERSCAN_STDEV", sigma));
            trimmedChunkMaskedImage.setMetadata(trimmedChunkMetadata);    

        } else if (constantMeth = "median"){
            // not yet implemented
        } else if (constantMeth = "mode"){
            // not yet implemented
        } else {
            throw lsst::pex::exceptions::InvalidParameter(std::string("Invalid method for computing the overscan value."));
        }
    }
    // Record final sub-stage provenance to the Image Metadata
    trimmedChunkMetadata->addProperty(lsst::daf::base::DataProperty("ISR_OVERSCAN", "Complete"));
    trimmedChunkMaskedImage.setMetadata(trimmedChunkMetadata);

    chunkExposure = trimmedChunkExposure;

    // Calculate additional SDQA Metrics here??

    //Issue a logging message if the sub-stage executes without issue to this point!
    lsst::pex::logging::TTrace<7>(std::string("ISR sub-stage") +__func__ + std::string("completed successfully."));

}


	
/************************************************************************/
// /* Explicit instantiations */

// template 
// lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> overscanCorrectAndTrimChunkExposure(
//     lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> const &chunkExposure,
//     lsst::pex::policy::Policy &isrPolicy,
//     lsst::pex::policy::Policy &datasetPolicy
//     );

// template 
// lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> overscanCorrectAndTrimChunkExposure(
//     lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> const &chunkExposure,
//     lsst::pex::policy::Policy &isrPolicy,
//     lsst::pex::policy::Policy &datasetPolicy
//     );

/************************************************************************/


