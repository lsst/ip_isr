// -*- LSST-C++ -*- // fixed format comment for emacs
/**
  * \file
  *
  * \ingroup isr
  *
  * \brief Implementation of the templated stage, Overscan Correct
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
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

#include "boost/cstdint.hpp"
#include "boost/format.hpp"
#include "boost/shared_ptr.hpp"
#include "vw/Math/BBox.h"

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

std::string between(std::string& s, char ldelim, char rdelim)
{
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


// std::string between(const std::string& s, char ldelim, char rdelim)
// {
//     static std::string::const_iterator b(s.begin());
//     static std::string::const_iterator e(s.end());
//     static std::string::const_iterator lp;
//     static std::string::const_iterator rp;
    
//     std::string result;
    
//     if((lp = std::find(b, e, ldelim)) != e)
//         if((rp = std::find(++lp, e, rdelim)) != e)
//             result = std::string(lp, rp);
    
//     return result;
// }

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

typedef double vectorType;
typedef double funcType;
std::string overStage = "lsst.ip.isr.overscanCorrectAndTrimChunkExposure";

template<typename ImageT, typename MaskT>
void lsst::ip::isr::overscanCorrectAndTrimChunkExposure(
    lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure, 
    lsst::pex::policy::Policy &isrPolicy, 
    lsst::pex::policy::Policy &datasetPolicy 
    ) {

    // Start with some setup and information gathering...
    // Get the Chunk MaskedImageand metadata from the Chunk Exposure

    lsst::pex::logging::TTrace<3>("Entering ISR stage: %s", overStage);

    lsst::afw::image::MaskedImage<ImageT, MaskT> chunkMaskedImage = chunkExposure.getMaskedImage();
    lsst::daf::base::DataProperty::PtrType chunkMetadata = chunkMaskedImage.getImage()->getMetaData();    

    // Make sure we haven't run this sub-stage previously

    lsst::daf::base::DataProperty::PtrType isrOverFlag = chunkMetadata->findUnique("OVER_END");  // overscan stage processing flag
    if (isrOverFlag) {
        lsst::pex::logging::TTrace<3>("ISR stage %s: Terminated.  Overscan correction has previously been performed on the Chunk Exposure.", overStage);
        throw lsst::pex::exceptions::Runtime(std::string("Overscan Correction has previously been applied to the Chunk Exposure - terminating stage."));
    }  
    
    // Parse the ISR policy file for the overscan and trim sub-stage information
   
    lsst::pex::policy::Policy::Ptr overscanPolicy;
    //std::string function;
    std::string constantMeth;
    std::string method;
    std::string funcForm;
    int funcOrder;
    double numParams;
    double stepSize;
    bool smooth;
    bool sigClip; 
    std::string datasecKey;
    std::string biassecKey;
    std::string trimsecKey;
    try {
        overscanPolicy = isrPolicy.getPolicy("overscanPolicy");
        method = overscanPolicy->getString("method");
        smooth = overscanPolicy->getBool("smooth");
        sigClip = overscanPolicy->getBool("sigClip");
        datasecKey = datasetPolicy.getString("datasec");
        biassecKey = datasetPolicy.getString("biassec");
        trimsecKey = datasetPolicy.getString("trimsec");
    } catch (std::exception e) {
        throw lsst::pex::exceptions::Runtime(std::string("Can not parse the main ISR policy file (for the overscan policy)."));    
    }

    try {
        if (method == "function"){
            funcForm = overscanPolicy->getString("funcForm");
            funcOrder = overscanPolicy->getInt("funcOrder");
            stepSize = overscanPolicy->getDouble("stepSize");
            numParams = overscanPolicy->getDouble("numParams");
        } else {
            constantMeth = overscanPolicy->getString("constantMeth");
        }
    } catch (std::exception e) {
        throw lsst::pex::exceptions::Runtime(std::string("Can not parse the overscan policy file for function information.")); 
    }

    int smoothSigma;
    int smoothOrder;
    std::string smoothOpt;
    try {
        if (smooth = true) {
            smoothOpt = overscanPolicy->getString("smoothOpt");
            if (smoothOpt == "gaussian") {
                smoothSigma = overscanPolicy->getInt("smoothSigma");
            } else {
                smoothOrder = overscanPolicy->getInt("smoothOrder");
            }
        }
    } catch (std::exception e) {
        throw lsst::pex::exceptions::Runtime(std::string("Can not parse the overscan policy file for smoothing information."));
    }
   
    int sigClipVal;
    try {
        if (sigClip = true) {
            sigClipVal = overscanPolicy->getInt("sigClipVal");
        }
    } catch (std::exception e) {
        throw lsst::pex::exceptions::Runtime(std::string("Can not parse the overscan policy file for sigma-clipping information."));
    }

    // Get the necessary metadata for the stage

    lsst::daf::base::DataProperty::PtrType datasecProp = chunkMetadata->findUnique(datasecKey);
    std::string datasec; 
    if (datasecProp) {
        datasec = boost::any_cast<std::string>(datasecProp->getValue());
        std::cout << "DATASEC: " << datasec << std::endl;
    } else {
        throw lsst::pex::exceptions::Runtime(std::string("Could not get DATASEC from the image metadata"));
    }
   
    lsst::daf::base::DataProperty::PtrType biassecProp = chunkMetadata->findUnique(biassecKey);
    std::string biassec;
    if (biassecProp) {
        biassec = boost::any_cast<std::string>(biassecProp->getValue());
        std::cout << "BIASSEC: " << biassec << std::endl;
    } else {
        throw lsst::pex::exceptions::Runtime(std::string("Could not get BIASSEC from the image metadata"));
    } 

    lsst::daf::base::DataProperty::PtrType trimsecProp = chunkMetadata->findUnique(trimsecKey);  
    std::string trimsec;
    if (trimsecProp) {
        trimsec = boost::any_cast<std::string>(trimsecProp->getValue());
        std::cout << "TRIMSEC: " << trimsec << std::endl;
    } else {
        throw lsst::pex::exceptions::Runtime(std::string("Could not get TRIMSEC from the image metadata"));
    } 

    // need to parse the biassec and trimsec into column and row numbers...
    // typical pattern for the data sections is: [####:####,####:####]
    
    const char begin('[');
    const char end(']');
    const char delim1(':');
    const char delim2(',');

   
    std::string temp1(between(biassec, delim2, end));
    std::size_t position1 = temp1.find(":");

    // note: atoi() needs to be passed a c_str() to get the int
    const int overscanColsStart = atoi(between(biassec, begin, delim1).c_str());
    std::cout << "OverscanColsStart: " << overscanColsStart << std::endl;
    const int overscanColsEnd = atoi(between(biassec, delim1, delim2).c_str());
    std::cout << "OverscanColsEnd: " << overscanColsEnd << std::endl;
    const int overscanRowsStart = atoi(between(biassec, delim2, delim1).c_str());
    std::cout << "OverscanRowsStart: " << overscanRowsStart << std::endl;
    const int overscanRowsEnd = atoi(temp1.substr(position1 + 1).c_str());
    std::cout << "OverscanRowsEnd: " << overscanRowsEnd << std::endl;
    const int overscanColSpan = overscanColsEnd - overscanColsStart;
    const int overscanRowSpan = overscanRowsEnd - overscanRowsStart;

    std::string temp2(between(trimsec, delim2, end));   
    std::size_t position2 = temp2.find(":");

    const int trimColsStart = atoi(between(trimsec, begin, delim1).c_str());
    std::cout << "TrimColsStart: " << trimColsStart << std::endl;
    const int trimColsEnd = atoi(between(trimsec, delim1, delim2).c_str());
    std::cout << "TrimColsEnd: " << trimColsEnd << std::endl;
    const int trimRowsStart = atoi(between(trimsec, delim2, delim1).c_str());
    std::cout << "TrimRowsStart: " << trimRowsStart << std::endl;
    const int trimRowsEnd = atoi(temp2.substr(position2 + 1).c_str()); 
    std::cout << "TrimRowsEnd: " << trimRowsEnd << std::endl;
    const int trimColSpan = trimColsEnd - trimColsStart;
    const int trimRowSpan = trimRowsEnd - trimRowsStart;

    
    // create a MaskedImage holding the overscan region to be subtracted and trimmed
    // QQQ: NEED TO TRIM RAMP UP HERE in addition to overscan region. How do I
    // determine appropriate ramp to trim?
    const vw::BBox2i overscanBbox = vw::BBox2i(overscanColsStart, 
                                               overscanRowsStart, 
                                               overscanColSpan, 
                                               overscanRowSpan);

    lsst::afw::image::Exposure<ImageT, MaskT> overscanExposure = chunkExposure.getSubExposure(overscanBbox);
    std::string overscanOut = "overscanStripExposure";
    overscanExposure.writeFits(overscanOut);
    lsst::afw::image::MaskedImage<ImageT, MaskT> overscanMaskedImage = overscanExposure.getMaskedImage();
    const int overscanCols = static_cast<int>(overscanMaskedImage.getCols());
    const int overscanRows = static_cast<int>(overscanMaskedImage.getRows());
    const int chunkCols = static_cast<int>(chunkMaskedImage.getCols());
    const int chunkRows = static_cast<int>(chunkMaskedImage.getRows());
   
    // trim the Chunk Exposure so we can now fit or apply the overscan
    // correction to only the pixels associated with the Chunk Exposures pixels

    const vw::BBox2i trimmedBBox = vw::BBox2i(trimColsStart,
                                              trimRowsStart,
                                              trimColSpan,
                                              trimRowSpan);

    lsst::afw::image::Exposure<ImageT, MaskT> trimmedChunkExposure = chunkExposure.getSubExposure(trimmedBBox);
    
    //get the new Chunk MaskedImage  
    lsst::afw::image::MaskedImage<ImageT, MaskT> trimmedChunkMaskedImage = trimmedChunkExposure.getMaskedImage();
    lsst::daf::base::DataProperty::PtrType trimmedChunkMetadata = trimmedChunkMaskedImage.getImage()->getMetaData();

    const int trimmedChunkCols = static_cast<int>(trimmedChunkMaskedImage.getCols());
    const int trimmedChunkRows = static_cast<int>(trimmedChunkMaskedImage.getRows());
    std::cout << "OrigCols: " << chunkCols<< std::endl;
    std::cout << "OrigRows: " << chunkRows<< std::endl;
    std::cout << "TrimmedCols_AmpA: " << trimmedChunkCols<< std::endl;
    std::cout << "TrimmedRows_AmpA: " << trimmedChunkRows<< std::endl;
    std::cout << "OverscanCols_AmpA: " << overscanCols<< std::endl;
    std::cout << "OverscanRows_AmpA: " << overscanRows<< std::endl;
 
    // Remove the datasec and biassec from the trimmed Chunk Exposure's
    // metadata...which currently lives in both the image and the variance
    // image.

    trimmedChunkMetadata->deleteAll(datasecKey);
    trimmedChunkMetadata->deleteAll(biassecKey);
    lsst::daf::base::DataProperty::PtrType trimmedVarianceMetadata = trimmedChunkMaskedImage.getVariance()->getMetaData();
    trimmedVarianceMetadata->deleteAll(datasecKey);
    trimmedVarianceMetadata->deleteAll(biassecKey);
    
    trimmedChunkMaskedImage.setMetadata(trimmedChunkMetadata);
    trimmedChunkMaskedImage.setMetadata(trimmedVarianceMetadata);

    if (method == "function"){

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

        if (funcForm == "polynomial") {

            lsst::afw::math::PolynomialFunction1<funcType> polyFunc1(funcOrder);
            //lsst::afw::math::Function1<funcType> function1(funcOrder);

            // find the best fit function
            lsst::afw::math::FitResults overscanFit  = lsst::afw::math::minimize( 
                polyFunc1, 
                parameterList, 
                stepSizeList, 
                overscanMeasurementList, ///< overscan values 
                overscanVarianceList,    ///< variance for each value 
                colPositionList,   ///< x position  
                //rowPositionList,   ///< y position (if Function2 function)
                sigma 
                );  

            std::vector<vectorType> parameters;
            for (unsigned int i = 0; i < overscanFit.parameterList.size(); ++i){
                parameters[i] = overscanFit.parameterList[i];      
            }
            unsigned int order = overscanFit.parameterList.size() - 1;
            lsst::afw::math::PolynomialFunction1<funcType> polyFunction(order);
            polyFunction.setParameters(parameters);

            lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> chunkRowAcc(trimmedChunkMaskedImage);
                
            for (int chunkRow = 0; chunkRow < trimmedChunkRows; chunkRow++, chunkRowAcc.nextRow()) {
                lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> chunkColAcc = chunkRowAcc;
                for (int chunkCol = 0; chunkCol < trimmedChunkCols; chunkCol++, chunkColAcc.nextCol()) {
                    *chunkColAcc.image = static_cast<ImageT>(polyFunction(*chunkColAcc.image));
                    *chunkColAcc.variance = static_cast<ImageT>(polyFunction(*chunkColAcc.image) * polyFunction(*chunkColAcc.image));
                }
            }
            // Record the provenance
            // QQQ: is there a way of getting the goodness-of-fit for the polynomial??

            lsst::daf::base::DataProperty::PtrType isrFunctionFlag(new lsst::daf::base::DataProperty("OVER_FN", funcForm));     
            trimmedChunkMetadata->addProperty(isrFunctionFlag);
            lsst::daf::base::DataProperty::PtrType isrOrderFlag(new lsst::daf::base::DataProperty("OVER_OR", funcOrder));        
            trimmedChunkMaskedImage.setMetadata(trimmedChunkMetadata);  

        } else if (funcForm == "spline") {
            // not yet implemented
            throw lsst::pex::exceptions::InvalidParameter(std::string("Spline is not yet implemented."));
        } else {
            throw lsst::pex::exceptions::InvalidParameter(std::string("Invalid functional form for overscan fit requested."));
        }
    } else {

        // Subtract a constant value.  For now, compute the mean for the
        // constant value to be subtracted. ADD THE OTHERS LATER...
        
        if (constantMeth == "mean") {

            // Compute the number of elements (n), mean (mu) and standard
            // deviation (sigma) borrowing some code from Russell Owen's
            // medianBinApprox

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
            // subtract off the overscan mean form the Chunk Masked Image
            trimmedChunkMaskedImage -= mu;

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
            lsst::daf::base::DataProperty::PtrType isrMeanFlag(new lsst::daf::base::DataProperty("OVER_MU", mu));
            trimmedChunkMetadata->addProperty(isrMeanFlag);
            lsst::daf::base::DataProperty::PtrType isrStdevFlag(new lsst::daf::base::DataProperty("OVER_SD", sigma));
            trimmedChunkMaskedImage.setMetadata(trimmedChunkMetadata);    

        } else if (constantMeth == "median"){
            // not yet implemented
            throw lsst::pex::exceptions::InvalidParameter(std::string("Median is not yet implemented."));
        } else if (constantMeth == "mode"){
            // not yet implemented
            throw lsst::pex::exceptions::InvalidParameter(std::string("Mode is not yet implemented."));
        } else {
            throw lsst::pex::exceptions::InvalidParameter(std::string("Invalid method for computing the overscan value requested."));
        }
    }
    // Record final sub-stage provenance to the Image Metadata
    lsst::daf::base::DataProperty::PtrType isrEndFlag(new lsst::daf::base::DataProperty("OVER_END", std::string("Completed Successfully"))); 
    trimmedChunkMetadata->addProperty(isrEndFlag);
    trimmedChunkMaskedImage.setMetadata(trimmedChunkMetadata);

    chunkMaskedImage = trimmedChunkMaskedImage;

    // Calculate additional SDQA Metrics here??

    //Issue a logging message if the sub-stage executes without issue to this point!
    lsst::pex::logging::TTrace<3>("ISR stage, %s, completed successfully.", overStage);
    lsst::pex::logging::TTrace<3>("Leaving ISR stage: %s", overStage);
}


	
/************************************************************************/
// /* Explicit instantiations */

template 
void lsst::ip::isr::overscanCorrectAndTrimChunkExposure(
    lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> &chunkExposure,
    lsst::pex::policy::Policy &isrPolicy,
    lsst::pex::policy::Policy &datasetPolicy
    );

template 
void lsst::ip::isr::overscanCorrectAndTrimChunkExposure(
    lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> &chunkExposure,
    lsst::pex::policy::Policy &isrPolicy,
    lsst::pex::policy::Policy &datasetPolicy
    );

/************************************************************************/


