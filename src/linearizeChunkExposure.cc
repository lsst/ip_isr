// -*- LSST-C++ -*- // fixed format comment for emacs
/**
  * \file
  *
  * \ingroup isr
  *
  * \brief Implementation of the templated subStage, Linearize Chunk 
  *  Exposure, of the Instrument Signature Removal stage forthe nightly LSST
  *  Image Processing Pipeline.
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

#include "boost/cstdint.hpp"
#include "boost/format.hpp"
#include "boost/shared_ptr.hpp"

#include <lsst/afw/image/Exposure.h>
#include <lsst/afw/image/Mask.h>
#include <lsst/afw/image/MaskedImage.h>
#include <lsst/afw/math.h>
#include <lsst/afw/image/PixelAccessors.h>
#include <lsst/afw/math/Function.h>
#include <lsst/afw/math/FunctionLibrary.h>
#include <lsst/daf/base/DataProperty.h>
#include <lsst/pex/exceptions/Exception.h>
#include <lsst/pex/logging/Trace.h>
#include <lsst/pex/policy/Policy.h>

#include "lsst/ip/isr/isr.h"

typedef double vectorType;
typedef double funcType;

/** \brief Correct for non-linear effects of mapping from electrons to ADU
 * (amplifier roll-off). Must be done AFTER bias subtraction.  Apply correction
 * as a function of pixel value from either a lookup table for fitted function
 * (polynomial or spline).
 *
 * \return chunkExposure properly linearized
 *
 * \throw Runtime if this sub-stage has been run previously
 * \throw NotFound if any Policy or metadata value can not be obtained
 * \throw InvalidParameter if functional form for the lineaization fit is invalid 
 */
template<typename ImageT, typename MaskT>
void lsst::ip::isr::linearizeChunkExposure(
    lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure,    
    lsst::pex::policy::Policy &isrPolicy, 
    lsst::pex::policy::Policy &datasetPolicy,
    std::vector<vectorType> &linearizeLookupTable
    ) { 

    // Get the Chunk MaskedImage from the Chunk Exposure
    lsst::afw::image::MaskedImage<ImageT, MaskT> chunkMaskedImage = chunkExposure.getMaskedImage();
    lsst::daf::base::DataProperty::PtrType chunkMetadata = chunkMaskedImage.getImage()->getMetadata();

    std::string subStage = "Linearize Chunk Exposure";
    
    // Check that this ISR sub-stage has not previously been run on the Exposure
    lsst::daf::base::DataProperty::PtrType isrLinearizeField = chunkMetadata->findUnique("ISR_LINCOR");
    if (isrLinearizeField) {   
        lsst::pex::logging::TTrace<3>("In %s: Exposure has already been corrected.  Terminating ISR sub-stage for this Chunk Exposure.", subStage);
        throw lsst::pex::exceptions::Runtime(std::string("Linearization previously performed"));
    }

    int numCols = chunkMaskedImage.getCols();
    int numRows = chunkMaskedImage.getRows();

// PARSE THE POLICY FILES

    // Get the type of linearization to be applied (function or lookup table)
    lsst::pex::policy::Policy::Ptr linearizePolicy = isrPolicy.getPolicy("linearizePolicy");
    std::string linearizeType = linearizePolicy->getString("linearizeType");

    // currently accepts as "FUNCTION" either a polynomial or spline
    if (linearizeType == "FUNCTION"){
        // Get the functional form and coefficients for the polynomial from the policy
        std::string funcForm = linearizePolicy->getString("funcForm");
        const int funcOrder = linearizePolicy->getInt("funcOrder");

        if (funcForm == "POLYNOMIAL"){                    
            std::vector<vectorType> parameters(funcOrder + 1);
            lsst::afw::math::PolynomialFunction1<funcType> polyFunction(funcOrder);
            for (unsigned int j = 0; j < parameters.size(); ++j) {
                parameters[j] = static_cast<vectorType>(1 + funcOrder - j);
            }
            polyFunction.setParameters(parameters);
            lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> chunkRowAcc(chunkMaskedImage);

            for (int chunkRow = 0; chunkRow < numRows; chunkRow++,chunkRowAcc.nextRow()) {
                lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> chunkColAcc = chunkRowAcc;
                for (int chunkCol = 0; chunkCol < numCols; chunkCol++, chunkColAcc.nextCol()) {
                    *chunkColAcc.image = static_cast<ImageT>(polyFunction(*chunkColAcc.image));
                    *chunkColAcc.variance = static_cast<ImageT>(polyFunction(*chunkColAcc.image) * polyFunction(*chunkColAcc.image));
                }
            }
        } else if (funcForm == "SPLINE"){
            // need to add a spline function to afw/math/FunctionLibrary
            // ptrFcn = FunctionPtr(new SplineFunction1(funcOrder));
        }
    }

    if (linearizeType == "LOOKUP"){

        // Assume lookup table is a vector of deltas (difference between the
        // image value and actual value)

        std::vector<vectorType>::iterator tableIter = linearizeLookupTable.begin();
        std::vector<vectorType>::iterator tableIterEnd = linearizeLookupTable.end();
        lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> chunkRowAcc(chunkMaskedImage);

        for (int chunkRow = 0; chunkRow < numRows; chunkRow++,chunkRowAcc.nextRow()) {
            lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> chunkColAcc = chunkRowAcc;
            for (int chunkCol = 0; chunkCol < numCols; chunkCol++, chunkColAcc.nextCol()) {
                if (*tableIter != *tableIterEnd) {
                    
                    double tableVal = linearizeLookupTable[*tableIter];     

                    *chunkColAcc.image += static_cast<ImageT>(tableVal);
                    *chunkColAcc.variance += static_cast<ImageT>(tableVal * tableVal);
 
                    *tableIter++;
                }
            }     
        }
    }

    //Record the sub-stage provenance to the Image Metadata
    chunkMetadata->addProperty(lsst::daf::base::DataProperty("ISR_LINCOR"));
    lsst::daf::base::DataProperty::PtrType linCorProp = chunkMetadata->findUnique("ISR_LINCOR");
    std::string exitTrue = "Completed Successfully";
    linCorProp->setValue(boost::any_cast<std::string>(exitTrue));
    chunkMaskedImage.setMetadata(chunkMetadata);

    //Calculate additional SDQA Metrics here. 


    //Issue a logging message if the sub-stage executes without issue
    lsst::pex::logging::TTrace<7>("ISR sub-stage, %s, completed successfully.", subStage);
         
}

/************************************************************************/
/* Explicit instantiations */

template
void lsst::ip::isr::linearizeChunkExposure(
    lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> &chunkExposure,    
    lsst::pex::policy::Policy &isrPolicy, 
    lsst::pex::policy::Policy &datasetPolicy,
    std::vector<vectorType> &linearizeLookupTable
    );

template
void lsst::ip::isr::linearizeChunkExposure(
    lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> &chunkExposure,    
    lsst::pex::policy::Policy &isrPolicy, 
    lsst::pex::policy::Policy &datasetPolicy,
    std::vector<vectorType> &linearizeLookupTable
    );

/************************************************************************/
