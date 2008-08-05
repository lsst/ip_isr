// -*- LSST-C++ -*- // fixed format comment for emacs
/**
  * \file
  *
  * \ingroup isr
  *
  * \brief Implementation of the templated Instrument Signature Removal 
  * stage of the nightly LSST Image Processing Pipeline.
  *
  * \author Nicole M. Silvestri, University of Washington
  *
  * Contact: nms@astro.washington.edu
  *
  * \version
  *
  * LSST Legalese here...
  */
	
#ifndef LSST_IP_ISR_ISR_H
#define LSST_IP_ISR_ISR_H
	
#include <string>
	
#include <boost/shared_ptr.hpp>

#include <lsst/daf/base.h>
#include <lsst/daf/data/LsstBase.h>	
#include <lsst/afw/image/Exposure.h>
#include <lsst/afw/image/Image.h>
#include <lsst/afw/image/Mask.h>
#include <lsst/afw/math/Function.h>
#include <lsst/daf/base/DataProperty.h>
#include <lsst/pex/policy/Policy.h>

/** \brief Remove all non-astronomical source counts from the Chunk Exposure's
  * pixels.
  * 
  * The sequence of sub-stages within the Instrument Signature Removal (ISR)
  * stage of the nightly IPP are as follows:
  * 
  * (00) Assemble Chunk Exposure
  * (01) Saturation Correction for Chunk Exposure
  * (02) Overscan Correct Chunk Exposure
  * (03) Trim Chunk Exposure
  * (04) Bias Correct Chunk Exposure
  * (05) Dark Current Correct Chunk Exposure
  * (06) Linearize Chunk Exposure
  * (07) Flat Correct Chunk Exposure
  * (7a) Pupil Image Correction
  * (7b) Geometric Distortion Correction
  * (7c) Scattered Light Correction
  * (7d) Illumination Correction
  * (08) Mask Additional Artifacts
  * (09) Defringe Chunk Exposure
  * (--) Additional Flat Correction
  * (--) Crosstalk Correct Chunk Exposure
  * (--) Interpolate Over Masked Pixels - Utility Function called by (01) and (08)
  * 
  * Crosstalk Correction will be incorporated as an individual sub-stage for the
  * Archive Center ISR stage.  It is not currently part of the nightly ISR stage
  * as we receive a crosstalk corrected image from the camera for the nightly
  * pipeline.
  * 
  */
	
namespace lsst {
namespace ip {
namespace isr {
	       
    typedef boost::uint16_t maskPixelType;
    
    template<typename ImageT, typename MaskT>
    lsst::afw::image::Exposure<ImageT, MaskT> assembleChunkExposure(
        lsst::afw::image::Image<ImageT> &chunkImage,
        lsst::daf::base::DataProperty::PtrType &chunkMetaData,
        lsst::afw::image::Mask<MaskT> &badPixelMask       
        );
    
    template<typename ImageT, typename MaskT>
    lsst::afw::image::Exposure<ImageT, MaskT> saturationCorrectionForChunkExposure(
        lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure,
        lsst::daf::base::DataProperty::PtrType &chunkMetaData,
        lsst::pex::policy::Policy &policy
        //lsst::daf::data::DataProperty::PtrType &saturationLookUpTable
        );
    
    template<typename ImageT, typename MaskT>
    lsst::afw::image::Exposure<ImageT, MaskT> overscanCorrectChunkExposure(
        lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure,
        lsst::daf::base::DataProperty::PtrType &chunkMetaData,
        lsst::pex::policy::Policy &policy
        );

    template<typename ImageT, typename MaskT>
    lsst::afw::image::Exposure<ImageT, MaskT> trimChunkExposure(
        lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure,
        lsst::daf::base::DataProperty::PtrType &chunkMetaData,
        lsst::pex::policy::Policy &policy
        );

    template<typename ImageT, typename MaskT>
    lsst::afw::image::Exposure<ImageT, MaskT> biasCorrectChunkExposure(
        lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure,
        lsst::afw::image::Exposure<ImageT, MaskT> &masterExposure,
        lsst::daf::base::DataProperty::PtrType &chunkMetaData,
        lsst::daf::base::DataProperty::PtrType &masterMetaData
        );

    template<typename ImageT, typename MaskT>
    lsst::afw::image::Exposure<ImageT, MaskT> darkCurrentCorrectChunkExposure(
        lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure,
        lsst::afw::image::Exposure<ImageT, MaskT> &masterExposure,
        lsst::daf::base::DataProperty::PtrType &chunkMetaData,
        lsst::daf::base::DataProperty::PtrType &masterMetaData
        );

    template<typename ImageT, typename MaskT>
    lsst::afw::image::Exposure<ImageT, MaskT> linearizeChunkExposure(
        lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure,
        lsst::daf::base::DataProperty::PtrType &chunkMetaData,
  //     lsst::afw::math::Function::Function2<ReturnT> const &function,
        lsst::pex::policy::Policy &policy
        );

    template<typename ImageT, typename MaskT>
    lsst::afw::image::Exposure<ImageT, MaskT> flatCorrectChunkExposure(
        lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure,
        lsst::afw::image::Exposure<ImageT, MaskT> &masterExposure,
        lsst::daf::base::DataProperty::PtrType &chunkMetaData,
        lsst::daf::base::DataProperty::PtrType &masterMetaData,
        lsst::pex::policy::Policy &policy
        );

    template<typename ImageT, typename MaskT>
    lsst::afw::image::Exposure<ImageT, MaskT> scatteredLightCorrection(
        
        
        );

    template<typename ImageT, typename MaskT>
    lsst::afw::image::Exposure<ImageT, MaskT> geometricDistortionCorrection(
        
        
        );

    template<typename ImageT, typename MaskT>
    lsst::afw::image::Exposure<ImageT, MaskT> pupilImageCorrection(
        
        
        );

    template<typename ImageT, typename MaskT>
    lsst::afw::image::Exposure<ImageT, MaskT> illuminationCorrection(


        );

    template<typename ImageT, typename MaskT>
    lsst::afw::image::Exposure<ImageT, MaskT> maskAdditionalArtifacts(
        
        
        );

    template<typename ImageT, typename MaskT>
    lsst::afw::image::Exposure<ImageT, MaskT> defringeChunkExposure(
        lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure,
        lsst::afw::image::Exposure<ImageT, MaskT> &masterExposure,
        lsst::daf::base::DataProperty::PtrType &chunkMetaData,
        lsst::daf::base::DataProperty::PtrType &masterMetaData,
        lsst::pex::policy::Policy &policy  
        );

// Probable Utility Function that will need to be moved to afw/image/math

    template<typename ImageT, typename MaskT>
    lsst::afw::image::Exposure<ImageT, MaskT> interpolateOverMaskedPixels(
        lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure,
        lsst::daf::base::DataProperty::PtrType &chunkMetaData,
        lsst::pex::policy::Policy &policy  
        );


// // The following are not necessary for the nightly ISR stage - but needed
// // for the Archive Center ISR stage.

//     template<typename ImageT, typename MaskT>
//     lsst::afw::image::Exposure<ImageT, MaskT> crosstalkCorrectChunkExposure(
//         );

// // Will need additional correction to images based on dataCube model from
// // dome flats and aux telescope + models frm satellites + etc.

//     template<typename ImageT, typename MaskT>
//     lsst::afw::image::Exposure<ImageT, MaskT> additionalFlatCorrection(
//         );

}}} // namespace lsst::ip::isr
	
#endif // !defined(LSST_IP_ISR_ISR_H)
