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
	
#include <vector>

#include <lsst/daf/base.h>
#include <lsst/daf/data/LsstBase.h>	
#include <lsst/afw/image/Exposure.h>
#include <lsst/pex/policy/Policy.h>

/** \brief Remove all non-astronomical source counts from the Chunk Exposure's
  * pixels.
  * 
  * The sequence of sub-stages within the Instrument Signature Removal (ISR)
  * stage of the nightly IPP are as follows:
  * 
  * Saturation Correction for Chunk Exposure
  * Overscan Correct and Trim Chunk Exposure
  * Bias Correct Chunk Exposure
  * Dark Current Correct Chunk Exposure
  * Linearize Chunk Exposure
  * Flat Field Correct Chunk Exposure
  *  (DC3 Stretch) - Pupil Image Correction
  *  - Illumination Correction (scattered light correction)
  * Defringe Chunk Exposure
  * (DC3 Stretch) Geometric Distortion Correction
  * (DC3 Stretch) Mask and Correct Additional Artifacts
  * (DC3 Stretch) Additional Flat Correction
  * (DC3 Stretch) Crosstalk Correct Chunk Exposure
  * (DC3 Stretch) Cosmic Ray Detection
  * 
  * Crosstalk Correction will be incorporated as an individual sub-stage for the
  * Data Release ISR stage.  It is not currently part of the nightly ISR stage
  * as we receive a crosstalk corrected image from the camera for the nightly
  * pipeline.
  * 
  */
	
namespace lsst {
namespace ip {
namespace isr {
	       
    typedef boost::uint16_t maskPixelType;
    
    template<typename ImageT, typename MaskT>
    lsst::afw::image::Exposure<ImageT, MaskT> saturationCorrectionForChunkExposure(
        lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure,
        lsst::pex::policy::Policy &isrPolicy,
	lsst::pex::policy::Policy &datasetPolicy
//        std::vector<float> &saturationLookupTable
        );
    
    template<typename ImageT, typename MaskT>
    lsst::afw::image::Exposure<ImageT, MaskT> overscanCorrectAndTrimChunkExposure(
        lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure,
        lsst::pex::policy::Policy &isrPolicy,
	lsst::pex::policy::Policy & datasetPolicy
//	lsst::afw::image::Exposure<ImageT, MaskT> &chunkOverscanExposure
        );

    template<typename ImageT, typename MaskT>
    lsst::afw::image::Exposure<ImageT, MaskT> biasCorrectChunkExposure(
        lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure,
        lsst::afw::image::Exposure<ImageT, MaskT> &masterChunkExposure,
	lsst::pex::policy::Policy &isrPolicy,
	lsst::pex::policy::Policy &datasetPolicy
        );

    template<typename ImageT, typename MaskT>
    lsst::afw::image::Exposure<ImageT, MaskT> darkCurrentCorrectChunkExposure(
        lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure,
        lsst::afw::image::Exposure<ImageT, MaskT> &masterChunkExposure,
	lsst::pex::policy::Policy &isrPolicy,
	lsst::pex::policy::Policy &datasetPolicy
        );

    template<typename ImageT, typename MaskT>
    lsst::afw::image::Exposure<ImageT, MaskT> linearizeChunkExposure(
        lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure,
        lsst::pex::policy::Policy &isrPolicy,
	lsst::pex::policy::Policy &datasetPolicy,
	std::vector<float> &linearizeLookupTable
        );

    template<typename ImageT, typename MaskT>
    lsst::afw::image::Exposure<ImageT, MaskT> flatFieldCorrectChunkExposure(
        lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure,
        lsst::afw::image::Exposure<ImageT, MaskT> &masterChunkExposure,
        lsst::pex::policy::Policy &isrPolicy,
	lsst::pex::policy::Policy &datasetPolicy
        );

    template<typename ImageT, typename MaskT>
    lsst::afw::image::Exposure<ImageT, MaskT> illuminationCorrection(
	lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure,
	lsst::afw::image::Exposure<ImageT, MaskT> &masterChunkExposure,
	lsst::pex::policy::Policy &isrPolicy,
	lsst::pex::policy::Policy &datasetPolicy
        );

    template<typename ImageT, typename MaskT>
    lsst::afw::image::Exposure<ImageT, MaskT> defringeChunkExposure(
        lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure,
        lsst::afw::image::Exposure<ImageT, MaskT> &masterChunkExposure,
        lsst::pex::policy::Policy &isrPolicy,
	lsst::pex::policy::Policy &datasetPolicy  
        );

// DC3 STRETCH GOALS

//    template<typename ImageT, typename MaskT>
//    lsst::afw::image::Exposure<ImageT, MaskT> geometricDistortionCorrection(
//	lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure,
//	lsst::pex::policy::Policy &isrPolicy
//      lsst::pex::policy::Policy &datasetPolicy
//      );

//    template<typename ImageT, typename MaskT>
//    lsst::afw::image::Exposure<ImageT, MaskT> pupilImageCorrection(
//	lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure,
//	lsst::afw::image::Exposure<ImageT, MaskT> &masterChunkExposure,
//	lsst::pex::policy::Policy &isrPolicy,
//	lsst::pex::policy::Policy &pupilPolicy        
//      );

//    template<typename ImageT, typename MaskT>
//    lsst::afw::image::Exposure<ImageT, MaskT> maskAndCorrectAdditionalArtifacts(
//	lsst::afw::image::Exposure(ImageT, MaskT) &chunkExposure,
//	lsst::pex::policy::Policy &isrPolicy
//      lsst::pex::policy::Policy &datasetPolicy
//      );

//    template<typename ImageT, typename MaskT>
//    lsst::afw::image::Exposure<ImageT, MaskT> crosstalkCorrectChunkExposure(
//      lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure,
//	lsst::pex::policy::Policy &isrPolicy,
//  	lsst::pex::policy::Policy &crosstalkPolicy,
//	std::vector &crosstalkLookupTable
//      );

//    template<typename ImageT, typename MaskT>
//    lsst::afw::image::Exposure<ImageT, MaskT> cosmicRayDetection(
//      lsst::afw::image::Exposure<ImageT, MaskT> &chunkExposure
//	lsst::pex::policy::Policy &isrPolicy,
//	lsst::pex::policy::Policy &cosmicRayPolicy
//      );

// LSST GOAL - not for DC3.  Requires tunable laser dome flate
//    template<typename ImageT, typename MaskT>
//    lsst::afw::image::Exposure<ImageT, MaskT> additionalFlatCorrection(         
//      );

}}} // namespace lsst::ip::isr
	
#endif // !defined(LSST_IP_ISR_ISR_H)
