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
#include <vector>

#include "vw/Math/BBox.h"

#include <lsst/daf/base.h> 
#include <lsst/daf/data/LsstBase.h>  	
#include <lsst/afw/image.h>
#include <lsst/afw/image/MaskedImage.h>
#include <lsst/afw/math.h>

/** \brief Remove all non-astronomical source counts from the Chunk Exposure's
  * pixels.
  * 
  */
	
namespace lsst {
namespace ip {
namespace isr {
    
    typedef double vectorType;
    typedef double funcType;    

    template<typename ImageT, typename MaskT>
    double easyMean(
        lsst::afw::image::MaskedImage<ImageT, MaskT> &maskedImage
        );

    template<typename ImageT, typename MaskT>
    void fitFunctionToImage(
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &maskedImage,
        lsst::afw::math::PolynomialFunction1<funcType> const &polyFunction
        );

    template<typename ImageT, typename MaskT> 
    void iterateTable(
        lsst::afw::image::MaskedImage<ImageT, MaskT> &maskedImage,
        std::vector<vectorType> &lookupTable
        );

    template<typename ImageT, typename MaskT> 
    lsst::afw::math::FitResults findBestFit(
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &maskedImage,
        std::string const &funcForm,
        int funcOrder,
        double stepSize
        );

    vw::BBox2i stringParse(
        std::string &section
        );

}}} // namespace lsst::ip::isr
	
#endif // !defined(LSST_IP_ISR_ISR_H)
