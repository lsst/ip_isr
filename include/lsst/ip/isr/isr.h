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

#include "lsst/afw/image.h"
#include "lsst/afw/math.h"

/** \brief Remove all non-astronomical source counts from the Chunk Exposure's
  * pixels.
  * 
  */
	
namespace lsst {
namespace ip {
namespace isr {
    
    template<typename ImagePixelT>
    void fitFunctionToImage(
        lsst::afw::image::MaskedImage<ImagePixelT> &maskedImage,
        lsst::afw::math::Function1<double> const &function
        );

    template<typename ImagePixelT>
    void iterateTable(
        lsst::afw::image::MaskedImage<ImagePixelT> &maskedImage,
        std::vector<double> const &lookupTable
        );


    template<typename ImagePixelT>
    lsst::afw::math::FitResults findBestFit(
        lsst::afw::image::MaskedImage<ImagePixelT> const &maskedImage,
        std::string const &funcForm,
        int funcOrder,
        double stepSize
        );

    lsst::afw::image::BBox stringParse(
        std::string &section
        );

}}} // namespace lsst::ip::isr
	
#endif // !defined(LSST_IP_ISR_ISR_H)
