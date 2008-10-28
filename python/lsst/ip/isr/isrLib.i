// -*- lsst-c++ -*-
%define isrLib_DOCSTRING
"
Python bindings for lsst::ip::isr code
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.ip.isr",docstring=isrLib_DOCSTRING) isrLib

// Suppress swig complaints; see afw/image/imageLib.i for more 
#pragma SWIG nowarn=362                 // operator=  ignored 


// Everything we will need in the _wrap.cc file
%{
#include <lsst/ip/isr/isr.h>
#include "boost/cstdint.hpp" 
#include "lsst/afw/image.h" 
#include "lsst/afw/math.h" 
%}

%init %{
%}

namespace boost {
    class bad_any_cast; // remove warning: Nothing known about 'boost::bad_any_cast'
}

// Everything whose bindings we will have to know about
%import "lsst/daf/data/LsstBase.h"  // avoid warning: Nothing known about base class 'lsst::daf::data::LsstBase' 
%import "lsst/afw/image/Mask.h" // needed so SWIG knows lsst::afw::image::maskPixelType = boost::uint16_t 
%include "lsst/p_lsstSwig.i"    // this needs to go first otherwise i do not know about e.g. boost
%include "lsst/afw/image/lsstImageTypes.i"  // vw and Image/Mask types and typedefs
%include "lsst/detection/detectionLib.i"    // need for Footprints

// handle C++ arguments that should be outputs in python
%apply int& OUTPUT { int& };
%apply float& OUTPUT { float& };
%apply double& OUTPUT { double& };

%pythoncode %{
def version(HeadURL = r"$HeadURL: svn+ssh://svn.lsstcorp.org/DMS/afw/trunk/python/lsst/ip/isr/isrLib.i $"):

    """Return a version given a HeadURL string; default: ip_isr's version"""
    return guessSvnVersion(HeadURL)

%}

// Here you give the names of the functions you want to Swig

%include "lsst/ip/isr/saturationCorrectionForChunkExposure.h"
%include "lsst/ip/isr/isr.h"
%include "lsst/ip/isr/interpolateOverMaskedPixels.h"

// %template(instrumentSignatureRemovalController)
//     lsst::ip::isr::instrumentSignatureRemovalController<float, lsst::afw::image::maskPixelType>;
// %template(instrumentSignatureRemovalController)
//     lsst::ip::isr::instrumentSignatureRemovalController<double, lsst::afw::image::maskPixelType>;

%template(saturationCorrectionForChunkExposure)
    lsst::ip::isr::saturationCorrectionForChunkExposure<float, lsst::afw::image::maskPixelType>;
%template(saturationCorrectionForChunkExposure)
    lsst::ip::isr::saturationCorrectionForChunkExposure<double, lsst::afw::image::maskPixelType>;

%template(overscanCorrectAndTrimChunkExposure)
    lsst::ip::isr::overscanCorrectAndTrimChunkExposure<float, lsst::afw::image::maskPixelType>;
%template(overscanCorrectAndTrimChunkExposure)
    lsst::ip::isr::overscanCorrectAndTrimChunkExposure<double, lsst::afw::image::maskPixelType>;

%template(overscanCorrectAndTrimChunkExposure)
    lsst::ip::isr::overscanCorrectAndTrimChunkExposure<float, lsst::afw::image::maskPixelType>;
%template(overscanCorrectAndTrimChunkExposure)
    lsst::ip::isr::overscanCorrectAndTrimChunkExposure<double, lsst::afw::image::maskPixelType>;

%template(biasCorrectChunkExposure)
    lsst::ip::isr::biasCorrectChunkExposure<float, lsst::afw::image::maskPixelType>;
%template(biasCorrectChunkExposure)
    lsst::ip::isr::biasCorrectChunkExposure<double, lsst::afw::image::maskPixelType>;

%template(darkCurrentCorrectChunkExposure)
    lsst::ip::isr::darkCurrentCorrectChunkExposure<float, lsst::afw::image::maskPixelType>;
%template(darkCurrentChunkExposure)
    lsst::ip::isr::darkCurrentCorrectChunkExposure<double, lsst::afw::image::maskPixelType>;

%template(linearizeChunkExposure)
    lsst::ip::isr::linearizeChunkExposure<float, lsst::afw::image::maskPixelType>;
%template(linearizeChunkExposure)
    lsst::ip::isr::linearizeChunkExposure<double, lsst::afw::image::maskPixelType>;

%template(flatFieldCorrectChunkExposure)
    lsst::ip::isr::flatFieldCorrectChunkExposure<float, lsst::afw::image::maskPixelType>;
%template(flatFieldCorrectChunkExposure)
    lsst::ip::isr::flatFieldCorrectChunkExposure<double, lsst::afw::image::maskPixelType>;

// %template(illuminationCorrection)
//     lsst::ip::isr::illuminationCorrection<float, lsst::afw::image::maskPixelType>;
// %template(illuminationCorrection)
//     lsst::ip::isr::illuminationCorrection<double, lsst::afw::image::maskPixelType>;

// %template(pupilImageCorrection)
//     lsst::ip::isr::pupilImageCorrection<float, lsst::afw::image::maskPixelType>;
// %template(pupilImageCorrection)
//     lsst::ip::isr::pupilImageCorrection<double, lsst::afw::image::maskPixelType>;

// %template(crosstalkCorrectChunkExposure)
//     lsst::ip::isr::ChunkExposure<float, lsst::afw::image::maskPixelType>;
// %template(defringeChunkExposure)
//     lsst::ip::isr::defringeChunkExposure<double, lsst::afw::image::maskPixelType>;

%template(defringeChunkExposure)
    lsst::ip::isr::defringeChunkExposure<float, lsst::afw::image::maskPixelType>;
%template(defringeChunkExposure)
    lsst::ip::isr::defringeChunkExposure<double, lsst::afw::image::maskPixelType>;

%template(interpolateOverMaskedPixels)
    lsst::ip::isr::interpolateOverMaskedPixels<float, lsst::afw::image::maskPixelType>;
%template(interpolateOverMaskedPixels)
    lsst::ip::isr::interpolateOverMaskedPixels<double, lsst::afw::image::maskPixelType>;

// %template(geometricDistortionCorrection)
//     lsst::ip::isr::geometricDistortionCorrection<float, lsst::afw::image::maskPixelType>;
// %template(geometricDistortionCorrection)
//     lsst::ip::isr::geometricDistortionCorrection<double, lsst::afw::image::maskPixelType>;

// %template(maskAndCorrectAdditionalArtifacts)
//     lsst::ip::isr::maskAndCorrectAdditionalArtifacts<float, lsst::afw::image::maskPixelType>;
// %template(maskAndCorrectAdditionalArtifacts)
//     lsst::ip::isr::maskAndCorrectAdditionalArtifacts<double, lsst::afw::image::maskPixelType>;

// %template(additionalFlatFieldCorrection)
//     lsst::ip::isr::additionalFlatFieldCorrection<float, lsst::afw::image::maskPixelType>;
// %template(additionalFlatFieldCorrection)
//     lsst::ip::isr::additionalFlatFieldCorrection<double, lsst::afw::image::maskPixelType>;


/******************************************************************************/
// Local Variables: ***
// eval: (setq indent-tabs-mode nil) ***
// End: ***
