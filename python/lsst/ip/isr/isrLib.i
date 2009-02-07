// -*- lsst-c++ -*-
%define isrLib_DOCSTRING
"
Python bindings for lsst::ip::isr Instrument Signature Removal code
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.ip.isr", docstring=isrLib_DOCSTRING) isrLib

// Suppress swig complaints; see afw/image/imageLib.i for more 
#pragma SWIG nowarn=362                 // operator=  ignored 

// Everything we will need in the _wrap.cc file
%{
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/pex/policy/Policy.h"
#include "lsst/afw/detection.h"
#include "lsst/ip/isr/isr.h"
%}

%init %{
%}

// namespace boost {
//     class bad_any_cast; // remove warning: Nothing known about 'boost::bad_any_cast'
// }

// Everything whose bindings we will have to know about
%import "lsst/p_lsstSwig.i"
%import "lsst/afw/image/imageLib.i"
%import "lsst/afw/math/mathLib.i"
//%import "lsst/afw/detection/detectionLib.i"

%lsst_exceptions();

// Here you give the names of the functions you want to Swig

// handle C++ arguments that should be outputs in python
// %apply int& OUTPUT { int& };
// %apply float& OUTPUT { float& };
// %apply double& OUTPUT { double& };

%include "lsst/ip/isr/isr.h"

%template(fitFunctionToImage) lsst::ip::isr::fitFunctionToImage<float>;
%template(iterateTable) lsst::ip::isr::fitFunctionToImage<double>;

%template(iterateTable) lsst::ip::isr::iterateTable<float>;
%template(iterateTable) lsst::ip::isr::iterateTable<double>;

%template(findBestFit) lsst::ip::isr::findBestFit<float>;
%template(findBestFit) lsst::ip::isr::findBestFit<double>;

// %template(saturationCorrectionForChunkExposure)
//     lsst::ip::isr::saturationCorrectionForChunkExposure<float>;
// %template(saturationCorrectionForChunkExposure)
//     lsst::ip::isr::saturationCorrectionForChunkExposure<double>;

// %template(overscanCorrectAndTrimChunkExposure)
//     lsst::ip::isr::overscanCorrectAndTrimChunkExposure<float>;
// %template(overscanCorrectAndTrimChunkExposure)
//     lsst::ip::isr::overscanCorrectAndTrimChunkExposure<double>;

// %template(biasCorrectChunkExposure)
//     lsst::ip::isr::biasCorrectChunkExposure<float>;
// %template(biasCorrectChunkExposure)
//     lsst::ip::isr::biasCorrectChunkExposure<double>;

// %template(darkCurrentCorrectChunkExposure)
//     lsst::ip::isr::darkCurrentCorrectChunkExposure<float>;
// %template(darkCurrentChunkExposure)
//     lsst::ip::isr::darkCurrentCorrectChunkExposure<double>;

// %template(linearizeChunkExposure)
//     lsst::ip::isr::linearizeChunkExposure<float>;
// %template(linearizeChunkExposure)
//     lsst::ip::isr::linearizeChunkExposure<double>;

// %template(flatFieldCorrectChunkExposure)
//     lsst::ip::isr::flatFieldCorrectChunkExposure<float>;
// %template(flatFieldCorrectChunkExposure)
//     lsst::ip::isr::flatFieldCorrectChunkExposure<double>;

// %template(illuminationCorrection)
//     lsst::ip::isr::illuminationCorrection<float>;
// %template(illuminationCorrection)
//     lsst::ip::isr::illuminationCorrection<double>;

// %template(illuminationCorrectionDR)
//     lsst::ip::isr::illuminationCorrectionDR<float>;
// %template(illuminationCorrectionDR)
//     lsst::ip::isr::illuminationCorrectionDR<double>;

// %template(pupilImageCorrection)
//     lsst::ip::isr::pupilImageCorrection<float>;
// %template(pupilImageCorrection)
//     lsst::ip::isr::pupilImageCorrection<double>;

// %template(crosstalkCorrectChunkExposure)
//     lsst::ip::isr::crosstalkCorrectChunkExposure<float>;
// %template(crosstalkCorrectChunkExposure)
//     lsst::ip::isr::crosstalkCorrectChunkExposure<double>;

// %template(defringeChunkExposure)
//     lsst::ip::isr::defringeChunkExposure<float>;
// %template(defringeChunkExposure)
//     lsst::ip::isr::defringeChunkExposure<double>;

// %template(geometricDistortionCorrection)
//     lsst::ip::isr::geometricDistortionCorrection<float>;
// %template(geometricDistortionCorrection)
//     lsst::ip::isr::geometricDistortionCorrection<double>;

// %template(maskAndCorrectAdditionalArtifacts)
//     lsst::ip::isr::maskAndCorrectAdditionalArtifacts<float>;
// %template(maskAndCorrectAdditionalArtifacts)
//     lsst::ip::isr::maskAndCorrectAdditionalArtifacts<double>;

// %template(additionalFlatFieldCorrection)
//     lsst::ip::isr::additionalFlatFieldCorrection<float>;
// %template(additionalFlatFieldCorrection)
//     lsst::ip::isr::additionalFlatFieldCorrection<double>;

// %include "lsst/ip/isr/interpolateOverMaskedPixels.h"

// %template(interpolateOverMaskedPixels)
//     lsst::ip::isr::interpolateOverMaskedPixels<float>;
// %template(interpolateOverMaskedPixels)
//     lsst::ip::isr::interpolateOverMaskedPixels<double>;

/******************************************************************************/
// Local Variables: ***
// eval: (setq indent-tabs-mode nil) ***
// End: ***
