// -*- lsst-c++ -*-
%define isrLib_DOCSTRING
"
Python bindings for lsst::ip::isr Instrument Signature Removal code
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.ip.isr", docstring=isrLib_DOCSTRING) isrLib

// Suppress swig complaints; see afw/image/imageLib.i for more 
#pragma SWIG nowarn=314                 // print is a python keyword (--> _print)
#pragma SWIG nowarn=362                 // operator=  ignored 

// Everything we will need in the _wrap.cc file
%{
#include <boost/shared_ptr.hpp>

#include <lsst/afw/image.h>
#include <lsst/afw/math.h>
#include <lsst/afw/math/Statistics.h>

#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/pex/policy/Policy.h"
#include "lsst/afw/detection.h"
%}

%include "lsst/p_lsstSwig.i"
%import  "lsst/afw/image/imageLib.i" 
%import  "lsst/afw/math/mathLib.i" 
%lsst_exceptions();

%{
#include "lsst/ip/isr/Isr.h"
%}

SWIG_SHARED_PTR(LookupTableMultiplicativeF, lsst::ip::isr::LookupTableMultiplicative<float>);
SWIG_SHARED_PTR(LookupTableMultiplicativeD, lsst::ip::isr::LookupTableMultiplicative<double>);

SWIG_SHARED_PTR(LookupTableReplaceF, lsst::ip::isr::LookupTableReplace<float>);
SWIG_SHARED_PTR(LookupTableReplaceD, lsst::ip::isr::LookupTableReplace<double>);

%include "lsst/ip/isr/Isr.h"

%template(LookupTableMultiplicativeF) lsst::ip::isr::LookupTableMultiplicative<float>;
%template(LookupTableMultiplicativeD) lsst::ip::isr::LookupTableMultiplicative<double>;

%template(LookupTableReplaceI) lsst::ip::isr::LookupTableReplace<int>;

%template(fitOverscanImage) lsst::ip::isr::fitOverscanImage<float, double>;
%template(fitOverscanImage) lsst::ip::isr::fitOverscanImage<double, double>;

/******************************************************************************/
// Local Variables: ***
// eval: (setq indent-tabs-mode nil) ***
// End: ***
