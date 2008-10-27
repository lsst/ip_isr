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





/******************************************************************************/
// Local Variables: ***
// eval: (setq indent-tabs-mode nil) ***
// End: ***
