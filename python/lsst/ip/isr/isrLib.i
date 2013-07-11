// -*- lsst-c++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
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

#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging.h"
#include "lsst/ip/isr.h"
#include "lsst/pex/policy/Policy.h"
#include "lsst/afw/detection.h"
#include "lsst/afw/geom.h" // work around ticket #1121
#include "lsst/afw/cameraGeom.h"
%}

%include "lsst/p_lsstSwig.i"
%import  "lsst/afw/image/imageLib.i" 
%import  "lsst/afw/math/mathLib.i" 
%lsst_exceptions();

%shared_ptr(lsst::ip::isr::CountMaskedPixels<float>);
%shared_ptr(lsst::ip::isr::CountMaskedPixels<double>);

%shared_ptr(lsst::ip::isr::LookupTableMultiplicative<float>);
%shared_ptr(lsst::ip::isr::LookupTableMultiplicative<double>);

%shared_ptr(lsst::ip::isr::LookupTableReplace<float>);
%shared_ptr(lsst::ip::isr::LookupTableReplace<double>);

%include "lsst/ip/isr.h"

%template(CountMaskedPixelsF) lsst::ip::isr::CountMaskedPixels<float>;
%template(CountMaskedPixelsD) lsst::ip::isr::CountMaskedPixels<double>;

%template(LookupTableMultiplicativeF) lsst::ip::isr::LookupTableMultiplicative<float>;
%template(LookupTableMultiplicativeD) lsst::ip::isr::LookupTableMultiplicative<double>;

%template(LookupTableReplaceI) lsst::ip::isr::LookupTableReplace<int>;
%template(LookupTableReplaceF) lsst::ip::isr::LookupTableReplace<float>;

%template(fitOverscanImage) lsst::ip::isr::fitOverscanImage<float, double>;
%template(fitOverscanImage) lsst::ip::isr::fitOverscanImage<double, double>;


%define %instantiateUNC(TYPE, PIXELTYPE)
%newobject makeUnmaskedNanCounter;
%template(UnmaskedNanCounter##TYPE) lsst::ip::isr::UnmaskedNanCounter<PIXELTYPE>;
%template(makeUnmaskedNanCounter) lsst::ip::isr::makeUnmaskedNanCounter<PIXELTYPE>;
%enddef

%instantiateUNC(F, float);
%instantiateUNC(D, double);
%instantiateUNC(I, int);
%instantiateUNC(U, boost::uint16_t);

/******************************************************************************/
// Local Variables: ***
// eval: (setq indent-tabs-mode nil) ***
// End: ***
