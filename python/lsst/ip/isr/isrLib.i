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
#include <memory>

#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging.h"
#include "lsst/ip/isr/applyLookupTable.h"
#include "lsst/ip/isr/isr.h"
#include "lsst/pex/policy/Policy.h"
#include "lsst/afw/detection.h"
#include "lsst/afw/cameraGeom.h"
%}

%include "lsst/p_lsstSwig.i"
%initializeNumPy(ip_isr)
%{
#include "ndarray/swig.h"
%}

%import  "lsst/afw/image/imageLib.i" 
%import  "lsst/afw/math/mathLib.i" 
%lsst_exceptions();

%shared_ptr(lsst::ip::isr::CountMaskedPixels<float>);
%shared_ptr(lsst::ip::isr::CountMaskedPixels<double>);

%include "lsst/ip/isr/applyLookupTable.h"
%include "lsst/ip/isr/isr.h"

%define %instantiateFloatLike(TYPE, PIXELTYPE)
%template(applyLookupTable) lsst::ip::isr::applyLookupTable<PIXELTYPE>;
%template(CountMaskedPixels##TYPE) lsst::ip::isr::CountMaskedPixels<PIXELTYPE>;
%template(fitOverscanImage) lsst::ip::isr::fitOverscanImage<PIXELTYPE, double>;
%template(maskNans) lsst::ip::isr::maskNans<PIXELTYPE>;
%enddef

%instantiateFloatLike(F, float);
%instantiateFloatLike(D, double);

/******************************************************************************/
// Local Variables: ***
// eval: (setq indent-tabs-mode nil) ***
// End: ***
