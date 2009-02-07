// -*- LSST-C++ -*- // fixed format comment for emacs
/**
  * \file
  *
  * \ingroup isr
  *
  * \brief Utility function to parse strings for fits header data sections.
  *
  * \author Nicole M. Silvestri, University of Washington
  *
  * Contact: nms@astro.washington.edu
  *
  * \version
  *
  * LSST Legalese here...
  */

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

#include "lsst/ip/isr/isr.h"

namespace afwImage = lsst::afw::image;
namespace ipIsr = lsst::ip::isr;

// Simple string iterator to help parse data sections.
// Finds the string between two delimiting characters
// Adapted from C-code snippet by Mike Wahler.

std::string between(std::string &s, char ldelim, char rdelim)
{
    std::string::iterator b(s.begin());
    std::string::iterator e(s.end());
    std::string::iterator lp;
    std::string::iterator rp;

    std::string result;

    if((lp = std::find(b, e, ldelim)) != e)
        if((rp = std::find(++lp, e, rdelim)) != e)
            result = std::string(lp, rp);

    return result;
}

afwImage::BBox ipIsr::stringParse(std::string &section)
{ 

    const char begin('[');
    const char end(']');
    const char delim1(':');
    const char delim2(',');

    std::string temp(between(section, delim2, end));
    std::size_t colonPos = temp.find(":");

    // NOTE: atoi() needs to be passed a c_str() to get the int out
    afwImage::PointI startPt(
        atoi(between(section, begin, delim1).c_str()),
        atoi(between(section, delim2, delim1).c_str()));
    std::cout << "start = " << startPt.getX() << ", " << startPt.getY() << std::endl;


    afwImage::PointI endPt(
        atoi(between(section, delim1, delim2).c_str()),
        atoi(temp.substr(colonPos + 1).c_str()));
    std::cout << "end = " << endPt.getX() << ", " << endPt.getY() << std::endl;

    return afwImage::BBox(startPt, endPt);
}
