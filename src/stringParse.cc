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

#include "vw/Math/BBox.h"
#include "lsst/ip/isr/stringParse.h"

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

vw::BBox2i lsst::ip::isr::stringParse(std::string &section)
{ 

    const char begin('[');
    const char end(']');
    const char delim1(':');
    const char delim2(',');

    std::string temp(between(section, delim2, end));
    std::size_t position = temp.find(":");

    // NOTE: atoi() needs to be passed a c_str() to get the int out
    const int colsStart = atoi(between(section, begin, delim1).c_str());
    std::cout << "colsStart: " << colsStart << std::endl;

    const int colsEnd = atoi(between(section, delim1, delim2).c_str());
    std::cout << "colsEnd: " << colsEnd << std::endl;

    const int rowsStart = atoi(between(section, delim2, delim1).c_str());
    std::cout << "rowsStart: " << rowsStart << std::endl;

    const int rowsEnd = atoi(temp.substr(position + 1).c_str());
    std::cout << "rowsEnd: " << rowsEnd << std::endl;

    const int colSpan = colsEnd - colsStart;
    const int rowSpan = rowsEnd - rowsStart;

    vw::BBox2i bBox = vw::BBox2i(colsStart, rowsStart, colSpan, rowSpan);

    return bBox;
}
