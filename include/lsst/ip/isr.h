// -*- LSST-C++ -*- 

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
 
/**
  * \file
  *
  * \ingroup isr
  *
  * \brief Implementation of the templated Instrument Signature Removal 
  * stage of the nightly LSST Image Processing Pipeline.
  *
  * \author Nicole M. Silvestri / ACB, University of Washington
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

#include <boost/shared_ptr.hpp>

#include <lsst/afw/math.h>
#include <lsst/afw/math/Statistics.h>
#include <lsst/afw/image.h>
#include <lsst/utils/ieee.h>
#include <lsst/pex/exceptions/Exception.h>

/** \brief Remove all non-astronomical counts from the Chunk Exposure's pixels.
  * 
  */
	
namespace lsst {
namespace ip {
namespace isr {

    /** Signature strings associated with each stage of the ISR
     * 
     * @note Added to Exposure metadata after stage processing
     *
     * @note As implementation detail, no more than 8 characters for fits
     * compliance?
     */
    std::string const& ISR_LIN     = "ISR_LIN";    ///< Linearization
    std::string const& ISR_OSCAN   = "ISR_OSCAN";  ///< Overscan
    std::string const& ISR_TRIM    = "ISR_TRIM";   ///< Trim
    std::string const& ISR_BIAS    = "ISR_BIAS";   ///< Bias 
    std::string const& ISR_DFLAT   = "ISR_DFLAT";  ///< Dome flat
    std::string const& ISR_ILLUM   = "ISR_ILLUM";  ///< Illumination correction
    std::string const& ISR_BADP    = "ISR_BADP";   ///< Bad pixel mask
    std::string const& ISR_SAT     = "ISR_SAT";    ///< Saturated pixels
    std::string const& ISR_FRING   = "ISR_FRING";  ///< Fringe correction
    std::string const& ISR_DARK    = "ISR_DARK";   ///< Dark correction
    std::string const& ISR_PUPIL   = "ISR_PUPIL";  ///< Pupil correction
    std::string const& ISR_CRREJ   = "ISR_CRREJ";  ///< Cosmic ray rejection
    std::string const& ISR_BACKSUB = "ISR_BACKSUB";  ///< Background subtraction

    enum StageId {
        ISR_LINid   = 0x1,   ///< Linearization
        ISR_OSCANid = 0x2,   ///< Overscan
        ISR_TRIMid  = 0x4,   ///< Trim
        ISR_BIASid  = 0x8,   ///< Bias 
        ISR_DFLATid = 0x10,  ///< Dome flat
        ISR_ILLUMid = 0x20,  ///< Illumination correction
        ISR_BADPid  = 0x40,  ///< Bad pixel mask
        ISR_SATid   = 0x80,  ///< Saturated pixels
        ISR_FRINid  = 0x100, ///< Fringe correction
        ISR_DARKid  = 0x200, ///< Dark correction
        ISR_PUPILid = 0x400, ///< Pupil correction
        ISR_CRREJid = 0x800, ///< Cosmic ray rejection
        ISR_BACKSUBid = 0x1000, ///< Cosmic ray rejection
    };

    /** Multiplicative linearization lookup table
     *
     * @ingroup isr
     */
    template <typename ImageT>
    class LookupTableMultiplicative {
    public:
        typedef typename lsst::afw::image::MaskedImage<ImageT>::x_iterator x_iterator;
        typedef typename lsst::afw::image::MaskedImage<ImageT>::Pixel PixelT;

        LookupTableMultiplicative(std::vector<double> table) : 
            _table(table), _max(table.size()) {};
        virtual ~LookupTableMultiplicative() {};

        void apply(lsst::afw::image::MaskedImage<ImageT> &image, float gain=1.0) {

            for (int y = 0; y != image.getHeight(); ++y) {
                for (x_iterator ptr = image.row_begin(y), end = image.row_end(y); ptr != end; ++ptr) {
                    int ind = static_cast<int>(ptr.image() + 0.5);  // Rounded pixel value
                    if (ind >= _max){
                        throw LSST_EXCEPT(lsst::pex::exceptions::Exception, 
                                          "Pixel value out of range in LookupTableMultiplicative::apply");
                    }
                    PixelT p = PixelT((*ptr).image() * _table[ind], 
                                      (*ptr).mask(), 
                                      (*ptr).variance() * _table[ind] * _table[ind]);
                    *ptr = p;
                }
            }
        }

        // Return the lookup table
        std::vector<double> getTable() const { return _table; }
    private:
        std::vector<double> _table;
        int _max;
    };

    /** Linearization lookup table with replacement
     *
     * @ingroup isr
     */
    template <typename ImageT>
    class LookupTableReplace {
    public:
        typedef typename lsst::afw::image::MaskedImage<ImageT>::x_iterator x_iterator;
        typedef typename lsst::afw::image::MaskedImage<ImageT>::Pixel PixelT;

        LookupTableReplace(std::vector<double> table) : 
            _table(table), _max(table.size()) {};
        virtual ~LookupTableReplace() {};

        void apply(lsst::afw::image::MaskedImage<ImageT> &image, float gain=1.0) const;

        // Return the lookup table
        std::vector<double> getTable() const { return _table; }
    private:
        std::vector<double> _table;
        int _max;
    };

    template <typename ImageT, typename MaskT=lsst::afw::image::MaskPixel>
    class CountMaskedPixels {
    public:
        typedef typename lsst::afw::image::MaskedImage<ImageT>::x_iterator x_iterator;
        CountMaskedPixels() : 
            _count(0) {} ;
        virtual ~CountMaskedPixels() {};

        // Clear the accumulator
        void reset() { _count = 0; }

        // Count pixels
        void apply(lsst::afw::image::MaskedImage<ImageT> const& image,
                   MaskT bitmask) {
            reset();
            for (int y = 0; y != image.getHeight(); ++y) {
                for (x_iterator ptr = image.row_begin(y); ptr != image.row_end(y); ++ptr) {
                    if ( ((*ptr).mask() & bitmask) == bitmask ) {
                        _count += 1;
                    }
                }
            }
        }

        // Return the total counts
        int getCount() const { return _count; }
        
    private:
        int _count;
    };

    /// Mask NANs in an image
    ///
    /// NANs in the image or variance that are not already masked by
    /// the 'allow' value are masked with the 'maskVal'.
    ///
    /// @return Number of pixels masked
    template <typename PixelT>
    size_t maskNans(
        afw::image::MaskedImage<PixelT> const& mi, ///< Input image
        afw::image::MaskPixel maskVal,  ///< Bit mask value to give a NaN
        afw::image::MaskPixel allow=0 ///< Allow NANs with this bit mask (0 to disallow all NANs)
        );


    template<typename ImagePixelT, typename FunctionT>
    void fitOverscanImage(
        boost::shared_ptr<lsst::afw::math::Function1<FunctionT> > &overscanFunction,
        lsst::afw::image::MaskedImage<ImagePixelT> const& overscan,
        double ssize=1.,
        int sigma=1
        );
    

}}} // namespace lsst::ip::isr
	
#endif // !defined(LSST_IP_ISR_ISR_H)
