// -*- LSST-C++ -*- 
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

#include <lsst/afw/math.h>
#include <lsst/afw/image.h>
#include <lsst/pex/exceptions/Exception.h>
#include <lsst/pex/logging/Trace.h>

/** \brief Remove all non-astronomical counts from the Chunk Exposure's pixels.
  * 
  */
	
namespace lsst {
namespace ip {
namespace isr {

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

        void apply(lsst::afw::image::MaskedImage<ImageT> &image) {

            for (int y = 0; y != image.getHeight(); ++y) {
                for (x_iterator ptr = image.row_begin(y); ptr != image.row_end(y); ++ptr) {
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

        void apply(lsst::afw::image::MaskedImage<ImageT> &image, double gain=1.0) {
            for (int y = 0; y != image.getHeight(); ++y) {
                for (x_iterator ptr = image.row_begin(y); ptr != image.row_end(y); ++ptr) {
                    int ind = static_cast<int>(ptr.image() + 0.5);  // Rounded pixel value
                    if (ind >= _max){
                        throw LSST_EXCEPT(lsst::pex::exceptions::Exception, 
                                          "Pixel value out of range in LookupTableReplace::apply");
                    }
                    PixelT p = PixelT(_table[ind], 
                                      (*ptr).mask(), 
                                      _table[ind] * gain);
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

    template<typename ImagePixelT>
    void fitFunctionToImage(
        lsst::afw::image::MaskedImage<ImagePixelT> &maskedImage,
        lsst::afw::math::Function1<double> const &function
        );

    

}}} // namespace lsst::ip::isr
	
#endif // !defined(LSST_IP_ISR_ISR_H)
