/*
 * LSST Data Management System
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 * See the COPYRIGHT file
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
 * see <https://www.lsstcorp.org/LegalNotices/>.
 */
#include <memory>

#include <pybind11/pybind11.h>

#include "lsst/ip/isr/isr.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace ip {
namespace isr {

namespace {

template <typename PixelT>
void declareCountMaskedPixels(py::module& mod, std::string const& suffix) {
    py::class_<CountMaskedPixels<PixelT>, std::shared_ptr<CountMaskedPixels<PixelT>>> cls(
        mod, ("CountMaskedPixels" + suffix).c_str());

    cls.def("reset", &CountMaskedPixels<PixelT>::reset);
    cls.def("apply", &CountMaskedPixels<PixelT>::apply, "image"_a, "bitmask"_a);
    cls.def("getCount", &CountMaskedPixels<PixelT>::getCount);
}

/**
 * Wrap all code in Isr.h for a given template parameter
 *
 * @tparam PixelT  Pixel type; typically `float` or `double` (potentially could also be
 *                  and integer class, but so far we have not needed those)
 * @param mod  pybind11 module to which to add the wrappers.
 * @param[in] suffix  Class name suffix associated with `PixelT`, e.g. "F" for `float` and "D" for `double`
 *
 * Note that the second (function type) template parameter of `fitOverscanImage` is always `double`.
 */
template <typename PixelT>
void declareAll(py::module& mod, std::string const& suffix) {
    declareCountMaskedPixels<PixelT>(mod, suffix);

    mod.def("maskNans", &maskNans<PixelT>, "maskedImage"_a, "maskVal"_a, "allow"_a = 0);
    mod.def("fitOverscanImage", &fitOverscanImage<PixelT, double>, "overscanFunction"_a, "overscan"_a,
            "stepSize"_a = 1.1, "sigma"_a = 1);
}

}  // namespace lsst::ip::isr::<anonymous>

PYBIND11_PLUGIN(_isr) {
    py::module mod("_isr", "Python wrapper for _isr library");

    declareAll<float>(mod, "F");
    declareAll<double>(mod, "D");

    return mod.ptr();
}
}
}
}  // lsst::ip::isr
