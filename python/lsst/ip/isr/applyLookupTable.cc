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
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "numpy/arrayobject.h"
#include "ndarray/pybind11.h"

#include "lsst/ip/isr/applyLookupTable.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace ip {
namespace isr {

namespace {

template <typename PixelT>
static void declareApplyLookupTable(py::module& mod) {
    mod.def("applyLookupTable", &applyLookupTable<PixelT>, "image"_a, "table"_a, "indOffset"_a);
}

}  // namespace lsst::ip::isr::<anonymous>

PYBIND11_PLUGIN(applyLookupTable) {
    py::module mod("applyLookupTable");

    // Need to import numpy for ndarray and eigen conversions
    if (_import_array() < 0) {
        PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import");
        return nullptr;
    }

    declareApplyLookupTable<float>(mod);
    declareApplyLookupTable<double>(mod);

    return mod.ptr();
}
}
}
}  // lsst::ip::isr
