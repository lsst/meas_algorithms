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
#include <pybind11/pybind11.h>
//#include <pybind11/stl.h>

#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/meas/algorithms/PSF.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace algorithms {

PYBIND11_PLUGIN(psf) {
    py::module mod("psf");

    py::class_<PsfAttributes> clsPsfAttributes(mod, "PsfAttributes");

    /* Member types and enums */
    py::enum_<PsfAttributes::Method>(clsPsfAttributes, "Method")
        .value("ADAPTIVE_MOMENT", PsfAttributes::Method::ADAPTIVE_MOMENT)
        .value("FIRST_MOMENT", PsfAttributes::Method::FIRST_MOMENT)
        .value("SECOND_MOMENT", PsfAttributes::Method::SECOND_MOMENT)
        .value("NOISE_EQUIVALENT", PsfAttributes::Method::NOISE_EQUIVALENT)
        .value("BICKERTON", PsfAttributes::Method::BICKERTON)
        .export_values();

    /* Constructors */
    clsPsfAttributes.def(py::init<std::shared_ptr<lsst::afw::detection::Psf const>, int const, int const>(),
                         "psf"_a, "iX"_a, "iY"_a);
    clsPsfAttributes.def(
        py::init<std::shared_ptr<lsst::afw::detection::Psf const>, lsst::afw::geom::Point2I const &>(),
        "psf"_a, "cen"_a);

    /* Members */
    clsPsfAttributes.def("computeGaussianWidth", &PsfAttributes::computeGaussianWidth,
                         "how"_a = PsfAttributes::ADAPTIVE_MOMENT);
    clsPsfAttributes.def("computeEffectiveArea", &PsfAttributes::computeEffectiveArea);

    return mod.ptr();
}
}
}
}  // lsst::meas::algorithms
