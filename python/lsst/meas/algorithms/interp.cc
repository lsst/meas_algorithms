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
#include <pybind11/stl.h>

#include "lsst/afw/detection/Psf.h"
#include "lsst/meas/algorithms/Interp.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace algorithms {

namespace {
template <typename PixelT>
void declareInterpolateOverDefects(py::module& mod) {
    mod.def("interpolateOverDefects",
            interpolateOverDefects<lsst::afw::image::MaskedImage<PixelT, lsst::afw::image::MaskPixel,
                                                                 lsst::afw::image::VariancePixel>>,
            "image"_a, "psf"_a, "badList"_a, "fallBackValue"_a = 0.0, "useFallbackValueAtEdge"_a = false);
}
}

PYBIND11_PLUGIN(interp) {
    py::module mod("interp");

    py::class_<Defect, std::shared_ptr<Defect>, lsst::afw::image::DefectBase> clsDefect(mod, "Defect");

    /* Member types and enums */
    py::enum_<Defect::DefectPosition>(clsDefect, "DefectPosition")
        .value("LEFT", Defect::DefectPosition::LEFT)
        .value("NEAR_LEFT", Defect::DefectPosition::NEAR_LEFT)
        .value("WIDE_LEFT", Defect::DefectPosition::WIDE_LEFT)
        .value("MIDDLE", Defect::DefectPosition::MIDDLE)
        .value("WIDE_NEAR_LEFT", Defect::DefectPosition::WIDE_NEAR_LEFT)
        .value("WIDE", Defect::DefectPosition::WIDE)
        .value("WIDE_NEAR_RIGHT", Defect::DefectPosition::WIDE_NEAR_RIGHT)
        .value("NEAR_RIGHT", Defect::DefectPosition::NEAR_RIGHT)
        .value("WIDE_RIGHT", Defect::DefectPosition::WIDE_RIGHT)
        .value("RIGHT", Defect::DefectPosition::RIGHT)
        .export_values();

    /* Constructors */
    clsDefect.def(py::init<const lsst::afw::geom::BoxI&>(), "bbox"_a = lsst::afw::geom::BoxI());

    /* Members */
    clsDefect.def("classify", &Defect::classify);
    clsDefect.def("getType", &Defect::getType);
    clsDefect.def("getPos", &Defect::getPos);

    declareInterpolateOverDefects<float>(mod);

    return mod.ptr();
}
}
}
}  // lsst::meas::algorithms
