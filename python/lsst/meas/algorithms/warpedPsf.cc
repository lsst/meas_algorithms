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

namespace py = pybind11;
using namespace pybind11::literals;

#include "lsst/meas/algorithms/WarpedPsf.h"

namespace lsst {
namespace meas {
namespace algorithms {

PYBIND11_PLUGIN(_warpedPsf) {
    py::module mod("_warpedPsf", "Python wrapper for afw _warpedPsf library");

    py::class_<WarpedPsf, std::shared_ptr<WarpedPsf>, ImagePsf> clsWarpedPsf(mod, "WarpedPsf");

    /* Constructors */
    clsWarpedPsf.def(
        py::init<std::shared_ptr<afw::detection::Psf const>, std::shared_ptr<afw::geom::XYTransform const>,
                 std::shared_ptr<afw::math::WarpingControl const>>(),
        "undistortedPsf"_a, "distortion"_a, "control"_a);
    clsWarpedPsf.def(
        py::init<std::shared_ptr<afw::detection::Psf const>, std::shared_ptr<afw::geom::XYTransform const>,
                 std::string const&, unsigned int>(),
        "undistortedPsf"_a, "distortion"_a, "kernelName"_a = "lanczos3", "cache"_a = 10000);

    /* Members */
    clsWarpedPsf.def("getAveragePosition", &WarpedPsf::getAveragePosition);
    clsWarpedPsf.def("clone", &WarpedPsf::clone);

    return mod.ptr();
}
}
}
}  // lsst::meas::algorithms
