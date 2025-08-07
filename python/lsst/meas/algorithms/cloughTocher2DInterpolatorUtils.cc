// -*- LSST-C++ -*-
/*
 * This file is part of meas_algorithms.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
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
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "pybind11/pybind11.h"
#include "ndarray/pybind11.h"
#include "lsst/cpputils/python.h"
#include "pybind11/stl.h"

#include "lsst/geom/Box.h"
#include "lsst/meas/algorithms/CloughTocher2DInterpolatorUtils.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace algorithms {
namespace {

void declareCloughTocher2DInterpolatorUtils(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyCloughTocher2DInterpolatorUtils =
            py::classh<CloughTocher2DInterpolatorUtils>;
    auto clsCloughTocher2DInterpolatorUtils = wrappers.wrapType(
            PyCloughTocher2DInterpolatorUtils(wrappers.module, "CloughTocher2DInterpolatorUtils"),
            [](auto &mod, auto &cls) {
                cls.def_static("findGoodPixelsAroundBadPixels",
                               &CloughTocher2DInterpolatorUtils::findGoodPixelsAroundBadPixels, "mimage"_a,
                               "badList"_a, "buffer"_a = 4);
                cls.def_static("updateArrayFromImage", &CloughTocher2DInterpolatorUtils::updateArrayFromImage,
                               "array"_a, "image"_a);
                cls.def_static("updateImageFromArray", &CloughTocher2DInterpolatorUtils::updateImageFromArray,
                               "image"_a, "array"_a);
            });
}
}  // namespace

void wrapCloughTocher2DInterpolatorUtils(lsst::cpputils::python::WrapperCollection &wrappers) {
    declareCloughTocher2DInterpolatorUtils(wrappers);
}

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst
