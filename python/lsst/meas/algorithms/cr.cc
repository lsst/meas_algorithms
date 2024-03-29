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
#include "lsst/cpputils/python.h"
#include "pybind11/stl.h"

#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/meas/algorithms/CR.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace algorithms {
namespace {

template <typename PixelT>
void declareFindCosmicRays(lsst::cpputils::python::WrapperCollection &wrappers) {
    wrappers.module.def("findCosmicRays", &findCosmicRays<afw::image::MaskedImage<PixelT>>, "image"_a, "psf"_a, "bkgd"_a,
            "policy"_a, "keep"_a = false);
}
}  // namespace
void wrapCr(lsst::cpputils::python::WrapperCollection &wrappers) {
    declareFindCosmicRays<float>(wrappers);
}

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst
