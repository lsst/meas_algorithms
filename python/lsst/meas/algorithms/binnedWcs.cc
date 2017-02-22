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

#include "lsst/meas/algorithms/BinnedWcs.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace algorithms {

PYBIND11_PLUGIN(binnedWcs) {
    py::module mod("binnedWcs");

    py::class_<BinnedWcs, std::shared_ptr<BinnedWcs>, afw::image::Wcs> clsBinnedWcs(mod, "BinnedWcs");

    /* Constructors */
    clsBinnedWcs.def(
        py::init<std::shared_ptr<afw::image::Wcs>, unsigned int, unsigned int, afw::geom::Point2I>(),
        "parent"_a, "xBin"_a, "yBin"_a, "xy0"_a);

    /* Members */
    clsBinnedWcs.def("clone", &BinnedWcs::clone);
    clsBinnedWcs.def("upcast", &BinnedWcs::upcast);
    clsBinnedWcs.def("getParent", &BinnedWcs::getParent);
    clsBinnedWcs.def("getXBin", &BinnedWcs::getXBin);
    clsBinnedWcs.def("getYBin", &BinnedWcs::getYBin);
    clsBinnedWcs.def("getXY0", &BinnedWcs::getXY0);
    clsBinnedWcs.def("hasDistortion", &BinnedWcs::hasDistortion);
    clsBinnedWcs.def("isPersistable", &BinnedWcs::isPersistable);
    clsBinnedWcs.def("flipImage", &BinnedWcs::flipImage);
    clsBinnedWcs.def("rotateImageBy90", &BinnedWcs::rotateImageBy90);
    clsBinnedWcs.def("getFitsMetadata", &BinnedWcs::getFitsMetadata);
    clsBinnedWcs.def("getBinnedToOriginal", &BinnedWcs::getBinnedToOriginal);
    clsBinnedWcs.def("getOriginalToBinned", &BinnedWcs::getOriginalToBinned);

    return mod.ptr();
}
}
}
}  // lsst::meas::algorithms
