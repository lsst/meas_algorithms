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

#include "lsst/afw/table/io/pybind11.h"
#include "lsst/meas/algorithms/CoaddPsf.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace algorithms {

PYBIND11_PLUGIN(_coaddPsf) {
    py::module mod("_coaddPsf", "Python wrapper for afw _coaddPsf library");

    afw::table::io::declarePersistableFacade<CoaddPsf>(mod, "CoaddPsf");

    py::class_<CoaddPsf, std::shared_ptr<CoaddPsf>, afw::table::io::PersistableFacade<CoaddPsf>, ImagePsf>
        clsCoaddPsf(mod, "CoaddPsf");

    /* Constructors */
    clsCoaddPsf.def(py::init<afw::table::ExposureCatalog const &, afw::image::Wcs const &,
                             std::string const &, std::string const &, int>(),
                    "catalog"_a, "coaddWcs"_a, "weightFieldName"_a = "weight",
                    "warpingKernelName"_a = "lanczos3", "cacheSize"_a = 10000);

    /* Members */
    clsCoaddPsf.def("clone", &CoaddPsf::clone);
    clsCoaddPsf.def("getAveragePosition", &CoaddPsf::getAveragePosition);
    clsCoaddPsf.def("getCoaddWcs", &CoaddPsf::getCoaddWcs);
    clsCoaddPsf.def("getComponentCount", &CoaddPsf::getComponentCount);
    clsCoaddPsf.def("getPsf", &CoaddPsf::getPsf);
    clsCoaddPsf.def("getWcs", &CoaddPsf::getWcs);
    clsCoaddPsf.def("getWeight", &CoaddPsf::getWeight);
    clsCoaddPsf.def("getId", &CoaddPsf::getId);
    clsCoaddPsf.def("getBBox", &CoaddPsf::getBBox);
    clsCoaddPsf.def("getValidPolygon", &CoaddPsf::getValidPolygon);
    clsCoaddPsf.def("isPersistable", &CoaddPsf::isPersistable);

    return mod.ptr();
}
}
}
}  // lsst::meas::algorithms
