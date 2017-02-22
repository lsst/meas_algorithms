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

#include "lsst/afw/table/io/python.h"
#include "lsst/meas/algorithms/KernelPsf.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace algorithms {

PYBIND11_PLUGIN(kernelPsf) {
    py::module mod("kernelPsf");

    afw::table::io::python::declarePersistableFacade<KernelPsf>(mod, "KernelPsf");

    py::class_<KernelPsf, std::shared_ptr<KernelPsf>, afw::table::io::PersistableFacade<KernelPsf>, ImagePsf>
            clsKernelPsf(mod, "KernelPsf");

    /* Constructors */
    clsKernelPsf.def(py::init<afw::math::Kernel const &, afw::geom::Point2D const &>(), "kernel"_a,
                     "averagePosition"_a = afw::geom::Point2D());

    /* Members */
    clsKernelPsf.def("getKernel", &KernelPsf::getKernel);
    clsKernelPsf.def("getAveragePosition", &KernelPsf::getAveragePosition);
    clsKernelPsf.def("clone", &KernelPsf::clone);
    clsKernelPsf.def("isPersistable", &KernelPsf::isPersistable);

    return mod.ptr();
}

}  // algorithms
}  // meas
}  // lsst

