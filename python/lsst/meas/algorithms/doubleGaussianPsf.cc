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
#include "lsst/meas/algorithms/DoubleGaussianPsf.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace algorithms {

PYBIND11_PLUGIN(doubleGaussianPsf) {
    py::module mod("doubleGaussianPsf");

    afw::table::io::python::declarePersistableFacade<DoubleGaussianPsf>(mod, "DoubleGaussianPsf");

    py::class_<DoubleGaussianPsf, std::shared_ptr<DoubleGaussianPsf>,
               afw::table::io::PersistableFacade<DoubleGaussianPsf>, KernelPsf>
            clsDoubleGaussianPsf(mod, "DoubleGaussianPsf");

    /* Constructors */
    clsDoubleGaussianPsf.def(py::init<int, int, double, double, double>(), "width"_a, "height"_a, "sigma1"_a,
                             "sigma2"_a = 0.0, "b"_a = 0.0);

    /* Members */
    clsDoubleGaussianPsf.def("clone", &DoubleGaussianPsf::clone);
    clsDoubleGaussianPsf.def("resized", &DoubleGaussianPsf::resized, "width"_a, "height"_a);
    clsDoubleGaussianPsf.def("getSigma1", &DoubleGaussianPsf::getSigma1);
    clsDoubleGaussianPsf.def("getSigma2", &DoubleGaussianPsf::getSigma2);
    clsDoubleGaussianPsf.def("getB", &DoubleGaussianPsf::getB);
    clsDoubleGaussianPsf.def("isPersistable", &DoubleGaussianPsf::isPersistable);

    return mod.ptr();
}

}  // algorithms
}  // meas
}  // lsst
