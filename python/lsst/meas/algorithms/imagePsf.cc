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

#include "lsst/cpputils/python/PySharedPtr.h"
#include "lsst/afw/table/io/python.h"
#include "lsst/meas/algorithms/ImagePsf.h"
#include "lsst/meas/algorithms/python.h"

namespace py = pybind11;
using namespace pybind11::literals;

using lsst::cpputils::python::PySharedPtr;

namespace lsst {
namespace meas {
namespace algorithms {
namespace {

void declareImagePsf(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyImagePsf = py::class_<ImagePsf, PySharedPtr<ImagePsf>, afw::detection::Psf, ImagePsfTrampoline<>>;
    auto clsImagePsf = wrappers.wrapType(PyImagePsf(wrappers.module, "ImagePsf"), [](auto &mod, auto &cls) {
        cls.def(py::init<bool>(), "init", "isFixed"_a = false);  // Ctor for pure python subclasses
    });
    afw::table::io::python::addPersistableMethods<ImagePsf>(clsImagePsf);
}
}  // namespace

void wrapImagePsf(lsst::cpputils::python::WrapperCollection &wrappers) {
    declareImagePsf(wrappers);
}
}  // namespace algorithms
}  // namespace meas
}  // namespace lsst
