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
#include "pybind11/stl.h"

#include "lsst/geom/Box.h"
#include "lsst/afw/table/io/python.h"
#include "lsst/meas/algorithms/CoaddBoundedField.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace algorithms {

PYBIND11_PLUGIN(coaddBoundedField) {
    py::module mod("coaddBoundedField");

    py::class_<CoaddBoundedFieldElement> clsCoaddBoundedFieldElement(mod, "CoaddBoundedFieldElement");

    clsCoaddBoundedFieldElement.def(
        "__init__",
        [](CoaddBoundedFieldElement &instance,
            std::shared_ptr<afw::math::BoundedField> field,
            std::shared_ptr<afw::geom::SkyWcs const> wcs,
            py::object polygon,
            double weight) {
            if (polygon == py::none()) {
                    new (&instance) CoaddBoundedFieldElement(field, wcs, nullptr, weight);
                } else {
                    auto pgon = py::cast<std::shared_ptr<afw::geom::polygon::Polygon const>>(polygon);
                    new (&instance) CoaddBoundedFieldElement(field, wcs, pgon, weight);
                }
            },
            "field"_a, "wcs"_a, "validPolygon"_a, "weight"_a = 1.0);

    clsCoaddBoundedFieldElement.def_readwrite("field", &CoaddBoundedFieldElement::field);
    clsCoaddBoundedFieldElement.def_readwrite("wcs", &CoaddBoundedFieldElement::wcs);
    clsCoaddBoundedFieldElement.def_readwrite("validPolygon", &CoaddBoundedFieldElement::validPolygon);
    clsCoaddBoundedFieldElement.def_readwrite("weight", &CoaddBoundedFieldElement::weight);

    afw::table::io::python::declarePersistableFacade<CoaddBoundedField>(mod, "CoaddBoundedField");

    py::class_<CoaddBoundedField, std::shared_ptr<CoaddBoundedField>,
               afw::table::io::PersistableFacade<CoaddBoundedField>, afw::math::BoundedField>
            clsCoaddBoundedField(mod, "CoaddBoundedField");

    clsCoaddBoundedField.attr("Element") = clsCoaddBoundedFieldElement;

    /* Constructors */
    clsCoaddBoundedField.def(py::init<geom::Box2I const &, std::shared_ptr<afw::geom::SkyWcs const>,
                                      typename CoaddBoundedField::ElementVector const &>(),
                             "bbox"_a, "coaddWcs"_a, "elements"_a);
    clsCoaddBoundedField.def(py::init<geom::Box2I const &, std::shared_ptr<afw::geom::SkyWcs const>,
                                      typename CoaddBoundedField::ElementVector const &, double>(),
                             "bbox"_a, "coaddWcs"_a, "elements"_a, "default"_a);

    /* Operators */
    clsCoaddBoundedField.def("__eq__", &CoaddBoundedField::operator==, py::is_operator());
    clsCoaddBoundedField.def("__ne__", &CoaddBoundedField::operator!=, py::is_operator());
    clsCoaddBoundedField.def("__imul__", &CoaddBoundedField::operator*);

    /* Members */
    clsCoaddBoundedField.def("evaluate", &CoaddBoundedField::evaluate);
    clsCoaddBoundedField.def("isPersistable", &CoaddBoundedField::isPersistable);

    return mod.ptr();
}

}  // algorithms
}  // meas
}  // lsst
