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
#include "lsst/cpputils/python.h"

#include "lsst/meas/algorithms/SpanSetMoments.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace algorithms {

void wrapSpanSetMoments(lsst::cpputils::python::WrapperCollection& wrappers) {
    using PySpanSetMoments = py::classh<SpanSetMoments>;
    auto clsSpanSetMoments =
            wrappers.wrapType(PySpanSetMoments(wrappers.module, "SpanSetMoments"), [](auto& mod, auto& cls) {
                cls.def_readonly("flux", &SpanSetMoments::flux);
                cls.def_readonly("variance", &SpanSetMoments::variance);
                // def_readonly doesn't ever seem to do the right thing w.r.t.
                // memory management, even with explicit return value policies.
                // Copying in a getter seems to be the only safe way to avoid
                // issues.
                cls.def_property_readonly(
                        "center", [](SpanSetMoments const& self) -> geom::Point2D { return self.center; });
                cls.def_property_readonly("shape",
                                          [](SpanSetMoments const& self) -> afw::geom::ellipses::Quadrupole {
                                              return self.shape;
                                          });
                cls.def_readonly("spans", &SpanSetMoments::spans);
                cls.def_readonly("too_many_bad_pixels", &SpanSetMoments::too_many_bad_pixels);
                cls.def_readonly("center_out_of_bounds", &SpanSetMoments::center_out_of_bounds);
                cls.def_readonly("bad_pixel_in_center", &SpanSetMoments::bad_pixel_in_center);
                cls.def_readonly("singular_second_moments", &SpanSetMoments::singular_second_moments);
                cls.def_property_readonly("any_flags_set", &SpanSetMoments::any_flags_set);
                cls.def("get_x_array", &SpanSetMoments::get_x_array);
                cls.def("get_y_array", &SpanSetMoments::get_y_array);
                cls.def_static("compute", &SpanSetMoments::compute, "spans"_a, "masked_image"_a,
                               "bad_bitmask"_a, "bad_pixel_max_fraction"_a, "bad_pixel_exclusion_radius"_a);
                cls.def_static("fit_shapelets", &SpanSetMoments::fit_shapelets, "masked_image"_a, "moments"_a,
                               "order"_a, "scale"_a, "circular"_a);
            });
}
}  // namespace algorithms
}  // namespace meas
}  // namespace lsst
