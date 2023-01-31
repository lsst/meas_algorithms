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

#include "lsst/meas/algorithms/SpatialModelPsf.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace algorithms {
namespace {

template <typename PixelT>
void declareFunctions(lsst::cpputils::python::WrapperCollection &wrappers) {
    using MaskedImageT = afw::image::MaskedImage<PixelT, afw::image::MaskPixel, afw::image::VariancePixel>;
    auto &mod = wrappers.module;
    mod.def("createKernelFromPsfCandidates", createKernelFromPsfCandidates<PixelT>, "psfCells"_a, "dims"_a,
            "xy0"_a, "nEigenComponents"_a, "spatialOrder"_a, "ksize"_a, "nStarPerCell"_a = -1,
            "constantWeight"_a = true, "border"_a = 3);
    mod.def("countPsfCandidates", countPsfCandidates<PixelT>, "psfCells"_a, "nStarPerCell"_a = -1);
    mod.def("fitSpatialKernelFromPsfCandidates",
            (std::pair<bool, double>(*)(afw::math::Kernel *, afw::math::SpatialCellSet const &, int const,
                                        double const, double const))fitSpatialKernelFromPsfCandidates<PixelT>,
            "kernel"_a, "psfCells"_a, "nStarPerCell"_a = -1, "tolerance"_a = 1e-5, "lambda"_a = 0.0);
    mod.def("fitSpatialKernelFromPsfCandidates",
            (std::pair<bool, double>(*)(afw::math::Kernel *, afw::math::SpatialCellSet const &, bool const,
                                        int const, double const,
                                        double const))fitSpatialKernelFromPsfCandidates<PixelT>,
            "kernel"_a, "psfCells"_a, "doNonLinearFit"_a, "nStarPerCell"_a = -1, "tolerance"_a = 1e-5,
            "lambda"_a = 0.0);
    mod.def("subtractPsf", subtractPsf<MaskedImageT>, "psf"_a, "data"_a, "x"_a, "y"_a,
            "psfFlux"_a = std::numeric_limits<double>::quiet_NaN());
    mod.def("fitKernelParamsToImage", fitKernelParamsToImage<MaskedImageT>, "kernel"_a, "image"_a, "pos"_a);
    mod.def("fitKernelToImage", fitKernelToImage<MaskedImageT>, "kernel"_a, "image"_a, "pos"_a);
}
}  // namespace

void wrapSpatialModelPsf(lsst::cpputils::python::WrapperCollection &wrappers) {
    declareFunctions<float>(wrappers);
}

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst
