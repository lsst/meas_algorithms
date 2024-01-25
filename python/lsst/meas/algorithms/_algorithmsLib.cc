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
#include "lsst/cpputils/python.h"

namespace py = pybind11;
using namespace pybind11::literals;
using lsst::cpputils::python::WrapperCollection;

namespace lsst {
namespace meas {
namespace algorithms {
void wrapCoaddTransmissionCurve(WrapperCollection &wrappers);
void wrapSingleGaussianPsf(WrapperCollection &wrappers);
void wrapKernelPsf(WrapperCollection &wrappers);
void wrapDoubleGaussianPsf(WrapperCollection &wrappers);
void wrapInterp(WrapperCollection &wrappers);
void wrapPcaPsf(WrapperCollection &wrappers);
void wrapWarpedPsf(WrapperCollection &wrappers);
void wrapPsfCandidate(WrapperCollection &wrappers);
void wrapImagePsf(WrapperCollection &wrappers);
void wrapSpatialModelPsf(WrapperCollection &wrappers);
void wrapCoaddPsf(WrapperCollection &wrappers);
void wrapCoaddBoundedField(WrapperCollection &wrappers);
void wrapCr(WrapperCollection &wrappers);
void wrapCloughTocher2DInterpolatorUtils(WrapperCollection &wrappers);

PYBIND11_MODULE(_algorithmsLib, mod) {
WrapperCollection wrappers(mod, "lsst.meas.algorithms");
    wrappers.addInheritanceDependency("lsst.afw.geom");
    wrappers.addInheritanceDependency("lsst.afw.image");
    wrappers.addInheritanceDependency("lsst.afw.table");
    wrappers.addInheritanceDependency("lsst.afw.detection");

    wrapImagePsf(wrappers);
    wrapKernelPsf(wrappers);
    wrapCoaddTransmissionCurve(wrappers);
    wrapSingleGaussianPsf(wrappers);
    wrapDoubleGaussianPsf(wrappers);
    wrapInterp(wrappers);
    wrapPcaPsf(wrappers);
    wrapWarpedPsf(wrappers);
    wrapPsfCandidate(wrappers);
    wrapSpatialModelPsf(wrappers);
    wrapCoaddPsf(wrappers);
    wrapCoaddBoundedField(wrappers);
    wrapCr(wrappers);
    wrapCloughTocher2DInterpolatorUtils(wrappers);
    wrappers.finish();
}

}  // algorithms
}  // meas
}  // lsst