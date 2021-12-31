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

#ifndef LSST_MEAS_ALGORITHMS_PYTHON_H
#define LSST_MEAS_ALGORITHMS_PYTHON_H

#include "pybind11/pybind11.h"
#include "lsst/meas/algorithms/ImagePsf.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/detection/python.h"

using lsst::afw::detection::PsfTrampoline;
using lsst::afw::typehandling::StorableHelper;

namespace lsst {
namespace meas {
namespace algorithms {

/**
  * "Trampoline" for ImagePsf to let it be used as a base class in Python.
  *
  * Subclasses of ImagePsf that are wrapped in %pybind11 should have a similar
  * helper that subclasses `ImagePsfTrampoline<subclass>`. This helper can be
  * skipped if the subclass neither adds any virtual methods nor implements
  * any abstract methods.
  *
  * @tparam Base the exact (most specific) class being wrapped
  *
  * @see [pybind11 documentation](https://pybind11.readthedocs.io/en/stable/advanced/classes.html)
  */
template <typename Base = ImagePsf>
class ImagePsfTrampoline : public PsfTrampoline<Base> {
public:
    /**
     * Delegating constructor for wrapped class.
     *
     * While we would like to simply inherit base class constructors, when doing so, we cannot
     * change their access specifiers.  One consequence is that it's not possible to use inheritance
     * to expose a protected constructor to python.  The alternative, used here, is to create a new
     * public constructor that delegates to the base class public or protected constructor with the
     * same signature.
     *
     * @tparam Args  Variadic type specification
     * @param ...args  Arguments to forward to the Base class constructor.
     */
    template<typename... Args>
    ImagePsfTrampoline<Base>(Args... args) : PsfTrampoline<Base>(args...) {}

    double doComputeApertureFlux(
         double radius, geom::Point2D const& position,
         afw::image::Color const& color
    ) const override {
        PYBIND11_OVERLOAD_NAME(
            double, Base, "_doComputeApertureFlux", doComputeApertureFlux, radius, position, color
        );
    }

    afw::geom::ellipses::Quadrupole doComputeShape(
        geom::Point2D const& position,
        afw::image::Color const& color
    ) const override {
        PYBIND11_OVERLOAD_NAME(
            afw::geom::ellipses::Quadrupole, Base, "_doComputeShape", doComputeShape, position, color
        );
    }
};

}  // algorithms
}  // meas
}  // lsst




#endif // LSST_MEAS_ALGORITHMS_PYTHON_H
