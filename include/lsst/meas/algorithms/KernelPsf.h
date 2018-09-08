// -*- LSST-C++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
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
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
#ifndef LSST_MEAS_ALGORITHMS_KernelPsf_h_INCLUDED
#define LSST_MEAS_ALGORITHMS_KernelPsf_h_INCLUDED

#include "lsst/geom/Box.h"
#include "lsst/meas/algorithms/ImagePsf.h"

namespace lsst {
namespace meas {
namespace algorithms {

/**
 *  @brief A Psf defined by a Kernel
 */
class KernelPsf : public afw::table::io::PersistableFacade<KernelPsf>, public ImagePsf {
public:
    /**
     *  @brief Construct a KernelPsf with a clone of the given kernel.
     *
     *  We clone the Kernel in the public constructor to ensure the Psf is immutable
     *  after construction (we don't want someone with another copy of the Kernel to
     *  be able to modify the one held by the Psf).
     *
     *  Derived classes may use the protected constructor, which takes a shared_ptr
     *  to Kernel and does not copy it.
     */
    explicit KernelPsf(afw::math::Kernel const& kernel,
                       geom::Point2D const& averagePosition = geom::Point2D());

    /// Return the Kernel used to define this Psf.
    PTR(afw::math::Kernel const) getKernel() const { return _kernel; }

    /// Return average position of stars; used as default position.
    geom::Point2D getAveragePosition() const override;

    /// Polymorphic deep copy.
    PTR(afw::detection::Psf) clone() const override;

    /// Return a clone with specified kernel dimensions
    PTR(afw::detection::Psf) resized(int width, int height) const override;

    /// Whether this object is persistable; just delegates to the kernel.
     bool isPersistable() const noexcept override;

protected:
    /// Construct a KernelPsf with the given kernel; it should not be modified afterwards.
    explicit KernelPsf(PTR(afw::math::Kernel) kernel, geom::Point2D const& averagePosition = geom::Point2D());

    // Name to use persist this object as (should be overridden by derived classes).
    std::string getPersistenceName() const override;

    // Python module used in persistence to ensure factory is loaded.
    std::string getPythonModule() const override;

    // Output persistence implementation (should be overridden by derived classes if they add data members).
    void write(OutputArchiveHandle& handle) const override;

    // For access to protected ctor; avoids unnecessary copies when loading
    template <typename T, typename K>
    friend class KernelPsfFactory;

private:
    PTR(Image)
    doComputeKernelImage(geom::Point2D const& position, afw::image::Color const& color) const override;

    geom::Box2I doComputeBBox(geom::Point2D const& position, afw::image::Color const& color) const override;

    PTR(afw::math::Kernel) _kernel;
    geom::Point2D _averagePosition;
};

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst

#endif  // !LSST_MEAS_ALGORITHMS_KernelPsf_h_INCLUDED
