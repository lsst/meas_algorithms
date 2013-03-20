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

#include "lsst/meas/algorithms/ImagePsf.h"

namespace lsst { namespace meas { namespace algorithms {

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
    explicit KernelPsf(
        afw::math::Kernel const & kernel,
        afw::geom::Point2D const & averagePosition=afw::geom::Point2D()
    );

    /// Return the Kernel used to define this Psf.
    PTR(afw::math::Kernel const) getKernel() const { return _kernel; }

    /// Return average position of stars; used as default position.
    virtual afw::geom::Point2D getAveragePosition() const;

    /// Polymorphic deep copy.
    virtual PTR(afw::detection::Psf) clone() const;

    /// Whether this object is persistable; just delegates to the kernel.
    virtual bool isPersistable() const;

protected:

    /// Construct a KernelPsf with the given kernel; it should not be modified afterwards.
    explicit KernelPsf(
        PTR(afw::math::Kernel) kernel,
        afw::geom::Point2D const & averagePosition=afw::geom::Point2D()
    );

    // Name to use persist this object as (should be overridden by derived classes).
    virtual std::string getPersistenceName() const;

    // Python module used in persistence to ensure factory is loaded.
    virtual std::string getPythonModule() const;

    // Output persistence implementation (should be overridden by derived classes if they add data members).
    virtual void write(OutputArchiveHandle & handle) const;

    // For access to protected ctor; avoids unnecessary copies when loading
    template <typename T, typename K> friend class KernelPsfFactory;

private:

    virtual PTR(Image) doComputeKernelImage(
        afw::geom::Point2D const & position,
        afw::image::Color const & color
    ) const;

    PTR(afw::math::Kernel) _kernel;
    afw::geom::Point2D _averagePosition;
};

}}} // namespace lsst::meas::algorithms

#endif // !LSST_MEAS_ALGORITHMS_KernelPsf_h_INCLUDED
