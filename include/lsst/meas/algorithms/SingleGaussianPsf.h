// -*- lsst-c++ -*-
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

#ifndef LSST_MEAS_ALGORITHMS_SingleGaussianPsf_h_INCLUDED
#define LSST_MEAS_ALGORITHMS_SingleGaussianPsf_h_INCLUDED

#include "lsst/base.h"
#include "lsst/meas/algorithms/KernelPsf.h"

namespace lsst {
namespace meas {
namespace algorithms {

/*!
 * @brief Represent a PSF as a circularly symmetrical Gaussian
 */
class SingleGaussianPsf : public afw::table::io::PersistableFacade<SingleGaussianPsf>, public KernelPsf {
public:
    /**
     *  @brief Constructor for a SingleGaussianPsf
     *
     *  @param[in] width   Number of columns in realizations of the PSF at a point.
     *  @param[in] height  Number of rows in realizations of the PSF at a point.
     *  @param[in] sigma   Radius of the Gaussian.
     *
     *  Additional arguments are historical and ignored, and maybe be removed in the future.
     */
    explicit SingleGaussianPsf(int width, int height, double sigma);

    /// Polymorphic deep copy; should usually unnecessary because Psfs are immutable.
    PTR(afw::detection::Psf) clone() const override;

    /// Return a clone with specified kernel dimensions
    PTR(afw::detection::Psf) resized(int width, int height) const override;

    /// Return the radius of the Gaussian.
    double getSigma() const { return _sigma; }

    /// Whether the Psf is persistable; always true.
    bool isPersistable() const noexcept override { return true; }

protected:
    std::string getPersistenceName() const override;

    void write(OutputArchiveHandle& handle) const override;

private:
    double _sigma;  ///< Width of Gaussian
};

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst

#endif  // !LSST_MEAS_ALGORITHMS_SingleGaussianPsf_h_INCLUDED
