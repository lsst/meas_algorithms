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

#ifndef LSST_MEAS_ALGORITHMS_DoubleGaussianPsf_h_INCLUDED
#define LSST_MEAS_ALGORITHMS_DoubleGaussianPsf_h_INCLUDED

#include "lsst/meas/algorithms/KernelPsf.h"

namespace lsst {
namespace meas {
namespace algorithms {

/// Represent a Psf as a circularly symmetrical double Gaussian
class DoubleGaussianPsf : public afw::table::io::PersistableFacade<DoubleGaussianPsf>, public KernelPsf {
public:
    /**
     *  Constructor for a DoubleGaussianPsf
     *
     *  @param[in] width    Number of columns in realisations of Psf
     *  @param[in] height   Number of rows in realisations of Psf
     *  @param[in] sigma1   Radius of inner Gaussian
     *  @param[in] sigma2   Radius of outer Gaussian
     *  @param[in] b        Ratio of Gaussian peak amplitudes: outer/inner
     */
    DoubleGaussianPsf(int width, int height, double sigma1, double sigma2 = 0.0, double b = 0.0);

    /// Polymorphic deep copy.  Usually unnecessary, as Psfs are immutable.
    PTR(afw::detection::Psf) clone() const override;

    /// Return a clone with specified kernel dimensions
    PTR(afw::detection::Psf) resized(int width, int height) const override;

    /// Return the radius of the inner Gaussian.
    double getSigma1() const { return _sigma1; }

    /// Return the radius of the outer Gaussian.
    double getSigma2() const { return _sigma2; }

    /// Return the ratio of Gaussian peak amplitudes: outer/inner
    double getB() const { return _b; }

    /// Whether this Psf is persistable (always true for DoubleGaussianPsf).
    bool isPersistable() const noexcept override { return true; }

protected:
    std::string getPersistenceName() const override;

    void write(OutputArchiveHandle& handle) const override;

private:
    double _sigma1;
    double _sigma2;
    double _b;
};

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst

#endif  // !LSST_MEAS_ALGORITHMS_DoubleGaussianPsf_h_INCLUDED
