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
#include "boost/serialization/nvp.hpp"
#include "boost/serialization/void_cast.hpp"

namespace lsst { namespace meas { namespace algorithms {

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
    virtual PTR(afw::detection::Psf) clone() const;

    /// Return the radius of the Gaussian.
    double getSigma() const { return _sigma; }

    /// Whether the Psf is persistable; always true.
    virtual bool isPersistable() const { return true; }

protected:

    virtual std::string getPersistenceName() const;

    virtual void write(OutputArchiveHandle & handle) const;

private:
    double _sigma;                     ///< Width of Gaussian

private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive&, unsigned int const) {
        boost::serialization::void_cast_register<SingleGaussianPsf, lsst::afw::detection::Psf>(
            static_cast<SingleGaussianPsf*>(0), static_cast<lsst::afw::detection::Psf*>(0));
    }
};

}}} // namespace lsst::meas::algorithms

namespace boost { namespace serialization {

template <class Archive>
inline void save_construct_data(
    Archive& ar, lsst::meas::algorithms::SingleGaussianPsf const* p,
    unsigned int const) {
    int width = p->getKernel()->getWidth();
    int height = p->getKernel()->getHeight();
    double sigma = p->getSigma();
    ar << make_nvp("width", width);
    ar << make_nvp("height", height);
    ar << make_nvp("sigma", sigma);
}

template <class Archive>
inline void load_construct_data(
    Archive& ar, lsst::meas::algorithms::SingleGaussianPsf* p,
    unsigned int const) {
    int width;
    int height;
    double sigma;
    ar >> make_nvp("width", width);
    ar >> make_nvp("height", height);
    ar >> make_nvp("sigma", sigma);
    ::new(p) lsst::meas::algorithms::SingleGaussianPsf(width, height, sigma);
}

}} // namespace boost::serialization

#endif // !LSST_MEAS_ALGORITHMS_SingleGaussianPsf_h_INCLUDED
