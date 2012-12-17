// -*- lsst-c++ -*-
/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
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
//!
// Describe an image's PSF
//
#include "lsst/base.h"
#include "lsst/afw/detection/Psf.h"
#include "boost/serialization/nvp.hpp"
#include "boost/serialization/void_cast.hpp"

// Forward declarations

namespace lsst { namespace meas { namespace algorithms {
class SingleGaussianPsf;
}}}

namespace boost {
namespace serialization {
    template <class Archive>
    void save_construct_data(
        Archive& ar, lsst::meas::algorithms::SingleGaussianPsf const* p,
        unsigned int const file_version);
}}

namespace lsst { namespace meas { namespace algorithms {
            
/*!
 * @brief Represent a PSF as a circularly symmetrical double Gaussian
 */
class SingleGaussianPsf : public lsst::afw::table::io::PersistableFacade<SingleGaussianPsf>, 
                          public lsst::afw::detection::KernelPsf
{
public:
    typedef PTR(SingleGaussianPsf) Ptr;
    typedef CONST_PTR(SingleGaussianPsf) ConstPtr;

    /**
     * @brief constructors for a SingleGaussianPsf
     *
     * Parameters:
     */
    explicit SingleGaussianPsf(int width, int height, double sigma, double=0, double=0);

    virtual lsst::afw::detection::Psf::Ptr clone() const {
        return boost::make_shared<SingleGaussianPsf>(
            getKernel()->getWidth(), getKernel()->getHeight(),
            _sigma
        );
    }

    double getSigma() const { return _sigma; }

private:
    double _sigma;                     ///< Width of Gaussian

private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive&, unsigned int const) {
        boost::serialization::void_cast_register<SingleGaussianPsf, lsst::afw::detection::Psf>(
            static_cast<SingleGaussianPsf*>(0), static_cast<lsst::afw::detection::Psf*>(0));
    }
    template <class Archive>
    friend void boost::serialization::save_construct_data(
            Archive& ar, SingleGaussianPsf const* p, unsigned int const file_version);
};

}}} // namespace lsst::meas::algorithms

namespace boost { namespace serialization {

template <class Archive>
inline void save_construct_data(
    Archive& ar, lsst::meas::algorithms::SingleGaussianPsf const* p,
    unsigned int const) {
    int width = p->getKernel()->getWidth();
    int height = p->getKernel()->getHeight();
    ar << make_nvp("width", width);
    ar << make_nvp("height", height);
    ar << make_nvp("sigma", p->_sigma);
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
