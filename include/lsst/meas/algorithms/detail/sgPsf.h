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
 
#if !defined(LSST_COADD_ALGORITHMS_DETAIL_SGPSF_H)
#define LSST_COADD_ALGORITHMS_DETAIL_SGPSF_H
//!
// Describe an image's PSF
//
#include "lsst/base.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/math/shapelets.h"
#include "boost/serialization/nvp.hpp"
#include "boost/serialization/void_cast.hpp"

// Forward declarations

namespace lsst { namespace meas { namespace algorithms {
    class sgPsf;
}}}

namespace boost {
namespace serialization {
    template <class Archive>
    void save_construct_data(
        Archive& ar, lsst::meas::algorithms::sgPsf const* p,
        unsigned int const file_version);
}}

namespace lsst { namespace meas { namespace algorithms {
            
/*!
 * @brief Represent a PSF as a circularly symmetrical double Gaussian
 */
class sgPsf : public lsst::afw::detection::KernelPsf
{
public:
    typedef PTR(sgPsf) Ptr;
    typedef CONST_PTR(sgPsf) ConstPtr;

    /**
     * @brief constructors for a sgPsf
     *
     * Parameters:
     */
    explicit sgPsf(int width, int height, double sigma, double=0, double=0);
    virtual lsst::afw::detection::Psf::Ptr clone() const {
        return boost::make_shared<sgPsf>(
            getKernel()->getWidth(), getKernel()->getHeight(),
            _sigma
        );
    }
private:
    double _sigma;                     ///< Width of Gaussian

private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive&, unsigned int const) {
        boost::serialization::void_cast_register<sgPsf, lsst::afw::detection::Psf>(
            static_cast<sgPsf*>(0), static_cast<lsst::afw::detection::Psf*>(0));
    }
    template <class Archive>
    friend void boost::serialization::save_construct_data(
            Archive& ar, sgPsf const* p, unsigned int const file_version);
};

}}} // namespace lsst::meas::algorithms

namespace boost {
namespace serialization {

template <class Archive>
inline void save_construct_data(
    Archive& ar, lsst::meas::algorithms::sgPsf const* p,
    unsigned int const) {
    int width = p->getKernel()->getWidth();
    int height = p->getKernel()->getHeight();
    ar << make_nvp("width", width);
    ar << make_nvp("height", height);
    ar << make_nvp("sigma", p->_sigma);
}

template <class Archive>
inline void load_construct_data(
    Archive& ar, lsst::meas::algorithms::sgPsf* p,
    unsigned int const) {
    int width;
    int height;
    double sigma;
    ar >> make_nvp("width", width);
    ar >> make_nvp("height", height);
    ar >> make_nvp("sigma", sigma);
    ::new(p) lsst::meas::algorithms::sgPsf(width, height, sigma);
}

}}


#endif
