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
 
#if !defined(LSST_DETECTION_PCAPSF_H)
#define LSST_DETECTION_PCAPSF_H
//!
// Describe an image's PSF
//
#include "lsst/meas/algorithms/PSF.h"

namespace lsst { namespace meas { namespace algorithms {
            
/*!
 * \brief Represent a PSF as a linear combination of PCA (== Karhunen-Loeve) basis functions
 */
class pcaPsf : public PSF {
public:
    typedef boost::shared_ptr<pcaPsf> Ptr;

    /**
     * @brief constructors for a pcaPsf
     *
     * Parameters:
     */
    explicit pcaPsf(lsst::afw::math::Kernel::Ptr kernel);

    lsst::afw::image::Image<PSF::Pixel>::Ptr getImage(double const x, double const y) const;
private:
    double doGetValue(double const dx, double const dy, int const xPositionInImage, int const yPositionInImage) const;

private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive&, unsigned int const) {
        boost::serialization::void_cast_register<pcaPsf, PSF>(
            static_cast<pcaPsf*>(0), static_cast<PSF*>(0));
    }
};

}}}

namespace boost {
namespace serialization {

template <class Archive>
inline void save_construct_data(
    Archive& ar, lsst::meas::algorithms::pcaPsf const* p,
    unsigned int const) {
    lsst::afw::math::Kernel const* kernel = p->getKernel().get();
    ar << make_nvp("kernel", kernel);
}

template <class Archive>
inline void load_construct_data(
    Archive& ar, lsst::meas::algorithms::pcaPsf* p,
    unsigned int const) {
    lsst::afw::math::Kernel* kernel;
    ar >> make_nvp("kernel", kernel);
    ::new(p) lsst::meas::algorithms::pcaPsf(
        lsst::afw::math::Kernel::Ptr(kernel));
}

}}

#endif
