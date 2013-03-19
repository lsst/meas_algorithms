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
 
#ifndef LSST_MEAS_ALGORITHMS_PcaPsf_h_INCLUDED
#define LSST_MEAS_ALGORITHMS_PcaPsf_h_INCLUDED

#include "lsst/meas/algorithms/KernelPsf.h"

namespace lsst { namespace meas { namespace algorithms {
            
/*!
 * @brief Represent a PSF as a linear combination of PCA (== Karhunen-Loeve) basis functions
 */
class PcaPsf : public lsst::afw::table::io::PersistableFacade<PcaPsf>, public KernelPsf {
public:
    typedef PTR(PcaPsf) Ptr;
    typedef CONST_PTR(PcaPsf) ConstPtr;

    /**
     *  @brief Constructor for a PcaPsf
     *
     *  @param[in] kernel           Kernel that defines the Psf.
     *  @param[in] averagePosition  Average position of stars used to construct the Psf.
     */
    explicit PcaPsf(
        PTR(afw::math::LinearCombinationKernel) kernel,
        afw::geom::Point2D const & averagePosition = afw::geom::Point2D()
    );

    /// Polymorphic deep copy.
    virtual PTR(afw::detection::Psf) clone() const;

    /// PcaPsf always has a LinearCombinationKernel, so we can override getKernel to make it more useful.
    PTR(afw::math::LinearCombinationKernel const) getKernel() const;

private:

    // Name used in table persistence; the rest of is implemented by KernelPsf.
    virtual std::string getPersistenceName() const { return "PcaPsf"; }

    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive&, unsigned int const) {
        boost::serialization::void_cast_register<PcaPsf,
            lsst::afw::detection::Psf>(static_cast<PcaPsf*>(0), static_cast<lsst::afw::detection::Psf*>(0));
    }

};

}}}

namespace boost {
namespace serialization {

template <class Archive>
inline void save_construct_data(
    Archive& ar, lsst::meas::algorithms::PcaPsf const* p,
    unsigned int const) {
    lsst::afw::math::LinearCombinationKernel const* kernel = p->getKernel().get();
    ar << make_nvp("kernel", kernel);
    lsst::afw::geom::Point2D averagePosition = p->getAveragePosition();
    ar << make_nvp("averagePositionX", averagePosition.getX());
    ar << make_nvp("averagePositionY", averagePosition.getY());
}

template <class Archive>
inline void load_construct_data(
    Archive& ar, lsst::meas::algorithms::PcaPsf* p,
    unsigned int const) {
    lsst::afw::math::LinearCombinationKernel* kernel;
    ar >> make_nvp("kernel", kernel);
    double x=0.0, y=0.0;
    ar >> make_nvp("averagePositionX", x);
    ar >> make_nvp("averagePositionY", y);
    ::new(p) lsst::meas::algorithms::PcaPsf(PTR(lsst::afw::math::LinearCombinationKernel)(kernel),
                                            lsst::afw::geom::Point2D(x, y));
}

}} // namespace boost::serialization

#endif // !LSST_MEAS_ALGORITHMS_PcaPsf_h_INCLUDED
