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
 
#if !defined(LSST_COADD_ALGORITHMS_DETAIL_PCAPSF_H)
#define LSST_COADD_ALGORITHMS_DETAIL_PCAPSF_H
//!
// Describe an image's PSF
//
#include <boost/make_shared.hpp>
#include "lsst/ndarray/eigen.h"
#include "lsst/base.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/detection/LocalPsf.h"
#include "lsst/afw/detection/PsfFormatter.h"

namespace lsst {
namespace afw {
    namespace math {
        class Kernel;
    }
}
namespace meas { namespace algorithms {
            
/*!
 * @brief Represent a PSF as a linear combination of PCA (== Karhunen-Loeve) basis functions
 */
class PcaPsf : public lsst::afw::detection::KernelPsf {
public:
    typedef PTR(PcaPsf) Ptr;
    typedef CONST_PTR(PcaPsf) ConstPtr;

    /**
     * @brief constructors for a PcaPsf
     *
     * Parameters:
     */
    explicit PcaPsf(PTR(lsst::afw::math::Kernel) kernel);
    virtual lsst::afw::detection::Psf::Ptr clone() const {
        return boost::make_shared<PcaPsf>(*this);
    } 
private:
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive&, unsigned int const) {
        boost::serialization::void_cast_register<PcaPsf,
            lsst::afw::detection::Psf>(static_cast<PcaPsf*>(0), static_cast<lsst::afw::detection::Psf*>(0));
    }

    virtual PTR(afw::detection::LocalPsf) doGetLocalPsf(
        lsst::afw::geom::Point2D const& point, 
        lsst::afw::image::Color const&
    ) const {
        //TODO make use of color
        Image::Ptr image = computeImage(point, false);
        
        double sum = image->getArray().asEigen().sum();
        if(sum != 0)
            *image /= sum;

        return boost::make_shared<afw::detection::ImageLocalPsf>(point, *image);
        
    };
};

}}}

BOOST_CLASS_EXPORT_GUID(lsst::meas::algorithms::PcaPsf, "lsst::meas::algorithms::pcaPsf") // lowercase initial for backward compatibility

lsst::daf::persistence::FormatterRegistration
lsst::afw::detection::PsfFormatter::pcaPsfRegistration("PcaPsf", typeid(lsst::meas::algorithms::PcaPsf),
                                                       lsst::afw::detection::PsfFormatter::createInstance);

namespace boost {
namespace serialization {

template <class Archive>
inline void save_construct_data(
    Archive& ar, lsst::meas::algorithms::PcaPsf const* p,
    unsigned int const) {
    lsst::afw::math::Kernel const* kernel = p->getKernel().get();
    ar << make_nvp("kernel", kernel);
}

template <class Archive>
inline void load_construct_data(
    Archive& ar, lsst::meas::algorithms::PcaPsf* p,
    unsigned int const) {
    lsst::afw::math::Kernel* kernel;
    ar >> make_nvp("kernel", kernel);
    ::new(p) lsst::meas::algorithms::PcaPsf(PTR(lsst::afw::math::Kernel)(kernel));
}

}}

#endif
