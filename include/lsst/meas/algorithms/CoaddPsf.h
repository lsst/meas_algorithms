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
 
#if !defined(LSST_MEAS_ALGORITHMS_COADDPSF_H)
#define LSST_MEAS_ALGORITHMS_COADDPSF_H
//!
// Describe an image's PSF
//
#include <boost/make_shared.hpp>
#include "ndarray/eigen.h"
#include "lsst/base.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/detection/PsfFormatter.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/math/Kernel.h"

namespace lsst { namespace meas { namespace algorithms {

class CoaddPsfKernel : public lsst::afw::math::Kernel {
public:
 
    explicit CoaddPsfKernel() {};

    virtual lsst::afw::math::Kernel::Ptr clone() const;

    virtual double computeImage(
        lsst::afw::image::Image<lsst::afw::math::Kernel::Pixel> &image,
        bool doNormalize,
        double x=0.0,
        double y=0.0
    ) const;

    void addPsfComponent(PTR(lsst::afw::detection::Psf) psf, lsst::afw::geom::Box2D bbox, double weight);

private:
    struct Component {
        PTR(lsst::afw::detection::Psf) psf;
        lsst::afw::geom::Box2D bbox;
        double weight;
    };

    typedef std::vector<Component> ComponentVector;

    ComponentVector _components;
};

/*!
 * @brief Represent a PSF as a linear combination of PCA (== Karhunen-Loeve) basis functions
 */
class CoaddPsf : public lsst::afw::detection::KernelPsf {
public:
    typedef PTR(CoaddPsf) Ptr;
    typedef CONST_PTR(CoaddPsf) ConstPtr;

    /**
     * @brief constructors for a CoadPsf
     *
     * Parameters:
     */
    explicit CoaddPsf(){};
    explicit CoaddPsf(PTR(CoaddPsfKernel) kernel);
    explicit CoaddPsf(PTR(lsst::afw::math::Kernel) kernel);
    virtual lsst::afw::detection::Psf::Ptr clone() const {
        return boost::make_shared<CoaddPsf>(*this); 
    }
};

}}}

//BOOST_CLASS_EXPORT_GUID(lsst::meas::algorithms::CoaddPsf, "lsst::meas::algorithms::coaddPsf") // lowercase initial for backward compatibility


#endif
