// -*- LSST-C++ -*-

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
 
/*!
 * @brief Represent a PSF as for a Coadd based on the James Jee stacking
 * algorithm which was extracted from Stackfit.
 *
 * Note that this Psf subclass only support computeImage, not the 
 * parameterization methodes defined on its super class.  In that sense,
 * it is not a true subclass.
 *
 * @file
 *
 * @ingroup algorithms
 */
#include <cmath>
#include <sstream>
#include <numeric>
#include "lsst/base.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/image/ImageUtils.h"
#include "lsst/afw/math/Statistics.h"
#include "lsst/meas/algorithms/detail/pcaPsf.h"
#include "lsst/meas/algorithms/CoaddPsf.h"
#include "lsst/coadd/utils/addToCoadd.h"

namespace afwDetection = lsst::afw::detection;
namespace coaddUtils = lsst::coadd::utils;
namespace afwMath = lsst::afw::math;
namespace afwGeom = lsst::afw::geom;
namespace afwImage = lsst::afw::image;
namespace lsst {
namespace meas {
namespace algorithms {

//
// Member Functions
//
/**
     *  @brief addComponent() used to set info about callexps
     *         as the information about individual calexps which went into
     *         the original Coadd.  This information is required to create     
     *         a similar stack of Psfs.
     *
     */

void ComponentVector::addComponent(PTR(lsst::afw::detection::Psf)  psf, PTR(lsst::afw::image::Wcs) wcs, lsst::afw::geom::Box2I bbox, double weight) {
    int size = _components.size();
    _components.resize(size + 1);
    _components[size].psf = psf;
    _components[size].wcs = wcs;
    _components[size].bbox = bbox;
    _components[size].weight = weight;
}

/**
     *  @brief set(ComponentVector) is used to set a vector of Components 
     *         as the information about individual calexps which went into
     *         the original Coadd.  This information is required to create     
     *         a similar stack of Psfs.
     *
     */

void ComponentVector::set(ComponentVector components) {
    _components.empty();
    for (int i = 0; i < components.size(); i++) {
        _components.push_back(components.at(i));
    }
}
 
    /**
     *  @brief  Forward a variety of different vector methods from
     *          the underlying component vector called '_components'
     */

int ComponentVector::size() const {
    return _components.size();
}

void ComponentVector::resize(int size) {
    _components.resize(size);
}

Component ComponentVector::ComponentVector::at(int i) const {
    return _components.at(i);
}

afwMath::Kernel::Ptr CoaddPsfKernel::clone() const {
    CoaddPsfKernel * clone = new CoaddPsfKernel();
    clone->setComponentVector(_components);
    afwMath::Kernel::Ptr retPtr(clone);
    return retPtr;
}

    /**
     *  @brief computeImage produces an estimate of the convolution Kernel at the given location
     *
     */

double CoaddPsfKernel::computeImage(afwImage::Image<double> &image, bool doNormalize, double x, double y) const {
    image *= 0.0;
    for (int i = 0; i < _components.size(); i++) {
        lsst::afw::geom::Box2I bbox = _components.at(i).bbox; 
        double xrel = x - bbox.getBeginX();
        double yrel = y - bbox.getBeginY();
        boost::shared_ptr<const lsst::meas::algorithms::PcaPsf> mypsf = boost::dynamic_pointer_cast<const lsst::meas::algorithms::PcaPsf>(_components.at(i).psf);
        afwGeom::Point2D point(xrel, yrel);
        PTR(afwImage::Image<double>) ii = mypsf->computeImage(point, true, true);
        image += *ii;
    }
   
   return 0;
}

    /**
     *  @brief setComponentVector is the only way to change the components in the underlying
     *         ComponentVector class.  A copy is made of the vector and the values and.
     *         shared pointers are copied into it
     */

void CoaddPsfKernel::setComponentVector(ComponentVector components) {

    _components.set(components);
}

int CoaddPsfKernel::getComponentCount() const {
    return _components.size();
}

/************************************************************************************************************/
/**
 * Constructors for a CoaddPsf
 */



CoaddPsf::CoaddPsf() : afwDetection::KernelPsf(boost::shared_ptr<lsst::afw::math::Kernel> (new lsst::meas::algorithms::CoaddPsfKernel()))
{
}

    /**
     *  @brief  null constructor with (CoaddPsfKernel argument) is supplied for class compatibility
     */

CoaddPsf::CoaddPsf(PTR(lsst::afw::math::Kernel) kernel ///< The desired Kernel
              ) : afwDetection::KernelPsf(kernel)
{
    
    // Check that it's a CoaddPsfKernel
    
    if (kernel.get() != NULL &&
        dynamic_cast<lsst::meas::algorithms::CoaddPsfKernel *>(kernel.get()) == NULL) {
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                          "CoaddPsf expects a CoaddPsfKernel");
    }
    setKernel(kernel);
}

CONST_PTR(lsst::meas::algorithms::CoaddPsfKernel) CoaddPsf::getCoaddPsfKernel() const {
   return boost::dynamic_pointer_cast<const lsst::meas::algorithms::CoaddPsfKernel>(getKernel());
}

//  Methods called from this class and related classes, all of which know this is a
//  CoaddPsf and must have a CoaddPsfKernel inside of it

PTR(lsst::meas::algorithms::CoaddPsfKernel) CoaddPsf::getCoaddPsfKernel() {
   return boost::dynamic_pointer_cast<lsst::meas::algorithms::CoaddPsfKernel>(getKernel());
}

//  The remaining routines are all forwarded to the CoaddPsfKernel class
double CoaddPsf::computeImage(afwImage::Image<double> &image, bool doNormalize, double x, double y) const {
    CoaddPsfKernel const * kernel = getCoaddPsfKernel().get();
    return kernel->computeImage(image, doNormalize, x, y);
}

void CoaddPsf::setComponentVector(ComponentVector components) {
    getCoaddPsfKernel().get()->setComponentVector(components);
}

int CoaddPsf::getComponentCount() const {
    return getCoaddPsfKernel().get()->getComponentCount();
}

//
// We need to make an instance here so as to register it with createPSF
//
// \cond
namespace {
    volatile bool isInstance =
        lsst::afw::detection::Psf::registerMe<CoaddPsf, PTR(lsst::afw::math::Kernel)>("COADD");
}

}}} // namespace lsst::meas::algorithms



// \endcond
