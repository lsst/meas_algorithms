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
#include "boost/iterator/iterator_adaptor.hpp"
#include "boost/iterator/transform_iterator.hpp"
#include "lsst/base.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/image/ImageUtils.h"
#include "lsst/afw/math/Statistics.h"
#include "lsst/meas/algorithms/PcaPsf.h"
#include "lsst/meas/algorithms/CoaddPsf.h"
#include "lsst/afw/table/types.h"

namespace lsst {
namespace meas {
namespace algorithms {

//
// Member Functions
 /** *  @brief addComponent() used to set info about callexps
     *         as the information about individual calexps which went into
     *         the original Coadd.  This information is required to create     
     *         a similar stack of Psfs.
     *
     */

void ComponentVector::addComponent(lsst::afw::table::RecordId id, CONST_PTR(lsst::afw::detection::Psf)  psf, CONST_PTR(lsst::afw::image::Wcs) wcs, const lsst::afw::geom::Box2I bbox, double weight) {
    int size = _components.size();
    _components.resize(size + 1);
    _components[size].id = id;
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

/**
  * @brief CoaddPsf class
  *
  */

CoaddPsf::CoaddPsf(afw::table::ExposureCatalog const & catalog ) {
    setExposures(catalog);
}
    /**
     *  @brief computeImage produces an estimate of the convolution Kernel at the given location
     *
     */

double CoaddPsf::computeImage(afw::image::Image<double> &image, bool doNormalize, double x, double y) const {
    image *= 0.0;
    for (int i = 0; i < _components.size(); i++) {
        lsst::afw::geom::Box2I bbox = _components.at(i).bbox; 
        double xrel = x - bbox.getBeginX();
        double yrel = y - bbox.getBeginY();
        boost::shared_ptr<const lsst::meas::algorithms::PcaPsf> mypsf = boost::dynamic_pointer_cast<const lsst::meas::algorithms::PcaPsf>(_components.at(i).psf);
        afw::geom::Point2D point(xrel, yrel);
        PTR(afw::image::Image<double>) ii = mypsf->computeImage(point, true, true);
        image += *ii;
    }
   
   return 0;
}

    /**
     *  @brief setComponentVector is the only way to change the components in the underlying
     *         ComponentVector class.  A copy is made of the vector and the values and.
     *         shared pointers are copied into it
     */

void CoaddPsf::setComponentVector(ComponentVector components) {

    _components.set(components);
}

int CoaddPsf::getComponentCount() const {
    return _components.size();
}

void CoaddPsf::setExposures(afw::table::ExposureCatalog const & catalog)
{
    _components.resize(0);
    for (lsst::afw::table::ExposureCatalog::const_iterator i = catalog.begin(); i != catalog.end(); ++i) {
         lsst::afw::table::ExposureRecord const & r = *i;
         lsst::afw::table::RecordId id = r.getId();
         CONST_PTR(lsst::afw::detection::Psf) psf = r.getPsf();       
         CONST_PTR(lsst::afw::image::Wcs) wcs = r.getWcs();     
         lsst::afw::geom::Box2I bbox = r.getBBox(); 
         _components.addComponent(id, psf, wcs, bbox, 1.0);
    }
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
