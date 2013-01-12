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
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/table/Exposure.h"
#include "lsst/afw/table/types.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/math/Kernel.h"

namespace lsst { namespace meas { namespace algorithms {

/**
 *  @brief CoaddPsf is the Psf descendent to be used for Coadd images.
 *  It incorporates the logic of James Jee's Stackfit algorithm
 *  for estimating the Psf of the stack of images (calexps)
 *  weighted by a given weighting vector
 *
 *  The user is expected to supply either an Exposure Catalog 
 *  which describes the images whose Psf's are to be stacked
 */

class CoaddPsf : public lsst::afw::detection::Psf {
public:
    typedef PTR(CoaddPsf) Ptr;
    typedef CONST_PTR(CoaddPsf) ConstPtr;

    /**
     * @brief constructors for a CoadPsf
     *
     * Parameters:
     */
    explicit CoaddPsf(afw::table::ExposureCatalog const & catalog) {
        setExposures(catalog);
    }; 

    explicit CoaddPsf(boost::shared_ptr<lsst::afw::math::Kernel>) {
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                          "CoaddPsf does not accept an lsst::afw::math::Kernel on its constructor");
    }


    virtual lsst::afw::detection::Psf::Ptr clone() const {
        return boost::make_shared<CoaddPsf>(*this); 
    }
    
    int getComponentCount() const;

    void setExposures(afw::table::ExposureCatalog const & catalog);


protected:
    lsst::afw::detection::Psf::Image::Ptr doComputeImage(lsst::afw::image::Color const& color,
                                  lsst::afw::geom::Point2D const& ccdXY,
                                  lsst::afw::geom::Extent2I const& size,
                                  bool normalizePeak,
                                  bool distort
                                 ) const; 

    lsst::afw::math::Kernel::Ptr doGetKernel(lsst::afw::image::Color const&) {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException, "CoaddPsf does not implement this method");
    }
        
    lsst::afw::math::Kernel::ConstPtr doGetKernel(lsst::afw::image::Color const&) const {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException, "CoaddPsf does not implement this method");
    }
        
    lsst::afw::math::Kernel::Ptr doGetLocalKernel(lsst::afw::geom::Point2D const&,
                                                          lsst::afw::image::Color const&) {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException, "CoaddPsf does not implement this method");
    }
        
    lsst::afw::math::Kernel::ConstPtr doGetLocalKernel(lsst::afw::geom::Point2D const&,
                                                               lsst::afw::image::Color const&) const {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException, "CoaddPsf does not implement this method");
    }

private:
    mutable lsst::afw::table::ExposureCatalog _catalog;
        
};

}}}


#endif
