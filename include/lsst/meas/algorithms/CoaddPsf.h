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
 *  @brief  Component is a structure used to pass the elements of the calexps needed
 *   to create a CoaddPsf
 *
 *  This structure contains information from the ISR and ImgChar tasks which is needed
 *  by the Psf stacking routine to find the optimal Psf for the stack of calexps.
 *  The bounding box is only an approximation of which pixels were contributed to by
 *  the calexp, in that it ignores masked pixels.
 */

struct Component {
    lsst::afw::table::RecordId id;
    CONST_PTR(lsst::afw::detection::Psf) psf;
    CONST_PTR(lsst::afw::image::Wcs) wcs;
    lsst::afw::geom::Box2I bbox;
    double weight;
}; 

/**
 *  @brief  ComponentVector is a vector-like collection of Components
 *
 *  Not all of the std::vector methods are supported, but they can be added as necessary
 *
 */

class ComponentVector {

public:

    explicit ComponentVector() {};
    
    Component operator[](int i) { return _components[i]; }

    void addComponent(lsst::afw::table::RecordId id, CONST_PTR(lsst::afw::detection::Psf)  psf, CONST_PTR(lsst::afw::image::Wcs) wcs, const lsst::afw::geom::Box2I bbox, double weight);

    void set(ComponentVector components);

    int size() const;
    
    void resize(int size);

    Component at(int i) const;

    friend class CoaddPsf;

private:

    std::vector<Component> _components;

};

/**
 *  @brief CoaddPsf is the Psf descendent to be used for Coadd images.
 *  It incorporates the logic of James Jee's Stackfit algorithm
 *  for estimating the Psf of the stack of images (calexps)
 *  weighted by a given weighting vector
 *
 *  The user is expected to supply either a ComponentVector or Exposure Catalog 
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
    explicit CoaddPsf(afw::table::ExposureCatalog const & catalog); 

    explicit CoaddPsf(boost::shared_ptr<lsst::afw::math::Kernel>) {
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                          "CoaddPsf does not accept an lsst::afw::math::Kernel on its constructor");
    }

    virtual lsst::afw::detection::Psf::Ptr clone() const {
        return boost::make_shared<CoaddPsf>(*this); 
    }
    
    double computeImage(
        afw::image::Image<double> &image,
        bool doNormalize,
        double x=0.0,
        double y=0.0
    ) const;

    int getComponentCount() const;

    void setComponentVector(ComponentVector components);

    void setExposures(afw::table::ExposureCatalog const & catalog);

    ComponentVector _components;

protected:

    Image::Ptr doComputeImage(lsst::afw::image::Color const& color,
                                      lsst::afw::geom::Point2D const& ccdXY,
                                      lsst::afw::geom::Extent2I const& size,
                                      bool normalizePeak,
                                      bool distort
                                     ) const {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException, "CoaddPsf does not implement this method");
    }

    lsst::afw::math::Kernel::Ptr doGetKernel(lsst::afw::image::Color const&) {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException, "CoaddPsf does not implement this method");
        return lsst::afw::math::Kernel::Ptr();
    }
        
    lsst::afw::math::Kernel::ConstPtr doGetKernel(lsst::afw::image::Color const&) const {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException, "CoaddPsf does not implement this method");
        return lsst::afw::math::Kernel::Ptr();
    }
        
    lsst::afw::math::Kernel::Ptr doGetLocalKernel(lsst::afw::geom::Point2D const&,
                                                          lsst::afw::image::Color const&) {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException, "CoaddPsf does not implement this method");
        return lsst::afw::math::Kernel::Ptr();
    }
        
    lsst::afw::math::Kernel::ConstPtr doGetLocalKernel(lsst::afw::geom::Point2D const&,
                                                               lsst::afw::image::Color const&) const {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException, "CoaddPsf does not implement this method");
        return lsst::afw::math::Kernel::Ptr();
    }

        
};

}}}

// pgee:  commented out, moved to cc file BOOST_CLASS_EXPORT_GUID(lsst::meas::algorithms::CoaddPsf, "lsst::meas::algorithms::coaddPsf") // lowercase initial for backward compatibility


#endif
