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
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/math/Kernel.h"

namespace afwImage = lsst::afw::image;
namespace afwGeom = lsst::afw::geom;
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
    void addComponent(PTR(lsst::afw::detection::Psf)  psf, PTR(lsst::afw::image::Wcs) wcs, lsst::afw::geom::Box2I bbox, double weight);
    void set(ComponentVector components);
    int size() const;
    
    void resize(int size);

    Component at(int i) const;

private:
    std::vector<Component> _components;

};

    /**
     *  @brief CoaddPsfKernel is the descendent of Kernel which is used to do computeImage for
     *         creating a convolution kernel when the corresponding image is a Coadd.
     *
     *  This kernel is really intended to be created privately by CoaddPsf.
     *  However, the class relationship of Psf and Kernel requires the constructor
     *  Psf(Kernel) to be available, so I have supported it.
     *
     *  Note that CoaddPsfKernel gets the info it needs to create a convolution image from
     *  a single call to setComponentVector.  It is expected that all the components
     *  will be set at the same time.
     */

class CoaddPsfKernel : public lsst::afw::math::Kernel {
public: 
    explicit CoaddPsfKernel() {
    };

    lsst::afw::math::Kernel::Ptr clone() const;

    friend class CoaddPsf;

private:

//  This is the critical piece of code to override.  The image needs to come from
//  the vector of PCA models, not from a spatially varying model itself

    double computeImage(
        afwImage::Image<double> &image,
        bool doNormalize,
        double x=0.0,
        double y=0.0
    ) const;

//  These methods are peculiar to the CoaddPsfKernel.  They are used to provide
//  the information about individual image psf's needed to supply a Psf at any point

    int getComponentCount() const;

    void setComponentVector(ComponentVector components);

//  Vector of components used to house the information from individual images

    ComponentVector _components;

};

    /**
     *  @brief CoaddPsf is the Psf descendent to be used for Coadd images.
     *  It incorporates the logic of James Jee's Stackfit algorithm
     *  for estimating the Psf of the stack of images (calexps)
     *  weighted by a given weighting vector
     *
     *  The user is expected to supply the ComponentVector which describe the 
     *  images which were added to the Coadd
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
    explicit CoaddPsf();
    explicit CoaddPsf(boost::shared_ptr<lsst::afw::math::Kernel>);
    virtual lsst::afw::detection::Psf::Ptr clone() const {
        return boost::make_shared<CoaddPsf>(*this); 
    }
    
    double computeImage(
        afwImage::Image<double> &image,
        bool doNormalize,
        double x=0.0,
        double y=0.0
    ) const;

    int getComponentCount() const;

    void setComponentVector(ComponentVector components);

    CONST_PTR(lsst::meas::algorithms::CoaddPsfKernel) getCoaddPsfKernel() const;
    PTR(lsst::meas::algorithms::CoaddPsfKernel) getCoaddPsfKernel();
};

}}}

//BOOST_CLASS_EXPORT_GUID(lsst::meas::algorithms::CoaddPsf, "lsst::meas::algorithms::coaddPsf") // lowercase initial for backward compatibility


#endif
