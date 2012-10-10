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
 * @brief Represent a PSF as a linear combination of PCA (== Karhunen-Loeve) basis functions
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

afwMath::Kernel::Ptr CoaddPsfKernel::clone() const {
    afwMath::Kernel::Ptr retPtr(new CoaddPsfKernel());
    return retPtr;
}

double CoaddPsfKernel::computeImage(afwImage::Image<double> &image, bool doNormalize, double x, double y) const {
//    psfstack = afwImg.ImageD(21,21,0.0)
//    for i in range(len(psfs)):
//        bbox = bboxes[i]
//        if x > bbox.getBeginX() and x < bbox.getEndX() and y > bbox.getBeginY() and y < bbox.getEndY():
//            xrel = x - bbox.getBeginX()
//            yrel = y - bbox.getBeginY()
//            #print "Box %d: (%.2f,%.2f) is in (%d,%d),(%d,%d) at (%.2f,%.2f)"%(i,x,y,bbox.getBeginX(),bbox.getBeginY(),bbox.getEndX(),bbox.getEndY(),xrel,yrel)
//            psf = psfs[i]
//            img = psf.computeImage(afwGeom.Point2D(xrel,yrel))
//            array = img.getArray()
//            sum_array = psfstack.getArray()
//            sum_array += weights[i] * array
//    return psfstack
//        PTR(lsst::afw::detection::Psf) psf;
//        PTR(lsst::afw::image::Wcs) wcs;
//        lsst::afw::geom::Box2I bbox;
//        double weight;
    std::cout << "starting computeImage \n";
    image *= 0.0;
    for (unsigned int i = 0; i < _components.size(); i++) {
        lsst::afw::geom::Box2I bbox = _components[i].bbox; 
        double xrel = x - bbox.getBeginX();
        double yrel = y - bbox.getBeginY();
        
        boost::shared_ptr<lsst::meas::algorithms::PcaPsf> mypsf = boost::dynamic_pointer_cast<lsst::meas::algorithms::PcaPsf>(_components[i].psf);
        afwGeom::Point2D point(xrel, yrel);
        std::cout << "computing image: \n";
        PTR(afwImage::Image<double>) ii = mypsf->computeImage(point, true, true);

        std::cout << "writing Fits \n";
        std::string s;
        std::stringstream out;
        out << i;
        s = out.str();
        ii->writeFits("kernel"+ s + ".fits");

        std::cout << "adding image \n";
        image += *ii;
    }
   
   return 0;
}

void CoaddPsfKernel::addPsfComponent(PTR(lsst::afw::detection::Psf) psf, PTR(lsst::afw::image::Wcs) wcs, lsst::afw::geom::Box2I bbox, double weight) {

    unsigned int size0 = _components.size();
    _components.resize(_components.size() + 1);
    unsigned int size1 = _components.size();
    std::cout << size0 << ", " << size1 << "\n";
    if (psf == NULL) std::cout << "psf is NULL\n";
    if (wcs == NULL) std::cout << "wcs is NULL\n";
    _components[size1-1].psf = psf;
    _components[size1-1].wcs = wcs;
    _components[size1-1].bbox = bbox;
    _components[size1-1].weight = weight;
    if (psf == NULL) 
        std::cout << " psf is null in add\n";
    else
    {
        std::cout << " computing image in add\n";
        afwGeom::Point2D point(0,0);
        psf->computeImage(point, true, true);
        std::cout << "done computing image in add\n";
    }
    if (wcs == NULL) 
        std::cout << "wcs is null in add\n";
}

int CoaddPsfKernel::getComponentCount() {
    return _components.size();
}

/************************************************************************************************************/
/**
 * Constructors for a CoaddPsf
 */

CoaddPsf::CoaddPsf() : afwDetection::KernelPsf(boost::shared_ptr<lsst::afw::math::Kernel> (new lsst::meas::algorithms::CoaddPsfKernel()))
{
}

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

CoaddPsf::CoaddPsf(PTR(lsst::meas::algorithms::CoaddPsfKernel) kernel ///< The desired Kernel
              ) : afwDetection::KernelPsf(boost::dynamic_pointer_cast<lsst::afw::math::Kernel>(kernel))
{
}

PTR(lsst::meas::algorithms::CoaddPsfKernel) CoaddPsf::getCoaddPsfKernel()
{
   return boost::dynamic_pointer_cast<lsst::meas::algorithms::CoaddPsfKernel>(getKernel());
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
