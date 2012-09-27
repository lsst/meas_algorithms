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
#include <numeric>
#include "lsst/base.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/image/ImageUtils.h"
#include "lsst/afw/math/Statistics.h"
#include "lsst/meas/algorithms/CoaddPsf.h"

namespace afwDetection = lsst::afw::detection;
namespace afwMath = lsst::afw::math;
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

double CoaddPsfKernel::computeImage(afwImage::Image<Pixel> &image, bool doNormalize, double x, double y) const {
   return 0;
}

//
void CoaddPsfKernel::addPsfComponent(PTR(lsst::afw::detection::Psf) psf, lsst::afw::geom::Box2D bbox, double weight) {
    _components.resize(_components.size() + 1);
    Component component = _components[_components.size()-1];
    component.psf = psf;
    component.bbox = bbox;
    component.weight = weight;
}
/************************************************************************************************************/
/**
 * Constructor for a CoaddPsf
 */
//#CoaddPsf::CoaddPsf(PTR(lsst::meas::algorithms::CoaddPsfKernel) kernel ///< The desired Kernel
CoaddPsf::CoaddPsf(PTR(lsst::meas::algorithms::CoaddPsfKernel) kernel ///< The desired Kernel
              ) : afwDetection::KernelPsf(kernel)
{
    //
    // Check that it's a LinearCombinationKernel
    //
    //if (kernel.get() != NULL &&
    //    dynamic_cast<lsst::meas::algorithms::CoaddPsfKernel *>(kernel.get()) == NULL) {
    //    throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
    //                      "CoaddPsf expects a LinearCombinationKernel");
    //}
}

CoaddPsf::CoaddPsf(PTR(lsst::afw::math::Kernel) kernel ///< The desired Kernel
              ) : afwDetection::KernelPsf(kernel)
{
    
    // Check that it's a LinearCombinationKernel
    
    if (kernel.get() != NULL &&
        dynamic_cast<lsst::meas::algorithms::CoaddPsfKernel *>(kernel.get()) == NULL) {
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                          "CoaddPsf expects a LinearCombinationKernel");
    }
}
// We need to make an instance here so as to register it with createPSF
//
// \cond
namespace {
    volatile bool isInstance =
        lsst::afw::detection::Psf::registerMe<CoaddPsf, PTR(lsst::afw::math::Kernel)>("COADD");
}

}}} // namespace lsst::meas::algorithms



// \endcond
