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
 * Represent a PSF as a circularly symmetrical single Gaussian
 *
 * @file
 *
 * @ingroup algorithms
 */
#include <cmath>
#include "lsst/afw/math/FunctionLibrary.h"
#include "lsst/pex/exceptions.h"
#include "lsst/meas/algorithms/detail/sgPsf.h"
#include "lsst/afw/image/ImageUtils.h"

namespace afwMath = lsst::afw::math;

namespace lsst {
namespace meas {
namespace algorithms {

/************************************************************************************************************/
/**
 * Constructor for a sgPsf
 */
sgPsf::sgPsf(int width,                         ///< Number of columns in realisations of PSF
             int height,                        ///< Number of rows in realisations of PSF
             double sigma,                       ///< Width of Gaussian
             double,        ///< needed to match PSF.h definition but unused
             double         ///< needed to match PSF.h definition but unused
            ) :
    KernelPsf(), _sigma(sigma)
{
    if (sigma <= 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::DomainErrorException,
                          (boost::format("sigma may not be 0: %g") % sigma).str());
    }
    
    if (width > 0) {
        afwMath::GaussianFunction1<double> sg(sigma);
        setKernel(afwMath::SeparableKernel::Ptr(new afwMath::SeparableKernel(width, height, sg, sg)));
    }
}

//
// We need to make an instance here so as to register it
//
// \cond
namespace {
    volatile bool isInstance =
        lsst::afw::detection::Psf::registerMe<sgPsf,
                                              boost::tuple<int, int, double,double,double> >("SingleGaussian");
}

// \endcond
}}}
