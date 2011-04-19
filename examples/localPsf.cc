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
#include <iostream>
#include <Eigen/Core>

#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/detection/FootprintArray.h"
#include "lsst/afw/detection/FootprintArray.cc"
#include "lsst/afw/detection/LocalPsf.h"
#include "lsst/ndarray.h"
#include "lsst/ndarray/eigen.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/math/shapelets.h"
#include "lsst/meas/algorithms/Shapelet.h"
#include "lsst/meas/algorithms/ShapeletKernel.h"

namespace ndarray = lsst::ndarray;
namespace measAlg = lsst::meas::algorithms;
namespace afwImage = lsst::afw::image;
namespace afwGeom = lsst::afw::geom;
namespace afwMath = lsst::afw::math;
namespace afwDet = lsst::afw::detection;

int main() {
    int order = 4;
    double sigma = 2.0;
    Eigen::VectorXd coeff = Eigen::VectorXd::Ones((order+1)*(order+2)/2); 

    measAlg::Shapelet::Ptr mjShapelet(
        new measAlg::Shapelet(order, sigma, coeff)
    );
    measAlg::LocalShapeletKernel mjKernel(
        mjShapelet, 
        afwImage::Wcs::Ptr(new afwImage::Wcs()),
        afwGeom::ExtentI(15)
    );
    afwImage::Image<double> mjImage(afwGeom::ExtentI(15));
    mjKernel.computeImage(mjImage, true);
    
    afwMath::shapelets::ShapeletFunction jbShapelet(
        mjShapelet->getOrder(), afwMath::shapelets::LAGUERRE,
        mjShapelet->getSigma(), afwGeom::Point2D(7,7)
    );
   
    jbShapelet.getCoefficients().deep() = ndarray::viewVectorAsArray(
        coeff
    );

    afwMath::shapelets::MultiShapeletFunction jbMsf(jbShapelet);
    afwDet::ShapeletLocalPsf jbLocalPsf(
        afwGeom::Point2D(7, 7), jbMsf, true
    );

    
    afwImage::Image<double> jbImage(mjImage.getBBox(afwImage::PARENT));
    afwDet::Footprint fp(mjImage.getBBox(afwImage::PARENT));
    ndarray::Array<double, 2,2> jbArray = 
        ndarray::static_dimension_cast<2>(jbImage.getArray());
    jbLocalPsf.evaluatePointSource(fp, ndarray::flatten<1>(jbArray));

    mjImage.writeFits("mjImage.fits");
    jbImage.writeFits("jbImage.fits");

}


