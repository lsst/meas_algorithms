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
 
//
// test a perfect Gaussian PSF and measure aperture photometry at different radii
//
#include <iostream>
#include <cmath>

#include "lsst/afw.h"
#include "lsst/meas/algorithms/SingleGaussianPsf.h"
#include "lsst/afw/geom/Angle.h"
#include "lsst/meas/algorithms/PSF.h"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE PsfAttributes

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#include "boost/test/unit_test.hpp"
#pragma clang diagnostic pop
#include "boost/test/floating_point_comparison.hpp"

namespace measAlg = lsst::meas::algorithms;
namespace afwDetection = lsst::afw::detection;
namespace afwImage = lsst::afw::image;
namespace afwGeom = lsst::afw::geom;

BOOST_AUTO_TEST_CASE(PsfAttributes) {

    double sigma0 = 5.0;
    double aEff0 = 4.0*afwGeom::PI*sigma0*sigma0;
    
    int xwid = static_cast<int>(12*sigma0);
    int ywid = xwid;

    std::shared_ptr<afwDetection::Psf> psf(new measAlg::SingleGaussianPsf(xwid, ywid, sigma0));

    measAlg::PsfAttributes psfAttrib(psf, xwid/2.0, ywid/2.0);
    double sigma = psfAttrib.computeGaussianWidth(measAlg::PsfAttributes::ADAPTIVE_MOMENT);
    double m1    = psfAttrib.computeGaussianWidth(measAlg::PsfAttributes::FIRST_MOMENT);
    double m2    = psfAttrib.computeGaussianWidth(measAlg::PsfAttributes::SECOND_MOMENT);
    double noise = psfAttrib.computeGaussianWidth(measAlg::PsfAttributes::NOISE_EQUIVALENT);
    double bick  = psfAttrib.computeGaussianWidth(measAlg::PsfAttributes::BICKERTON);
    double aEff  = psfAttrib.computeEffectiveArea();
    
    std::cout << sigma0 << " " << sigma << std::endl;
    std::cout << sigma0 << " " << m1 << std::endl;
    std::cout << sigma0 << " " << m2 << std::endl;
    std::cout << sigma0 << " " << noise << std::endl;
    std::cout << sigma0 << " " << bick << std::endl;
    std::cout << aEff0 << " " << aEff << std::endl;
    
    BOOST_CHECK_CLOSE(sigma0, sigma, 1.0e-2);
    BOOST_CHECK_CLOSE(sigma0, m1, 3.0e-2);
    BOOST_CHECK_CLOSE(sigma0, m2, 1.0e-2);
    BOOST_CHECK_CLOSE(sigma0, noise, 1.0e-2);
    BOOST_CHECK_CLOSE(sigma0, bick, 1.0e-2);
    BOOST_CHECK_CLOSE(aEff0, aEff, 1.0e-2);

}
