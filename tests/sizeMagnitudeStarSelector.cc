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

#include "lsst/meas/algorithms/SizeMagnitudeStarSelector.h"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SizeMagnitudeStarSelector

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#include "boost/test/unit_test.hpp"
#pragma clang diagnostic pop
#include "boost/test/floating_point_comparison.hpp"

namespace pexPolicy = lsst::pex::policy;
namespace measAlg = lsst::meas::algorithms;

BOOST_AUTO_TEST_CASE(SizeMagnitudeStarSelector) {

    pexPolicy::Policy pol;
    pol.set("minsize",0.0);
    pol.set("maxsize",1.0e100);
    pol.set("minmag",0.0);
    pol.set("maxmag",1.0e100);
    pol.set("logsize",false);
    pol.set("starfrac",0.5);
    pol.set("startn1",0.1);
    pol.set("fitorder",1);
    pol.set("fitsigclip",4.);
    pol.set("starsperbin",30);
    pol.set("purityratio",0.05);
    pol.set("aperture",5.0);
    std::cout<<"Pol = \n"<<pol<<std::endl;
    measAlg::SizeMagnitudeStarSelector selector(pol);
    std::cout<<"Successfully created selector\n";
}
