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
 
// Demonstrate how to measure centroids
//
#include <iostream>
#include "lsst/afw.h"
#include "lsst/meas/algorithms/Measure.h"

using namespace std;
namespace afwDetection = lsst::afw::detection;
namespace afwImage = lsst::afw::image;
namespace afwGeom = lsst::afw::geom;
namespace measAlgorithms = lsst::meas::algorithms;

typedef afwImage::Exposure<float> Exposure;

namespace {
    void getCentroid(std::string const& algorithm)
    {
        Exposure::Ptr exposure(new Exposure(afwGeom::ExtentI(100, 100)));
        
        int const ix = 10;
        int const iy = 20;
        (*exposure->getMaskedImage().getImage())(ix, iy) = 1000;

        measAlgorithms::MeasureAstrometry<Exposure> measureAstrom(*exposure, 
                                                                  lsst::pex::policy::Policy());
        measureAstrom.addAlgorithm(algorithm);

        afwDetection::Source source(0);
        afwDetection::Footprint::Ptr foot = boost::make_shared<afwDetection::Footprint>(exposure->getBBox());
        source.setFootprint(foot);
        afwDetection::Measurement<afwDetection::Astrometry>::Ptr meas = 
            measureAstrom.measure(source, exposure, afwGeom::Point2D(ix, iy));
        
        double const xcen = meas->find(algorithm)->getX();
        double const ycen = meas->find()->getY(); // you may omit "algorithm" if you specified exactly one
        
        cout << algorithm << ": (x, y) = " << xcen << ", " << ycen << endl;
    }
    //
    // Here's the same code, but written slightly differently
    //
    void getCentroid2(std::string const& algorithm)
    {
        Exposure::Ptr exposure(new Exposure(afwGeom::ExtentI(100, 100)));
        Exposure::MaskedImageT mi = exposure->getMaskedImage();
        mi.setXY0(100, 200);            // just to be nasty
        
        int const ix = mi.getX0() + 10;
        int const iy = mi.getY0() + 20;
        (*mi.getImage())(ix - mi.getX0(), iy - mi.getY0()) = 1000;

        lsst::pex::policy::Policy policy;
        policy.add(algorithm, lsst::pex::policy::Policy::Ptr(new lsst::pex::policy::Policy));

        afwDetection::Source source(0);
        afwDetection::Footprint::Ptr foot = boost::make_shared<afwDetection::Footprint>(exposure->getBBox());
        source.setFootprint(foot);
        afwGeom::Point2D center(ix, iy);

        afwDetection::Astrometry::Ptr centroid =
            measAlgorithms::MeasureAstrometry<Exposure>(*exposure, policy).measure(source, exposure,
                                                                                   center)->find();
        float const xcen = centroid->getX();
        float const ycen = centroid->getY();
        
        cout << algorithm << ": (x, y) = " << xcen << ", " << ycen << endl;
    }
}

int main() {
    getCentroid("NAIVE");
    getCentroid("SDSS");
    getCentroid2("GAUSSIAN");
}
