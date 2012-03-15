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
#include "lsst/meas/algorithms/CentroidControl.h"

using namespace std;
namespace afwDetection = lsst::afw::detection;
namespace afwImage = lsst::afw::image;
namespace afwGeom = lsst::afw::geom;
namespace afwTable = lsst::afw::table;
namespace measAlgorithms = lsst::meas::algorithms;

typedef afwImage::Exposure<float> Exposure;

namespace {
    void getCentroid(measAlgorithms::AlgorithmControl const& ctrl)
    {
        Exposure::Ptr exposure(new Exposure(afwGeom::ExtentI(100, 100)));
        
        int const ix = 10;
        int const iy = 20;
        (*exposure->getMaskedImage().getImage())(ix, iy) = 1000;

        afwTable::Schema schema = afwTable::SourceTable::makeMinimalSchema();
        measAlgorithms::MeasureSources ms = measAlgorithms::MeasureSourcesBuilder()
            .addAlgorithm(ctrl)
            .build(schema);
        
        PTR(afwTable::SourceTable) table = afwTable::SourceTable::make(schema);
        PTR(afwTable::SourceRecord) source = table->makeRecord();
        afwDetection::Footprint::Ptr foot = boost::make_shared<afwDetection::Footprint>(exposure->getBBox());
        source->setFootprint(foot);

        ms.apply(*source, *exposure, afwGeom::Point2D(ix, iy));
        
        afwGeom::Point2D cen = source->get(schema.find< afwTable::Point<double> >(ctrl.name).key);
        
        cout << "(x, y) = " << cen.getX() << ", " << cen.getY() << endl;
    }
}

int main() {
    getCentroid(measAlgorithms::NaiveCentroidControl());
    getCentroid(measAlgorithms::SdssCentroidControl());
    getCentroid(measAlgorithms::GaussianCentroidControl());
}
