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
 
// Demonstrate use of MeasureCentroid
//
#include <iostream>
#include "lsst/afw.h"
#include "lsst/meas/algorithms/Centroid.h"

using namespace std;
namespace algorithms = lsst::meas::algorithms;
namespace afwImage = lsst::afw::image;

namespace {
    template<typename ImageT>
    void computeCentroid(algorithms::MeasureCentroid<ImageT> const* cc) {
        algorithms::Centroid cen;
        try {
            cen = cc->apply(10, 20);
        } catch(lsst::pex::exceptions::InvalidParameterException &e) {
            std::cerr << e << std::endl;
        }

        typename ImageT::Ptr image(new ImageT(100, 100));

        (*image)(10, 20) = 1000;

        cc->setImage(image);
        cen = cc->apply(10, 20);

        cen = cc->apply(*image, 10, 20);

        cout << "(x, y) = " << cen.getX() << ", " << cen.getY() << endl;
    }
}

int main() {
    typedef afwImage::Image<float> Image;
    Image::Ptr image(new Image(100, 100));
    (*image)(10, 20) = 1000;
    algorithms::MeasureCentroid<Image> *nc = algorithms::createMeasureCentroid<Image>("NAIVE", image);

    computeCentroid(nc);

    algorithms::MeasureCentroid<Image> *sdssc = algorithms::createMeasureCentroid<Image>("SDSS");
    computeCentroid(sdssc);
}
