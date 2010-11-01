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
 
// Compute a coefficient image for weighting of Sinc Aperture photometry
//
#include <iostream>
#include "lsst/afw.h"
#include "lsst/meas/algorithms/detail/SincPhotometry.h"

using namespace std;
namespace algorithms = lsst::meas::algorithms;
namespace image = lsst::afw::image;

typedef image::MaskedImage<float> MaskedImage;

int main(int argc, char *argv[]) {
    
    double radius = 4.0;
    if (argc == 2) {
        radius = atof(argv[1]);
    }
    
    MaskedImage::ImagePtr cimage =
        algorithms::detail::getCoeffImage<MaskedImage::Image::Pixel>(radius);
    cimage->writeFits("cimage.fits");
    
}
