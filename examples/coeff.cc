// -*- LSST-C++ -*-
// Compute a coefficient image for weighting of Sinc Aperture photometry
//
#include <iostream>
#include "lsst/afw.h"
#include "lsst/meas/algorithms/Photometry.h"
#include "lsst/meas/algorithms/detail/SincPhotometry.h"

using namespace std;
namespace algorithms = lsst::meas::algorithms;
namespace image = lsst::afw::image;

typedef image::MaskedImage<float> MaskedImage;

int main(int argc, char *argv[]) {
    
    double const xcen = 0;
    double const ycen = 0;
    double radius = 4.0;
    if (argc == 2) {
        radius = atof(argv[1]);
    }
    
    MaskedImage::ImagePtr cimage =
        algorithms::getCoeffImage<MaskedImage::Image::Pixel>(xcen, ycen, radius);
    cimage->writeFits("cimage.fits");
    
}
