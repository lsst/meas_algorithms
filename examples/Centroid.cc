// -*- LSST-C++ -*-
// Demonstrate how to measure centroids
//
#include <iostream>
#include "lsst/afw.h"
#include "lsst/meas/algorithms/Measure.h"

using namespace std;
namespace afwDetection = lsst::afw::detection;
namespace afwImage = lsst::afw::image;
namespace measAlgorithms = lsst::meas::algorithms;

typedef afwImage::Exposure<float> Exposure;

namespace {
    void getCentroid(std::string const& algorithm)
    {
        Exposure::Ptr exposure(new Exposure(100, 100));
        
        int const ix = 10;
        int const iy = 20;
        (*exposure->getMaskedImage().getImage())(ix, iy) = 1000;

        measAlgorithms::MeasureAstrometry<Exposure>::Ptr measureAstrom = measAlgorithms::makeMeasureAstrometry(exposure);
        measureAstrom->addAlgorithm(algorithm);
        
#if 1
        afwDetection::Measurement<afwDetection::Astrometry> photom =
            measureAstrom->measure(afwDetection::Peak(ix, iy));
        
        double const xcen = photom.find(algorithm)->getX();
        double const ycen = photom.find(algorithm)->getY();
#else
        double const xcen = measureAstrom->measure(afwDetection::Peak(ix, iy)).find(algorithm)->getX();
#endif
        
        cout << algorithm << ": (x, y) = " << xcen << ", " << ycen << endl;
    }
}

int main() {
    getCentroid("NAIVE");
    getCentroid("SDSS");
}
