// -*- LSST-C++ -*-
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
        try {
            algorithms::Centroid cen = cc->apply(10, 20);
        } catch(lsst::pex::exceptions::InvalidParameterException &e) {
            std::cerr << e << std::endl;
        }

        typename ImageT::Ptr image(new ImageT(100, 100));

        (*image)(10, 20) = 1000;

        cc->setImage(image);
        algorithms::Centroid cen = cc->apply(10, 20);

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
