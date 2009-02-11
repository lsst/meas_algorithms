//
// Demonstrate use of Centroider
//
#include <iostream>
#include "lsst/afw.h"
#include "lsst/meas/algorithms/Centroid.h"

using namespace std;
namespace algorithms = lsst::meas::algorithms;
namespace image = lsst::afw::image;

template<typename ImageT>
void computeCentroid(algorithms::Centroider<ImageT> const* cc) {
    ImageT image(100, 100);
    image = 0; image(10, 20) = 1000;
    
    algorithms::Centroid cen = cc->apply(image, 10, 20);

    cout << "(x, y) = " << cen.getX() << ", " << cen.getY() << endl;
}

int main() {
    typedef image::Image<float> ImageT;
    algorithms::Centroider<ImageT> *nc = algorithms::createCentroider<ImageT>("NAIVE");

    computeCentroid(nc);

    algorithms::Centroider<ImageT> *sdssc = algorithms::createCentroider<ImageT>("SDSS");
    computeCentroid(sdssc);
}
