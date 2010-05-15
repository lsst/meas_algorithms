// -*- LSST-C++ -*-
//
// test a perfect Gaussian PSF and measure aperture photometry at different radii
//
#include <iostream>
#include <cmath>

#include "lsst/afw.h"
#include "lsst/meas/algorithms/PSF.h"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE PsfAttributes

#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"

namespace measAlg = lsst::meas::algorithms;
namespace afwImage = lsst::afw::image;

BOOST_AUTO_TEST_CASE(PsfAttributes) {

    double sigma0 = 5.0;
    double aEff0 = 4.0*M_PI*sigma0*sigma0;
    
    int xwid = static_cast<int>(12*sigma0);
    int ywid = xwid;

    // set the peak of the outer guassian to 0 so this is really a single gaussian.
    measAlg::PSF::Ptr psf = measAlg::createPSF("DoubleGaussian", xwid, ywid, sigma0, sigma0, 0.0);

    measAlg::PsfAttributes psfAttrib(psf, xwid/2.0, ywid/2.0);
    double sigma = psfAttrib.computeGaussianWidth();
    double aEff = psfAttrib.computeEffectiveArea();
    
    std::cout << sigma0 << " " << sigma << std::endl;
    std::cout << aEff0 << " " << aEff << std::endl;
    
    BOOST_CHECK_CLOSE(sigma0, sigma, 1.0e-2);
    BOOST_CHECK_CLOSE(aEff0, aEff, 1.0e-2);

}
