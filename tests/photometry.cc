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
#include <limits>
#include <cmath>
#include "lsst/pex/policy/Policy.h"
#include "lsst/afw.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/image/ImageAlgorithm.h"
#include "lsst/meas/algorithms/Measure.h"
#include "lsst/meas/algorithms/Algorithm.h"
#include "lsst/meas/algorithms/FluxControl.h"
#include "lsst/afw/math/Integrate.h"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Photometry

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#include "boost/test/unit_test.hpp"
#pragma clang diagnostic pop
#include "boost/test/floating_point_comparison.hpp"

using namespace std;
namespace pexPolicy = lsst::pex::policy;
namespace measAlgorithms = lsst::meas::algorithms;
namespace afwDet = lsst::afw::detection;
namespace afwImage = lsst::afw::image;
namespace afwMath = lsst::afw::math;
namespace afwGeom = lsst::afw::geom;
namespace afwTable = lsst::afw::table;

typedef afwImage::Exposure<float, short unsigned int, float> ExposureT;
typedef ExposureT::MaskedImageT MImage;

/* =====================================================================
 * a functor for the PSF
 */
class Gaussian: public afwImage::pixelOp1XY<float> {
public:
    Gaussian(double const xcen, double const ycen, double const sigma, double const a) :
        _xcen(xcen), _ycen(ycen), _sigma(sigma), _a(a) {}
    float operator() (int const x, int const y, float) const {
        double const xx = x - _xcen;
        double const yy = y - _ycen;
        return _a * (1.0/(afwGeom::TWOPI*_sigma*_sigma)) *
            std::exp( -(xx*xx + yy*yy) / (2.0*_sigma*_sigma)  );
    }
private:
    double const _xcen;
    double const _ycen;
    double const _sigma;
    double const _a;
};


/* =====================================================================
 * a radial functor for the PSF
 */
class RGaussian: public std::unary_function<double, double> {
public:
    RGaussian(double const sigma, double const a, double const apradius, double const aptaper) :
        _sigma(sigma), _a(a), _apradius(apradius), _aptaper(aptaper) {}
    double operator() (double const r) const {
        double const gauss = _a * (1.0/(afwGeom::TWOPI*_sigma*_sigma)) *
            std::exp( -(r*r) / (2.0*_sigma*_sigma)  );
        double aperture;
        if ( r <= _apradius ) {
            aperture = 1.0;
        } else if ( r > _apradius && r < _apradius + _aptaper ) {
            aperture = 0.5*(1.0 + std::cos(afwGeom::PI*(r - _apradius)/_aptaper));
        } else {
            aperture = 0.0;
        }
        return aperture*gauss * (r * afwGeom::TWOPI);
    }
private:
    double const _sigma;
    double const _a;
    double const _apradius;
    double const _aptaper;
};



/**
 * This test performs a crude comparison between a Sinc-integrated aperture flux for a perfect Gaussian
 *   and the theoretical analytic flux integrated over the same Gaussian and aperture.
 * The Sinc method is expected to be in error by a small amount as the Gaussian psf is
 *   not band-limited (a requirement of the method)
 * The code is an abbreviation of the example "examples/growthcurve.cc" 
 */

BOOST_AUTO_TEST_CASE(PhotometrySinc) {

    // select the radii to test
    std::vector<double> radius;
    double r1 = 3.0;
    double r2 = 4.0;
    double dr = 1.0;
    int nR = static_cast<int>( (r2 - r1)/dr + 1 );
    for (int iR = 0; iR < nR; iR++) {
        radius.push_back(r1 + iR*dr);
    }

    double expectedError = 2.0;  // in percent
    
    // make an image big enough to hold the largest requested aperture
    int const xwidth = 2*(0 + 128);
    int const ywidth = xwidth;

    MImage mimg(afwGeom::ExtentI(xwidth, ywidth));
    ExposureT::Ptr exposure(new ExposureT(mimg));

    double const a = 100.0;
    double const aptaper = 2.0;
    float const xcen = xwidth/2;
    float const ycen = ywidth/2;
    afwGeom::Point2D center(xcen, ycen);
    //
    // The PSF widths that we'll test
    //
    std::vector<double> sigmas;
    sigmas.push_back(1.5);
    sigmas.push_back(2.5);
    int const nS = sigmas.size();

    measAlgorithms::SincFluxControl photomControl;

    for (int iS = 0; iS < nS; ++iS) {
        double const sigma = sigmas[iS];

        Gaussian gpsf(xcen, ycen, sigma, a);

        for_each_pixel(*mimg.getImage(), gpsf); // Set the image to a perfect Gaussian PSF

        double const psfH = 2.0*(r2 + 2.0);
        double const psfW = 2.0*(r2 + 2.0);
        
        afwDet::Psf::Ptr psf = afwDet::createPsf("DoubleGaussian", psfW, psfH, sigma);
        
        pexPolicy::Policy policy;
        for (int iR = 0; iR < nR; ++iR) {
            photomControl.radius2 = radius[iR];
            // Create the object that'll measure sinc aperture fluxes
            afwTable::Schema schema = afwTable::SourceTable::makeMinimalSchema();
            measAlgorithms::MeasureSources ms = 
                measAlgorithms::MeasureSourcesBuilder()
                .addAlgorithm(photomControl)
                .build(schema);
            PTR(afwTable::SourceTable) table = afwTable::SourceTable::make(schema);
            PTR(afwTable::SourceRecord) source = table->makeRecord();
            source->setFootprint(boost::make_shared<afwDet::Footprint>(exposure->getBBox()));
            ms.apply(*source, *exposure, center);
            afwTable::Flux::MeasKey key = table->getSchema()[photomControl.name];
            double const fluxSinc = source->get(key);
            // get the exact flux for the theoretical smooth PSF
            RGaussian rpsf(sigma, a, radius[iR], aptaper);
            double const fluxInt = afwMath::integrate(rpsf, 0, radius[iR], 1.0e-8);

            BOOST_CHECK_CLOSE(fluxSinc, fluxInt, expectedError);
        }
    }
}
    
