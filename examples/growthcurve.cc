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
// make a perfect PSF and measure aperture photometry at different radii
//
#include <iostream>
#include "lsst/afw.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/image/ImageAlgorithm.h"
#include "lsst/afw/math/Integrate.h"
#include "lsst/meas/algorithms/Measure.h"

using namespace std;
namespace pexPolicy = lsst::pex::policy;
namespace afwImage = lsst::afw::image;
namespace afwGeom = lsst::afw::geom;
namespace afwDetection = lsst::afw::detection;
namespace afwMath = lsst::afw::math;
namespace algorithms = lsst::meas::algorithms;

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
        return _a * (1.0/(2.0*M_PI*_sigma*_sigma)) *
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
        double const gauss = _a * (1.0/(2.0*M_PI*_sigma*_sigma)) *
            std::exp( -(r*r) / (2.0*_sigma*_sigma)  );
        double aperture;
        if ( r <= _apradius ) {
            aperture = 1.0;
        } else if ( r > _apradius && r < _apradius + _aptaper ) {
            aperture = 0.5*(1.0 + std::cos(M_PI*(r - _apradius)/_aptaper));
        } else {
            aperture = 0.0;
        }
        return aperture*gauss * (r * 2.0 * M_PI);
    }
private:
    double const _sigma;
    double const _a;
    double const _apradius;
    double const _aptaper;
};


/* =====================================================================
 *  MAIN
 */
int main(int argc, char *argv[]) {
    

    // select the radii to test
    std::vector<double> radius;
    double r1 = 3.0;
    double r2 = 3.0;
    double dr = 0.5;
    if (argc == 4) {
        r1 = atof(argv[1]);
        r2 = atof(argv[2]);
        dr = atof(argv[3]);
    }
    int nR = static_cast<int>( (r2 - r1)/dr + 1 );
    for (int iR = 0; iR < nR; iR++) {
        radius.push_back(r1 + iR*dr);
    }

    // make an image big enough to hold the largest requested aperture
    int const xwidth = 2*(0 + 128);
    int const ywidth = xwidth;

    MImage mimg(afwGeom::ExtentI(xwidth, ywidth));
    ExposureT::Ptr exposure(new ExposureT(mimg));

    std::vector<double> sigmas;
    sigmas.push_back(1.5);
    sigmas.push_back(2.5);

    double const a = 100.0;
    double const aptaper = 2.0;
    float const xcen = xwidth/2.0;
    float const ycen = ywidth/2.0;
    afwGeom::Point2D center(xcen, ycen);
    
    for (unsigned iS = 0; iS != sigmas.size(); ++iS) {
        double const sigma = sigmas[iS];

        Gaussian gpsf(xcen, ycen, sigma, a);

        // make a perfect Gaussian PSF in an image
        for_each_pixel(*mimg.getImage(), gpsf);
        //
        // Create the measuring object
        //
        algorithms::MeasurePhotometry<ExposureT> measurePhotom(*exposure, pexPolicy::Policy());
        measurePhotom.addAlgorithm("NAIVE");
        measurePhotom.addAlgorithm("PSF");
        measurePhotom.addAlgorithm("SINC");
        //
        // And the PSF
        //
        double const psfH = 2.0*(r2 + 2.0);
        double const psfW = 2.0*(r2 + 2.0);
        afwDetection::Psf::Ptr psf = afwDetection::createPsf("DoubleGaussian", psfW, psfH, sigma);
        exposure->setPsf(psf);

        pexPolicy::Policy policy = pexPolicy::Policy();
        
        for (int iR = 0; iR < nR; iR++) {
            policy.set("NAIVE.radius", radius[iR]);
            policy.set("SINC.radius",  radius[iR]);

            measurePhotom.configure(policy);

            afwDetection::Source source(0);
            source.setFootprint(boost::make_shared<afwDetection::Footprint>(exposure->getBBox()));
            afwDetection::Measurement<afwDetection::Photometry>::Ptr photom = 
                measurePhotom.measure(source, exposure, center);

            double const fluxNaive = photom->find("NAIVE")->getFlux();
            double const fluxPsf =   photom->find("PSF")->getFlux();
            double const fluxSinc =  photom->find("SINC")->getFlux();

            // get the exact flux for the theoretical smooth PSF
            RGaussian rpsf(sigma, a, radius[iR], aptaper);
            double const fluxInt = afwMath::integrate(rpsf, 0, radius[iR] + aptaper, 1.0e-8);

            // output
            cout << sigma << " " << radius[iR] << " " <<
                fluxInt << " " << fluxNaive << " " << fluxSinc << " " << fluxPsf << endl;

        }
    }
}
    
