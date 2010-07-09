// -*- LSST-C++ -*-
//
// make a perfect PSF and measure aperture photometry at different radii
//
#include <iostream>
#include "lsst/afw.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/image/ImageAlgorithm.h"
#include "lsst/meas/algorithms/Photometry.h"
#include "lsst/afw/math/Integrate.h"

using namespace std;
namespace algorithms = lsst::meas::algorithms;
namespace afwImage = lsst::afw::image;
namespace afwDetection = lsst::afw::detection;
namespace afwMath = lsst::afw::math;

typedef afwImage::MaskedImage<float, short unsigned int, float> MImage;

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

    MImage::Ptr mimg(new MImage(xwidth, ywidth));

    std::vector<double> sigmas;
    sigmas.push_back(1.5);
    sigmas.push_back(2.5);

    double const a = 100.0;
    double const aptaper = 2.0;
    double const xcen = xwidth/2.0;
    double const ycen = ywidth/2.0;
    
    for (unsigned iS = 0; iS != sigmas.size(); ++iS) {
        double const sigma = sigmas[iS];

        Gaussian gpsf(xcen, ycen, sigma, a);

        // make a perfect Gaussian PSF in an image
        for_each_pixel(*mimg->getImage(), gpsf);
        //
        // Create the measuring objects
        //
        algorithms::MeasurePhotometry<MImage> const *mpSinc =
            algorithms::createMeasurePhotometry<MImage>("SINC", mimg);
        
        algorithms::MeasurePhotometry<MImage> const *mpNaive =
            algorithms::createMeasurePhotometry<MImage>("NAIVE", mimg);
        //
        // And the PSF
        //
        double const psfH = 2.0*(r2 + 2.0);
        double const psfW = 2.0*(r2 + 2.0);
        afwDetection::Psf::Ptr psf = afwDetection::createPsf("DoubleGaussian", psfW, psfH, sigma);
        
        for (int iR = 0; iR < nR; iR++) {
            mpNaive->setRadius(radius[iR]);
            mpSinc->setRadius(radius[iR]);

            // get the Naive aperture flux
            algorithms::Photometry photNaive = mpNaive->apply(xcen, ycen, psf.get(), 0.0);
            double const fluxNaive = photNaive.getApFlux();
            
            // get the Sinc aperture flux
            algorithms::Photometry photSinc = mpSinc->apply(xcen, ycen, psf.get(), 0.0);
            double const fluxSinc = photSinc.getApFlux();
            double const fluxPsf = photSinc.getPsfFlux();
            
            // get the exact flux for the theoretical smooth PSF
            RGaussian rpsf(sigma, a, radius[iR], aptaper);
            double const fluxInt = afwMath::integrate(rpsf, 0, radius[iR] + aptaper, 1.0e-8);

            // output
            cout << sigma << " " << radius[iR] << " " <<
                fluxInt << " " << fluxNaive << " " << fluxSinc << " " << fluxPsf << endl;

        }
    }
}
    
