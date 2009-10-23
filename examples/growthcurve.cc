// -*- LSST-C++ -*-
//
// make a perfect PSF and measure aperture photometry at different radii
//
#include <iostream>
#include "lsst/afw.h"
#include "lsst/meas/algorithms/Photometry.h"
#include "lsst/afw/math/Integrate.h"

using namespace std;
namespace algorithms = lsst::meas::algorithms;
namespace image = lsst::afw::image;
namespace math = lsst::afw::math;

typedef image::MaskedImage<float, short unsigned int, float> MImage;

/* =====================================================================
 * a functor for the PSF
 */
class Gaussian: public std::binary_function<double, double, double> {
public:
    Gaussian(double const xcen, double const ycen, double const sigma, double const a) :
        _xcen(xcen), _ycen(ycen), _sigma(sigma), _a(a) {}
    double operator() (double const x, double const y) const {
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
    double pixOffset = 0.0;
    double errMult = 1.0;
    if (argc == 6) {
        r1 = atof(argv[1]);
        r2 = atof(argv[2]);
        dr = atof(argv[3]);
        pixOffset = atof(argv[4]);
        errMult = atof(argv[5]);
    }
    int nR = static_cast<int>( (r2 - r1)/dr + 1 );
    for (int iR = 0; iR < nR; iR++) {
        radius.push_back(r1 + iR*dr);
    }

    // make an image big enough to hold the largest requested aperture
    int const xwidth = 2*(0 + 128);
    int const ywidth = xwidth;

    int const nS = 2;
    std::vector<double> sigmas(2);
    sigmas[0] = 1.5;
    sigmas[1] = 2.5;
    double const a = 100.0;
    double const aptaper = errMult*2.0 + pixOffset;
    double const xcen = xwidth/2;
    double const ycen = ywidth/2;


    for (int iS = 0; iS < nS; ++iS) {

        double const sigma = sigmas[iS];

        Gaussian gpsf(xcen, ycen, sigma, a);

        // make a perfect Gaussian PSF in an image
        MImage const mimg(xwidth, ywidth);
        double xBcen = 0.0, yBcen = 0.0; // barycenters - crude centroids
        double fluxBarySum = 0.0;
        for (int iY = 0; iY != mimg.getHeight(); ++iY) {
            int iX = 0;
            for (MImage::x_iterator ptr = mimg.row_begin(iY), end = mimg.row_end(iY);
                 ptr != end; ++ptr, ++iX) {
                double const flux = gpsf(iX, iY);
                ptr.image() = flux;
                if (flux > 0.01) {
                    xBcen += flux*iX;
                    yBcen += flux*iY;
                    fluxBarySum += flux;
                }
            }
        }
        xBcen /= fluxBarySum;
        yBcen /= fluxBarySum;
        
        char outfits[20];
        sprintf(outfits, "fakestar_%3.1f.fits", sigma);
        mimg.getImage()->writeFits(outfits);
        
        for (int iR = 0; iR < nR; iR++) {

            double const psfH = 2.0*(r2 + 2.0);
            double const psfW = 2.0*(r2 + 2.0);
            
            // get the aperture flux
            algorithms::measurePhotometry<MImage> const *mp =
                algorithms::createMeasurePhotometry<MImage>("SINC", radius[iR]);
            algorithms::PSF::Ptr psf = algorithms::createPSF("DoubleGaussian", psfW, psfH, sigma);
            algorithms::Photometry phot = mp->apply(mimg, xcen, ycen, psf.get(), 0.0);
            double const fluxSinc = phot.getApFlux();
            double const fluxPsf = phot.getPsfFlux();
            
            // get the exact flux for the theoretical smooth PSF
            RGaussian rpsf(sigma, a, radius[iR], aptaper);
            double const fluxInt = math::integrate(rpsf, 0, radius[iR] + aptaper, 1.0e-8);
            cout << sigma << " " << radius[iR] << " " <<
                fluxInt << " " << fluxSinc << " " << fluxPsf << endl;
            //mimg.writeFits("mimg.fits");
        }
    }
}
    
