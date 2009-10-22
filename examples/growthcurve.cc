// -*- lsst-c++ -*-
//
// make a perfect PSF and measure aperture photometry at different radii
//
#include <iostream>
#include "lsst/afw.h"
#include "lsst/meas/algorithms/Photometry.h"
#include "lsst/afw/math/Integrate.h"
#include "../src/photometry/SincPhotometry.h"

using namespace std;
namespace algorithms = lsst::meas::algorithms;
namespace image = lsst::afw::image;
namespace math = lsst::afw::math;

typedef image::MaskedImage<float, short unsigned int, float> MImageT;

typedef double FTYPE;


/* =====================================================================
 * a functor for the PSF
 */
class Gaussian: public std::binary_function<FTYPE,FTYPE,FTYPE> {
public:
    Gaussian(FTYPE const xcen, FTYPE const ycen, FTYPE const sigma, FTYPE const A) :
        _xcen(xcen), _ycen(ycen), _sigma(sigma), _A(A) {}
    FTYPE operator() (FTYPE const x, FTYPE const y) const {
        FTYPE const xx = x - _xcen;
        FTYPE const yy = y - _ycen;
        return _A * (1.0/(2.0*M_PI*_sigma*_sigma)) *
            std::exp( -(xx*xx + yy*yy) / (2.0*_sigma*_sigma)  );
    }
private:
    FTYPE const _xcen;
    FTYPE const _ycen;
    FTYPE const _sigma;
    FTYPE const _A;
};


/* =====================================================================
 * a radial functor for the PSF
 */
class RGaussian: public std::unary_function<FTYPE,FTYPE> {
public:
    RGaussian(FTYPE const sigma, FTYPE const A, FTYPE const apradius, FTYPE const aptaper) :
        _sigma(sigma), _A(A), _apradius(apradius), _aptaper(aptaper) {}
    FTYPE operator() (FTYPE const r) const {
        FTYPE const gauss = _A * (1.0/(2.0*M_PI*_sigma*_sigma)) *
            std::exp( -(r*r) / (2.0*_sigma*_sigma)  );
        FTYPE aperture;
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
    FTYPE const _sigma;
    FTYPE const _A;
    FTYPE const _apradius;
    FTYPE const _aptaper;
};


/* =====================================================================
 *  MAIN
 */
int main(int argc, char *argv[]) {
    

    // select the radii to test
    std::vector<FTYPE> radius;
    FTYPE r1 = 3.0;
    FTYPE r2 = 3.0;
    FTYPE dr = 0.5;
    FTYPE pix_offset = 0.0;
    FTYPE err_mult = 1.0;
    if (argc == 6) {
        r1 = atof(argv[1]);
        r2 = atof(argv[2]);
        dr = atof(argv[3]);
        pix_offset = atof(argv[4]);
        err_mult = atof(argv[5]);
    }
    int n_r = (int) (r2 - r1)/dr + 1;
    for (int i_r = 0; i_r < n_r; i_r++) {
        radius.push_back(r1 + i_r*dr);
    }

    // make an image big enough to hold the largest requested aperture
    int const xwidth = 2*(0 + 128);
    int const ywidth = xwidth;

    int const n_s = 2;
    FTYPE const sigmas[n_s] = {1.5, 2.5};
    FTYPE const A = 100.0;
    FTYPE const aptaper = err_mult*2.0 + pix_offset;
    FTYPE const xcen = xwidth/2;
    FTYPE const ycen = ywidth/2;


    for (int i_s = 0; i_s < n_s; ++i_s) {

        FTYPE const sigma = sigmas[i_s];

        Gaussian gpsf(xcen, ycen, sigma, A);

        // make a perfect Gaussian PSF in an image
        MImageT const mimg(xwidth, ywidth);
        FTYPE x_bcen = 0.0, y_bcen = 0.0; // barycenters - crude centroids
        FTYPE flux_bary_sum = 0.0;
        for (int i_y = 0; i_y != mimg.getHeight(); ++i_y) {
            int i_x = 0;
            for (MImageT::x_iterator ptr = mimg.row_begin(i_y), end = mimg.row_end(i_y);
                 ptr != end; ++ptr, ++i_x) {
                FTYPE const flux = gpsf(i_x, i_y);
                ptr.image() = flux;
                if (flux > 0.01) {
                    x_bcen += flux*i_x;
                    y_bcen += flux*i_y;
                    flux_bary_sum += flux;
                }
            }
        }
        x_bcen /= flux_bary_sum;
        y_bcen /= flux_bary_sum;
        
        char outfits[20];
        sprintf(outfits, "fakestar_%3.1f.fits", sigma);
        mimg.getImage()->writeFits(outfits);
        
        for (int i_r = 0; i_r < n_r; i_r++) {
            

            // get the aperture flux
            algorithms::measurePhotometry<MImageT> const *mp =
                algorithms::createMeasurePhotometry<MImageT>("SINC",radius[i_r]);
            FTYPE const FWHM = 5.0;
            algorithms::PSF::Ptr psf =
                algorithms::createPSF("DoubleGaussian", 2*(r2+2), 2*(r2+2), FWHM/(2*sqrt(2*log(2))));
            algorithms::Photometry phot = mp->apply(mimg, xcen, ycen, &(*psf), 0.0);
            FTYPE flux00 = phot.getApFlux();
            FTYPE psfFlux00 = phot.getPsfFlux();
            
            // get the exact flux for the theoretical smooth PSF
            RGaussian rpsf(sigma, A, radius[i_r], aptaper);
            FTYPE const flux0 = math::integrate(rpsf, 0, radius[i_r] + aptaper, 1.0e-8);
            cout << sigma << " " << radius[i_r] << " " <<
                flux0 << " " << flux00 << " " << psfFlux00 << endl;
            //mimg.writeFits("mimg.fits");
        }
    }
}
    
