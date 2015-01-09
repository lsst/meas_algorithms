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
 
#include <numeric>
#include <cmath>
#include <functional>
#include <map>

#include <complex>
#include "boost/math/special_functions/bessel.hpp"
#include <fftw3.h>

#include "boost/shared_array.hpp"
#include "boost/lambda/lambda.hpp"
#include "boost/lambda/bind.hpp"
#include "boost/regex.hpp"
#include "boost/shared_array.hpp"

#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/image.h"
#include "lsst/afw/math/Integrate.h"
#include "lsst/meas/algorithms/Measure.h"
#include "lsst/meas/algorithms/Photometry.h"
#include "lsst/meas/algorithms/detail/SincPhotometry.h"
#include "lsst/afw/math/offsetImage.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/detection/FootprintFunctor.h"
#include "lsst/meas/algorithms/FluxControl.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace afwDet = lsst::afw::detection;
namespace afwImage = lsst::afw::image;
namespace afwMath = lsst::afw::math;
namespace afwGeom = lsst::afw::geom;

namespace lsst {
namespace meas {
namespace algorithms {
namespace {

// Convenient wrapper for a Bessel function
inline double J1(double const x) {
    return boost::math::cyl_bessel_j(1, x);
}

// sinc function
template<typename T>
inline T sinc(T const x) {
    return (x != 0.0) ? (std::sin(x) / x) : 1.0;
}

/******************************************************************************/
/**
 * @brief Define a circular aperture function object g_i, cos-tapered?
 */
template<typename CoordT>            
class CircularAperture {
public:
    
    CircularAperture(
        CoordT const radius1,    ///< inner radius of the aperture
        CoordT const radius2,    ///< outer radius of the aperture
        CoordT const taperwidth ///< width to cosine taper from 1.0 to 0.0 (ie. 0.5*cosine period)
    ):
        _radius1(radius1),
        _radius2(radius2),
        _taperwidth1(taperwidth),
        _taperwidth2(taperwidth),
        _k1(1.0/(2.0*taperwidth)),
        _k2(1.0/(2.0*taperwidth)),
        _taperLo1(radius1 - 0.5*taperwidth),
        _taperHi1(radius1 + 0.5*taperwidth),
        _taperLo2(radius2 - 0.5*taperwidth),
        _taperHi2(radius2 + 0.5*taperwidth) {

        // if we're asked for a radius smaller than our taperwidth,
        // adjust the taper width smaller so it fits exactly
        // with smooth derivative=0 at r=0

        if (_radius1 > _radius2) {
            throw LSST_EXCEPT(pexExceptions::InvalidParameterException,
                              (boost::format("rad2 less than rad1: (rad1=%.2f, rad2=%.2f) ") %
                               _radius1 % _radius2).str());
        }
        if (_radius1 < 0.0 || _radius2 < 0.0) {
            throw LSST_EXCEPT(pexExceptions::InvalidParameterException,
                              (boost::format("radii must be >= 0 (rad1=%.2f, rad2=%.2f) ") %
                               _radius1 % _radius2).str());
        }
        
        if (_radius1 == 0) {
            _taperwidth1 = 0.0;
            _k1 = 0.0;
        }
        
        // if we don't have room to taper at r=0
        if ( _radius1 < 0.5*_taperwidth1) {
            _taperwidth1 = 2.0*_radius1;
            _k1 = 1.0/(2.0*_taperwidth1);
        }
            
        // if we don't have room to taper between r1 and r2
        if ((_radius2 - _radius1) < 0.5*(_taperwidth1+_taperwidth2)) {
            
            // if we *really* don't have room ... taper1 by itself is too big
            // - set taper1,2 to be equal and split the r2-r1 range
            if ((_radius2 - _radius2) < 0.5*_taperwidth1) {
                _taperwidth1 = _taperwidth2 = 0.5*(_radius2 - _radius1);
                _k1 = _k2 = 1.0/(2.0*_taperwidth1);
                
                // if there's room for taper1, but not taper1 and 2
            } else {
                _taperwidth2 = _radius2 - _radius1 - _taperwidth1;
                _k2 = 1.0/(2.0*_taperwidth2);
            }                
                
            _taperLo1 = _radius1 - 0.5*_taperwidth1; 
            _taperHi1 = _radius1 + 0.5*_taperwidth1;
            _taperLo2 = _radius2 - 0.5*_taperwidth2; 
            _taperHi2 = _radius2 + 0.5*_taperwidth2;
        }
    }
    

    // When called, return the throughput at the requested x,y
    // todo: replace the sinusoid taper with a band-limited
    CoordT operator() (CoordT const x, CoordT const y) const {
        CoordT const xyrad = std::sqrt(x*x + y*y);
        if ( xyrad < _taperLo1 ) {
            return 0.0;
        } else if (xyrad >= _taperLo1 && xyrad <= _taperHi1 ) {
            return 0.5*(1.0 + std::cos(  (afwGeom::TWOPI*_k1)*(xyrad - _taperHi1)));
        } else if (xyrad > _taperHi1 && xyrad <= _taperLo2 ) {
            return 1.0;
        } else if (xyrad > _taperLo2 && xyrad <= _taperHi2 ) {
            return 0.5*(1.0 + std::cos(  (afwGeom::TWOPI*_k2)*(xyrad - _taperLo2)));
        } else {
            return 0.0;
        }
    }
    
    CoordT getRadius1() { return _radius1; }
    CoordT getRadius2() { return _radius2; }
    
private:
    CoordT _radius1, _radius2;
    CoordT _taperwidth1, _taperwidth2;
    CoordT _k1, _k2;      // the angular wavenumber corresponding to a cosine with wavelength 2*taperwidth
    CoordT _taperLo1, _taperHi1;
    CoordT _taperLo2, _taperHi2;
};


template<typename CoordT>            
class CircApPolar : public std::unary_function<CoordT, CoordT> {
public:
    CircApPolar(double radius, double taperwidth) : _ap(CircularAperture<CoordT>(0.0, radius, taperwidth)) {}
    CoordT operator() (double r) const {
        return r*_ap(r, 0.0);
    }
private:
    CircularAperture<CoordT> _ap;
};
    
/******************************************************************************/

/**
 * Define a Sinc functor to be integrated over for Sinc interpolation
 */
template<typename IntegrandT>
class SincAperture : public std::binary_function<IntegrandT, IntegrandT, IntegrandT> {
public:
    
    SincAperture(
                 CircularAperture<IntegrandT> const &ap,
                 int const ix,        // sinc center x
                 int const iy         // sinc center y
                )
        : _ap(ap), _ix(ix), _iy(iy) {}
    
    IntegrandT operator() (IntegrandT const x, IntegrandT const y) const {
        double const fourierConvention = afwGeom::PI;
        double const dx = fourierConvention*(x - _ix);
        double const dy = fourierConvention*(y - _iy);
        double const fx = sinc<double>(dx);
        double const fy = sinc<double>(dy);
        return (1.0 + _ap(x, y)*fx*fy);
    }
    
private: 
    CircularAperture<IntegrandT> const &_ap;
    double _ix, _iy;
};
    


/******************************************************************************/
    
template <typename MaskedImageT, typename WeightImageT>
class FootprintWeightFlux : public afwDet::FootprintFunctor<MaskedImageT> {
public:
    FootprintWeightFlux(
                        MaskedImageT const& mimage, ///< The image the source lives in
                        PTR(WeightImageT const) wimage    ///< The weight image
                       ) :
        afwDet::FootprintFunctor<MaskedImageT>(mimage),
        _wimage(wimage),
        _sum(0.0), _sumVar(0.0),
        _x0(wimage->getX0()), _y0(wimage->getY0()) {}
    
    /// @brief Reset everything for a new Footprint
    void reset() {}        
    void reset(afwDet::Footprint const& foot) {
        _sum = 0.0;
        _sumVar = 0.0;

        afwGeom::BoxI const& bbox(foot.getBBox());
        _x0 = bbox.getMinX();
        _y0 = bbox.getMinY();

        if (bbox.getDimensions() != _wimage->getDimensions()) {
            throw LSST_EXCEPT(pexExceptions::LengthErrorException,
                              (boost::format("Footprint at %d,%d -- %d,%d is wrong size "
                                             "for %d x %d weight image") %
                               bbox.getMinX() % bbox.getMinY() % bbox.getMaxX() % bbox.getMaxY() %
                               _wimage->getWidth() % _wimage->getHeight()).str());
        }
    }
    
    /// @brief method called for each pixel by apply()
    void operator()(typename MaskedImageT::xy_locator iloc, ///< locator pointing at the image pixel
                    int x,                                 ///< column-position of pixel
                    int y                                  ///< row-position of pixel
                   ) {
        typename MaskedImageT::Image::Pixel ival = iloc.image(0, 0);
        typename MaskedImageT::Image::Pixel vval = iloc.variance(0, 0);
        typename WeightImageT::Pixel wval = (*_wimage)(x - _x0, y - _y0);
        _sum    += wval*ival;
        _sumVar += wval*wval*vval;
    }

    /// Return the Footprint's flux
    double getSum() const { return _sum; }
    double getSumVar() const { return _sumVar; }
    
private:
    PTR(WeightImageT const) _wimage;        // The weight image
    double _sum;                                      // our desired sum
    double _sumVar;                                   // sum of the variance
    int _x0, _y0;                                     // the origin of the current Footprint
};




class GaussPowerFunctor : public std::binary_function<double, double, double> {
public:
    GaussPowerFunctor(double sigma) : _sigma(sigma) {}

    double operator()(double kx, double ky) const {
        double const k = ::sqrt(kx*kx + ky*ky);
        double const gauss = ::exp(-0.5*k*k*_sigma*_sigma);
        return gauss*gauss;
    }
private:
    double _sigma;
};

std::pair<double, double> computeGaussLeakage(double const sigma) {

    GaussPowerFunctor gaussPower(sigma);
    
    double lim = afwGeom::PI;

    // total power: integrate GaussPowerFunctor -inf<x<inf, -inf<y<inf (can be done analytically) 
    double powerInf = afwGeom::PI/(sigma*sigma);

    // true power: integrate GaussPowerFunctor -lim<x<lim, -lim<y<lim (must be done numerically) 
    double truePower = afwMath::integrate2d(gaussPower, -lim, lim, -lim, lim, 1.0e-8);
    double trueLeak = (powerInf - truePower)/powerInf;

    // estimated power: function is circular, but coords are cartesian
    // - true power does the actual integral numerically, but we can estimate it by integrating
    //   in polar coords over lim <= radius < infinity.  The integral is analytic.
    double estLeak = ::exp(-sigma*sigma*afwGeom::PI*afwGeom::PI)/powerInf;
    
    return std::pair<double, double>(trueLeak, estLeak);
    
}
    

} // end of anonymous namespace

    
/************************************************************************************************************/

namespace detail {


    template<typename PixelT>
    typename afwImage::Image<PixelT>::Ptr
    calcImageRealSpace(double const rad1, double const rad2, double const taperwidth) {
        
        PixelT initweight = 0.0; // initialize the coeff values

        int log2   = static_cast<int>(::ceil(::log10(2.0*rad2)/log10(2.0)));
        if (log2 < 3) { log2 = 3; }
        int hwid = pow(2, log2);
        int width  = 2*hwid - 1;

        int const xwidth     = width;
        int const ywidth     = width;

        int const x0 = -xwidth/2;
        int const y0 = -ywidth/2;
    
        // create an image to hold the coefficient image
        typename afwImage::Image<PixelT>::Ptr coeffImage =
            boost::make_shared<afwImage::Image<PixelT> >(afwGeom::ExtentI(xwidth, ywidth), initweight);
        coeffImage->setXY0(x0, y0);

        // create the aperture function object
        // determine the radius to use that makes 'radius' the effective radius of the aperture
        double tolerance = 1.0e-12;
        double dr = 1.0e-6;
        double err = 2.0*tolerance;
        double apEff = afwGeom::PI*rad2*rad2;
        double radIn = rad2;
        int maxIt = 20;
        int i = 0;
        while (err > tolerance && i < maxIt) {
            CircApPolar<double> apPolar1(radIn, taperwidth);
            CircApPolar<double> apPolar2(radIn+dr, taperwidth); 
            double a1 = afwGeom::TWOPI * afwMath::integrate(apPolar1, 0.0, radIn+taperwidth, tolerance);
            double a2 = afwGeom::TWOPI * afwMath::integrate(apPolar2, 0.0, radIn+dr+taperwidth, tolerance);
            double dadr = (a2 - a1)/dr;
            double radNew = radIn - (a1 - apEff)/dadr;
            err = (a1 - apEff)/apEff;
            radIn = radNew;
            i++;
        }
        CircularAperture<double> ap(rad1, rad2, taperwidth);

        
        /* ******************************* */
        // integrate over the aperture
        
        // the limits of the integration over the sinc aperture
        double const limit = rad2 + taperwidth;
        double const x1 = -limit;
        double const x2 =  limit;
        double const y1 = -limit;
        double const y2 =  limit;
        
        for (int iY = y0; iY != y0 + coeffImage->getHeight(); ++iY) {
            int iX = x0;
            typename afwImage::Image<PixelT>::x_iterator end = coeffImage->row_end(iY-y0);
            for (typename afwImage::Image<PixelT>::x_iterator ptr = coeffImage->row_begin(iY-y0); ptr != end; ++ptr) {
                
                // create a sinc function in the CircularAperture at our location
                SincAperture<double> sincAp(ap, iX, iY);
                
                // integrate the sinc
                PixelT integral = afwMath::integrate2d(sincAp, x1, x2, y1, y2, 1.0e-8);
                
                // we actually integrated function+1.0 and now must subtract the excess volume
                // - just force it to zero in the corners
                double const dx = iX;
                double const dy = iY;
                *ptr = (std::sqrt(dx*dx + dy*dy) < xwidth/2) ?
                    integral - (x2 - x1)*(y2 - y1) : 0.0;
                ++iX;
            }
        }


        double sum = 0.0;
        for (int iY = y0; iY != y0 + coeffImage->getHeight(); ++iY) {
            typename afwImage::Image<PixelT>::x_iterator end = coeffImage->row_end(iY-y0);
            for (typename afwImage::Image<PixelT>::x_iterator ptr = coeffImage->row_begin(iY-y0); ptr != end; ++ptr) {
                sum += *ptr;
            }
        }
        
#if 0                           // debugging
        coeffImage->writeFits("cimage.fits");
#endif

        return coeffImage;
    }



    class FftShifter {
    public:
        FftShifter(int xwid) : _xwid(xwid) {}
        int shift(int x) {
            if (x >= _xwid/2) {
                return x - _xwid/2;
            } else {
                return x + _xwid/2 + 1;
            }
        }
    private:
        int _xwid;
    };

    std::pair<double, double> rotate(double x, double y, double angle) {
        double c = ::cos(angle);
        double s = ::sin(angle);
        return std::pair<double, double>(x*c + y*s, -x*s + y*c);
    }
    
    /** todo
     * - try sub pixel shift if it doesn't break even symmetry
     * - put values directly in an Image
     * - precompute the plan
     */

    template<typename PixelT>
    typename afwImage::Image<PixelT>::Ptr calcImageKSpaceCplx(double const rad1, double const rad2,
                                                              double const posAng, double const ellipticity
                                                             ) {
        
        // we only need a half-width due to symmetry
        // make the hwid 2*rad2 so we have some buffer space and round up to the next power of 2
        int log2   = static_cast<int>(::ceil(::log10(2.0*rad2)/log10(2.0)));
        if (log2 < 3) { log2 = 3; }
        int hwid = pow(2, log2);
        int wid  = 2*hwid - 1;
        int xcen = wid/2, ycen = wid/2;
        FftShifter fftshift(wid);
        
        boost::shared_array<std::complex<double> > cimg(new std::complex<double>[wid*wid]);
        std::complex<double> *c = cimg.get();
        // fftplan args: nx, ny, *in, *out, direction, flags
        // - done in-situ if *in == *out
        fftw_plan plan = fftw_plan_dft_2d(wid, wid,
                                          reinterpret_cast<fftw_complex*>(c),
                                          reinterpret_cast<fftw_complex*>(c),
                                          FFTW_BACKWARD, FFTW_ESTIMATE);
        
        // compute the k-space values and put them in the cimg array
        double const twoPiRad1 = afwGeom::TWOPI*rad1;
        double const twoPiRad2 = afwGeom::TWOPI*rad2;
        double const scale = (1.0 - ellipticity);
        for (int iY = 0; iY < wid; ++iY) {
            int const fY = fftshift.shift(iY);
            double const ky = (static_cast<double>(iY) - ycen)/wid;
            
            for (int iX = 0; iX < wid; ++iX) {
                
                int const fX = fftshift.shift(iX);
                double const kx = static_cast<double>(iX - xcen)/wid;

                // rotate
                std::pair<double, double> coo = rotate(kx, ky, posAng);
                double kxr = coo.first;
                double kyr = coo.second;
                // rescale
                double const k = ::sqrt(kxr*kxr + scale*scale*kyr*kyr);
                
                double const airy1 = (rad1 > 0 ? rad1*J1(twoPiRad1*k) : 0.0)/k;
                double const airy2 = rad2*J1(twoPiRad2*k)/k;
                double const airy = airy2 - airy1;
                
                c[fY*wid + fX] = std::complex<double>(scale*airy, 0.0);
            }
        }
        c[0] = scale*afwGeom::PI*(rad2*rad2 - rad1*rad1);
        
        // perform the fft and clean up after ourselves
        fftw_execute(plan);
        fftw_destroy_plan(plan);
        
        // put the coefficients into an image
        typename afwImage::Image<PixelT>::Ptr coeffImage =
            boost::make_shared<afwImage::Image<PixelT> >(afwGeom::ExtentI(wid, wid), 0.0);
        
        for (int iY = 0; iY != coeffImage->getHeight(); ++iY) {
            int iX = 0;
            typename afwImage::Image<PixelT>::x_iterator end = coeffImage->row_end(iY);
            for (typename afwImage::Image<PixelT>::x_iterator ptr = coeffImage->row_begin(iY); ptr != end; ++ptr) {
                int fX = fftshift.shift(iX);
                int fY = fftshift.shift(iY);
                *ptr = static_cast<PixelT>(c[fY*wid + fX].real()/(wid*wid));
                iX++;
            }
        }
        
        // reset the origin to be the middle of the image
        coeffImage->setXY0(-wid/2, -wid/2);
        return coeffImage;
    }


    
    // I'm not sure why this doesn't work with DCT-I (REDFT00), the DCT should take advantage of symmetry
    // and be much faster than the DFT.  It runs but the numbers are slightly off ...
    // but I have to do it as real-to-halfcplx (R2HC) to get the correct numbers.

    template<typename PixelT>
    typename afwImage::Image<PixelT>::Ptr calcImageKSpaceReal(double const rad1, double const rad2) {
        
        // we only need a half-width due to symmertry
        // make the hwid 2*rad2 so we have some buffer space and round up to the next power of 2
        int log2   = static_cast<int>(::ceil(::log10(2.0*rad2)/log10(2.0)));
        if (log2 < 3) { log2 = 3; }
        int hwid = pow(2, log2);
        int wid  = 2*hwid - 1;
        int xcen = wid/2, ycen = wid/2;
        FftShifter fftshift(wid);
        
        boost::shared_array<double> cimg(new double[wid*wid]);
        double *c = cimg.get();
        // fftplan args: nx, ny, *in, *out, kindx, kindy, flags
        // - done in-situ if *in == *out
        fftw_plan plan = fftw_plan_r2r_2d(wid, wid, c, c, FFTW_R2HC, FFTW_R2HC, FFTW_ESTIMATE);
        
        // compute the k-space values and put them in the cimg array
        double const twoPiRad1 = afwGeom::TWOPI*rad1;
        double const twoPiRad2 = afwGeom::TWOPI*rad2;
        for (int iY = 0; iY < wid; ++iY) {
            
            int const fY = fftshift.shift(iY);
            double const ky = (static_cast<double>(iY) - ycen)/wid;
            
            for (int iX = 0; iX < wid; ++iX) {
                int const fX = fftshift.shift(iX);
                
                // emacs indent breaks if this isn't separte
                double const iXcen = static_cast<double>(iX - xcen);
                double const kx = iXcen/wid;

                double const k = ::sqrt(kx*kx + ky*ky);
                double const airy1 = (rad1 > 0 ? rad1*J1(twoPiRad1*k) : 0.0)/k;
                double const airy2 = rad2*J1(twoPiRad2*k)/k;
                double const airy = airy2 - airy1;
                c[fY*wid + fX] = airy;

            }
        }
        int fxy = fftshift.shift(wid/2);
        c[fxy*wid + fxy] = afwGeom::PI*(rad2*rad2 - rad1*rad1);
        
        // perform the fft and clean up after ourselves
        fftw_execute(plan);
        fftw_destroy_plan(plan);
        
        // put the coefficients into an image
        typename afwImage::Image<PixelT>::Ptr coeffImage =
            boost::make_shared<afwImage::Image<PixelT> >(afwGeom::ExtentI(wid, wid), 0.0);
        
        for (int iY = 0; iY != coeffImage->getHeight(); ++iY) {
            int iX = 0;
            typename afwImage::Image<PixelT>::x_iterator end = coeffImage->row_end(iY);
            for (typename afwImage::Image<PixelT>::x_iterator ptr = coeffImage->row_begin(iY); ptr != end; ++ptr) {
                
                // now need to reflect the quadrant we solved to the other three
                int fX = iX < hwid ? hwid - iX - 1 : iX - hwid + 1;
                int fY = iY < hwid ? hwid - iY - 1 : iY - hwid + 1;
                *ptr = static_cast<PixelT>(c[fY*wid + fX]/(wid*wid));
                iX++;
            }
        }
    
        // reset the origin to be the middle of the image
        coeffImage->setXY0(-wid/2, -wid/2);
        return coeffImage;
    }
    
} // end of 'detail' namespace


namespace photometry {

template<typename PixelT>
SincCoeffs<PixelT>& SincCoeffs<PixelT>::getInstance()
{
    static SincCoeffs<PixelT> instance;
    return instance;
}

template<typename PixelT>
void SincCoeffs<PixelT>::cache(float r1, float r2)
{
    if (r1 < 0.0 || r2 < r1) {
        throw LSST_EXCEPT(pexExceptions::InvalidParameterException,
                          (boost::format("Invalid r1,r2 = %f,%f") % r1 % r2).str());
    }
    double const innerFactor = r1/r2;
    afw::geom::ellipses::Axes axes(r2, r2, 0.0);
    if (!getInstance()._lookup(axes, innerFactor)) {
        PTR(typename SincCoeffs<PixelT>::CoeffT) coeff = calculate(axes, innerFactor);
        coeff->markPersistent();
        getInstance()._cache[r2][innerFactor] = coeff;
    }
}

template<typename PixelT>
CONST_PTR(typename SincCoeffs<PixelT>::CoeffT)
SincCoeffs<PixelT>::get(afw::geom::ellipses::Axes const& axes, float const innerFactor)
{
    CONST_PTR(CoeffT) coeff = getInstance()._lookup(axes, innerFactor);
    return coeff ? coeff : calculate(axes, innerFactor);
}

template<typename PixelT>
CONST_PTR(typename SincCoeffs<PixelT>::CoeffT)
SincCoeffs<PixelT>::_lookup(afw::geom::ellipses::Axes const& axes, double const innerFactor) const
{
    if (innerFactor < 0.0 || innerFactor > 1.0) {
        throw LSST_EXCEPT(pexExceptions::InvalidParameterException,
                          (boost::format("innerFactor = %f is not between 0 and 1") % innerFactor).str());
    }

    CONST_PTR(typename SincCoeffs<PixelT>::CoeffT) const null = CONST_PTR(SincCoeffs<PixelT>::CoeffT)();

    // We only cache circular apertures
    if (!fuzzyCompare<float>().isEqual(axes.getA(), axes.getB())) {
        return null;
    }
    typename CoeffMapMap::const_iterator iter1 = _cache.find(axes.getA());
    if (iter1 == _cache.end()) {
        return null;
    }
    typename CoeffMap::const_iterator iter2 = iter1->second.find(innerFactor);
    return (iter2 == iter1->second.end()) ? null : iter2->second;
}

template<typename PixelT>
PTR(typename SincCoeffs<PixelT>::CoeffT)
SincCoeffs<PixelT>::calculate(afw::geom::ellipses::Axes const& axes, double const innerFactor)
{
    if (innerFactor < 0.0 || innerFactor > 1.0) {
        throw LSST_EXCEPT(pexExceptions::InvalidParameterException,
                          (boost::format("innerFactor = %f is not between 0 and 1") % innerFactor).str());
    }

    // Kspace-real is fastest, but only slightly faster than kspace cplx
    // but real won't work for elliptical apertures due to symmetries assumed for real transform

    double const rad1 = axes.getA() * innerFactor;
    double const rad2 = axes.getA();
    // if there's no angle and no ellipticity ... cache it/see if we have it cached
    if (fuzzyCompare<float>().isEqual(axes.getA(), axes.getB())) {
        // here we call the real transform
        return detail::calcImageKSpaceReal<PixelT>(rad1, rad2);
    } else {
        // here we call the complex transform
        double const ellipticity = 1.0 - axes.getB()/axes.getA();
        return detail::calcImageKSpaceCplx<PixelT>(rad1, rad2, axes.getTheta(), ellipticity);
    }
}


/************************************************************************************************************/
/**
 * Workhorse routine to calculate elliptical aperture fluxes
 */
template<typename MaskedImageT>
std::pair<double, double>
calculateSincApertureFlux(MaskedImageT const& mimage, afw::geom::ellipses::Ellipse const& ellipse,
                          double const innerFactor)
{
    double flux = std::numeric_limits<double>::quiet_NaN();
    double fluxErr = std::numeric_limits<double>::quiet_NaN();
    
    typedef typename MaskedImageT::Image Image;
    typedef typename Image::Pixel Pixel;
    typedef typename Image::Ptr ImagePtr;
    
    // BBox for data image
    afwGeom::BoxI imageBBox(mimage.getBBox(afwImage::PARENT));


    // make the coeff image
    // compute c_i as double integral over aperture def g_i(), and sinc()
    CONST_PTR(Image) cimage0 = SincCoeffs<Pixel>::get(ellipse.getCore(), innerFactor);
        
    // as long as we're asked for the same radius, we don't have to recompute cimage0
    // shift to center the aperture on the object being measured
    ImagePtr cimage = afwMath::offsetImage(*cimage0, ellipse.getCenter().getX(), ellipse.getCenter().getY());
    afwGeom::BoxI bbox(cimage->getBBox(afwImage::PARENT));
#if 0
    // I (Steve Bickerton) think this should work, but doesn't.
    // For the time being, I'll do the bounds check here
    // ... should determine why bbox/image behaviour not as expected.
    afwGeom::BoxI mbbox(mimage.getBBox(afwImage::PARENT));
    bbox.clip(mbbox);
    afwGeom::Point2I cimXy0(cimage->getXY0());
    bbox.shift(-cimage->getX0(), -cimage->getY0());
    cimage = typename Image::Ptr(new Image(*cimage, bbox, false));
    cimage->setXY0(cimXy0);
#else
    int x1 = (cimage->getX0() < mimage.getX0()) ? mimage.getX0() : cimage->getX0();
    int y1 = (cimage->getY0() < mimage.getY0()) ? mimage.getY0() : cimage->getY0();
    int x2 = (cimage->getX0() + cimage->getWidth() > mimage.getX0() + mimage.getWidth()) ?
        mimage.getX0() + mimage.getWidth() - 1 : cimage->getX0() + cimage->getWidth() - 1;
    int y2 = (cimage->getY0() + cimage->getHeight() > mimage.getY0() + mimage.getHeight()) ?
        mimage.getY0() + mimage.getHeight() - 1 : cimage->getY0() + cimage->getHeight() - 1; 
    
    // if the dimensions changed, put the image in a smaller bbox
    if ( (x2 - x1 + 1 != cimage->getWidth()) || (y2 - y1 + 1 != cimage->getHeight()) ) {
        bbox = afwGeom::BoxI(afwGeom::Point2I(x1 - cimage->getX0(), y1 - cimage->getY0()),
                             afwGeom::Extent2I(x2 - x1 + 1, y2 - y1 + 1));
        cimage = ImagePtr(new Image(*cimage, bbox, afwImage::LOCAL, false));
        
        // shift back to correct place
        cimage = afwMath::offsetImage(*cimage, x1, y1);
        bbox = afwGeom::BoxI(afwGeom::Point2I(x1, y1), 
                             afwGeom::Extent2I(x2-x1+1, y2-y1+1));
    }
#endif
        
    // pass the image and cimage into the wfluxFunctor to do the sum
    FootprintWeightFlux<MaskedImageT, Image> wfluxFunctor(mimage, cimage);
    
    afwDet::Footprint foot(bbox, imageBBox);
    wfluxFunctor.apply(foot);
    flux = wfluxFunctor.getSum();
    fluxErr = ::sqrt(wfluxFunctor.getSumVar());

    return std::make_pair(flux, fluxErr);
}
}
    
/************************************************************************************************************/

/**
 * @brief A class that knows how to calculate fluxes using the SINC photometry algorithm
 * @ingroup meas/algorithms
 */
class SincFlux : public FluxAlgorithm {
public:

    SincFlux(SincFluxControl const & ctrl, afw::table::Schema & schema) :
        FluxAlgorithm(ctrl, schema, "elliptical aperture photometry using sinc interpolation")
    {
        // calculate the needed coefficients 
        if (photometry::fuzzyCompare<float>().isEqual(ctrl.ellipticity, 0.0)) {
            photometry::SincCoeffs<float>::cache(ctrl.radius1, ctrl.radius2);
        }
    }

private:
    
    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(SincFlux);

};

/**
 * Calculate the desired aperture flux using the sinc algorithm
 */

template <typename PixelT>
void SincFlux::_apply(
    afw::table::SourceRecord & source, 
    afw::image::Exposure<PixelT> const& exposure,
    afw::geom::Point2D const & center
) const {
    source.set(getKeys().flag, true); // say we've failed so that's the result if we throw
    SincFluxControl const & ctrl = static_cast<SincFluxControl const &>(getControl());

    afw::geom::ellipses::Axes const axes(ctrl.radius2, ctrl.radius2*(1.0 - ctrl.ellipticity), ctrl.angle);
    std::pair<double, double> fluxes =
        photometry::calculateSincApertureFlux(exposure.getMaskedImage(),
                                              afw::geom::ellipses::Ellipse(axes, center),
                                              ctrl.radius1/ctrl.radius2);
    double flux = fluxes.first;
    double fluxErr = fluxes.second;
    source.set(getKeys().meas, flux);
    source.set(getKeys().err, fluxErr);
    source.set(getKeys().flag, false);
}

#define INSTANTIATE(T) \
    template lsst::afw::image::Image<T>::Ptr detail::calcImageRealSpace<T>(double const, double const, \
                                                                           double const); \
    template lsst::afw::image::Image<T>::Ptr detail::calcImageKSpaceReal<T>(double const, double const); \
    template lsst::afw::image::Image<T>::Ptr detail::calcImageKSpaceCplx<T>(double const, double const, \
                                                                            double const, double const); \
    template std::pair<double, double> \
    photometry::calculateSincApertureFlux<lsst::afw::image::MaskedImage<T> >( \
        lsst::afw::image::MaskedImage<T> const&, lsst::afw::geom::ellipses::Ellipse const&, double const); \
    template class photometry::SincCoeffs<T>;

INSTANTIATE(float);
INSTANTIATE(double);
   
LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(SincFlux);
 
PTR(AlgorithmControl) SincFluxControl::_clone() const {
    return boost::make_shared<SincFluxControl>(*this);
}

PTR(Algorithm) SincFluxControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const &
) const {
    return boost::make_shared<SincFlux>(*this, boost::ref(schema));
}


}}}
