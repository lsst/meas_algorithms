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
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/geom/Extent.h"
#include "lsst/afw/image.h"
#include "lsst/afw/math/Integrate.h"
#include "lsst/meas/algorithms/Measure.h"

#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/detection/Photometry.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace afwDetection = lsst::afw::detection;
namespace afwImage = lsst::afw::image;
namespace afwMath = lsst::afw::math;

namespace lsst {
namespace meas {
namespace algorithms {

/// primarily for debug
template<typename PixelT>
typename lsst::afw::image::Image<PixelT>::Ptr getCoeffImage(double const xcen0,
                                                            double const ycen0,
                                                            double const radius);

/**
 * @brief A class that knows how to calculate fluxes using the SINC photometry algorithm
 * @ingroup meas/algorithms
 */
class SincPhotometry : public afwDetection::Photometry
{
public:
    typedef boost::shared_ptr<SincPhotometry> Ptr;
    typedef boost::shared_ptr<SincPhotometry const> ConstPtr;

    /// Ctor
    SincPhotometry(double flux, double fluxErr=std::numeric_limits<double>::quiet_NaN()) :
        afwDetection::Photometry(flux, fluxErr) {}

    /// Add desired fields to the schema
    virtual void defineSchema(afwDetection::Schema::Ptr schema ///< our schema; == _mySchema
                     ) {
        Photometry::defineSchema(schema);
    }

    static bool doConfigure(lsst::pex::policy::Policy const& policy);

    template<typename ImageT>
    static Photometry::Ptr doMeasure(typename ImageT::ConstPtr im, afwDetection::Peak const*);

    /// Set the aperture radius to use
    static void setRadius(double radius) { _radius = radius; }

    /// Return the aperture radius to use
    static double getRadius() { return _radius; }

private:
    static double _radius;
};

double SincPhotometry::_radius = 0;      // radius to use for sinc photometry


    
/************************************************************************************************************/
namespace {

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
                     CoordT const xcen,      ///< xcenter coord
                     CoordT const ycen,      ///< ycenter coord
                     CoordT const radius,    ///< radius of the aperture
                     CoordT const taperwidth ///< width to cosine taper from 1.0 to 0.0 (ie. 0.5*cosine period)
                    ):
        _xcen(xcen), _ycen(ycen),
        _radius(radius),
        _taperwidth(taperwidth),
        _k(1.0/(2.0*taperwidth)),
        _taperLo(radius-0.5*taperwidth),
        _taperHi(radius+0.5*taperwidth)
        {

        // if we're asked for a radius smaller than our taperwidth,
        // adjust the taper width smaller so it fits exactly
        // with smooth derivative=0 at r=0
        if (radius < 0.5*_taperwidth) {
            _taperwidth = 2.0*radius;
            _k = 1.0/(2.0*_taperwidth); 
            _taperLo = _radius - 0.5*_taperwidth; 
            _taperHi = _radius + 0.5*_taperwidth;
        }
        
    }
    

    // When called, return the throughput at the requested x,y
    // todo: replace the sinusoid taper with a band-limited
    CoordT operator() (CoordT const x, CoordT const y) const {
        CoordT xx = (x - _xcen)*(x - _xcen);
        CoordT yy = (y - _ycen)*(y - _ycen);
        CoordT const xyrad = std::sqrt(xx + yy);
        if ( xyrad <= _taperLo ) {
            return 1.0;
        } else if (xyrad > _taperLo && xyrad <= _taperHi ) {
            return 0.5*(1.0 + std::cos(  (2.0*M_PI*_k)*(xyrad - _taperLo)) );
        } else {
            return 0.0;
        }
    }

    CoordT getRadius() { return _radius; }
    CoordT getXcen() { return _xcen; }
    CoordT getYcen() { return _ycen; }
    
private:
    CoordT _xcen, _ycen;
    CoordT _radius;
    CoordT _taperwidth;
    CoordT _k;            // the angular wavenumber corresponding to a cosine with wavelength 2*taperwidth
    CoordT _taperLo;
    CoordT _taperHi;
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
        : _ap(ap), _ix(ix), _iy(iy) {
        _xtaper = 10.0;
        _ytaper = 10.0;
    }
    
    IntegrandT operator() (IntegrandT const x, IntegrandT const y) const {
        double const fourierConvention = 1.0*M_PI;
        double const dx = fourierConvention*(x - _ix);
        double const dy = fourierConvention*(y - _iy);
        double const fx = 0.5*(1.0 + std::cos(dx/_xtaper)) * sinc<double>(dx);
        double const fy = 0.5*(1.0 + std::cos(dy/_ytaper)) * sinc<double>(dy);
        return (1.0 + _ap(x, y)*fx*fy);
    }
    
private: 
    CircularAperture<IntegrandT> const &_ap;
    double _ix, _iy;
    double _xtaper, _ytaper; // x,y distances over which to cos-taper the sinc to zero
};
    


/******************************************************************************/
    
template <typename MaskedImageT, typename WeightImageT>
class FootprintWeightFlux : public afwDetection::FootprintFunctor<MaskedImageT> {
public:
    FootprintWeightFlux(
                        MaskedImageT const& mimage, ///< The image the source lives in
                        typename WeightImageT::Ptr wimage    ///< The weight image
                       ) :
        afwDetection::FootprintFunctor<MaskedImageT>(mimage),
        _wimage(wimage),
        _sum(0.0), _sumVar(0.0),
        _x0(wimage->getX0()), _y0(wimage->getY0()) {}
    
    /// @brief Reset everything for a new Footprint
    void reset() {}        
    void reset(afwDetection::Footprint const& foot) {
        _sum = 0.0;
        _sumVar = 0.0;

        afwImage::BBox const& bbox(foot.getBBox());
        _x0 = bbox.getX0();
        _y0 = bbox.getY0();

        if (bbox.getDimensions() != _wimage->getDimensions()) {
            throw LSST_EXCEPT(pexExceptions::LengthErrorException,
                              (boost::format("Footprint at %d,%d -- %d,%d is wrong size "
                                             "for %d x %d weight image") %
                               bbox.getX0() % bbox.getY0() % bbox.getX1() % bbox.getY1() %
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
    typename WeightImageT::Ptr const& _wimage;        // The weight image
    double _sum;                                      // our desired sum
    double _sumVar;                                   // sum of the variance
    int _x0, _y0;                                     // the origin of the current Footprint
};

    
} // end of anonymous namespace



/************************************************************************************************************/
    
/**
 * @brief get the coefficient image for sinc aperture photometry
 */
template<typename PixelT>
typename afwImage::Image<PixelT>::Ptr getCoeffImage(
                                                 double const xcen,
                                                 double const ycen,
                                                 double const radius
                                                ) {
    // @todo this should be in a .paf file with radius
    double const taperwidth = 2.0;    // for the cosine taper
    double const bufferWidth = 10.0;  // extra pixels to provide a border
    
    PixelT initweight = 0.0;          // initialize the coeff values to this

    // widths as doubles
    double const xdwidth = 2.0*(radius + taperwidth + bufferWidth);
    double const ydwidth = 2.0*(radius + taperwidth + bufferWidth);

    // widths as ints
    int const xwidth = static_cast<int>(xdwidth) + 1;
    int const ywidth = static_cast<int>(ydwidth) + 1;

    // the origin -
    // xwidth is guaranteed to be odd, so this always places x0,y0 correctly
    //   to ensure the center is actually at xcen, ycen
    int const x0 = static_cast<int>(xcen) - xwidth/2;
    int const y0 = static_cast<int>(ycen) - ywidth/2;
    
    // create an image to hold the coefficient image
    typename afwImage::Image<PixelT>::Ptr cimage =
        typename afwImage::Image<PixelT>::Ptr(new afwImage::Image<PixelT>(xwidth, ywidth, initweight));
    cimage->setXY0(x0, y0);

    
    // create the aperture function object
    CircularAperture<double> ap(xcen, ycen, radius, taperwidth);
    

    /* ******************************* */
    // integrate over the aperture
    
    // the limits of the integration over the sinc aperture
    double const limit = radius + taperwidth;
    double const x1 = xcen - limit;
    double const x2 = xcen + limit;
    double const y1 = ycen - limit;
    double const y2 = ycen + limit;
    
    for (int iY = y0; iY != y0 + cimage->getHeight(); ++iY) {
        int iX = x0;
        typename afwImage::Image<PixelT>::x_iterator end = cimage->row_end(iY-y0);
        for (typename afwImage::Image<PixelT>::x_iterator ptr = cimage->row_begin(iY-y0); ptr != end; ++ptr) {

            // create a sinc function in the CircularAperture at our location
            SincAperture<double> sincAp(ap, iX, iY);

            // integrate the sinc
            PixelT integral = afwMath::integrate2d(sincAp, x1, x2, y1, y2, 1.0e-8);
            
            // we actually integrated function+1.0 and now must subtract the excess volume
            // - just force it to zero in the corners
            double const dx = iX - xcen;
            double const dy = iY - ycen;
            *ptr = (std::sqrt(dx*dx + dy*dy) < xwidth/2) ?
                integral - (x2 - x1)*(y2 - y1) : 0.0;
            
            ++iX;
        }
    }
    
    return cimage;
}

    
/************************************************************************************************************/
/**
 * Set parameters controlling how we do measurements
 */
bool SincPhotometry::doConfigure(lsst::pex::policy::Policy const& policy)
{
    if (policy.isDouble("radius")) {
        setRadius(policy.getDouble("radius"));
    } 

    return true;
}
    
/************************************************************************************************************/
/**
 * Calculate the desired aperture flux using the sinc algorithm
 */
template<typename ExposureT>
afwDetection::Photometry::Ptr SincPhotometry::doMeasure(typename ExposureT::ConstPtr exposure,
                                                        afwDetection::Peak const* peak
                                                       ) {
    
    double flux = std::numeric_limits<double>::quiet_NaN();
    double fluxErr = std::numeric_limits<double>::quiet_NaN();
    if (!peak) {
        return boost::make_shared<SincPhotometry>(flux, fluxErr);
    }
    
    typedef typename ExposureT::MaskedImageT MaskedImageT;
    typedef typename MaskedImageT::Image Image;
    typedef typename Image::Pixel Pixel;
    typedef typename Image::Ptr ImagePtr;

    MaskedImageT const& mimage = exposure->getMaskedImage();
    
    double const xcen = peak->getFx();   ///< object's column position
    double const ycen = peak->getFy();   ///< object's row position
    
    afwImage::BBox imageBBox(afwImage::PointI(mimage.getX0(), mimage.getY0()),
                             mimage.getWidth(), mimage.getHeight()); // BBox for data image
    
    static double last_radius = getRadius();

    /* ********************************************************** */
    // Aperture photometry
    {
        // make the coeff image
        // compute c_i as double integral over aperture def g_i(), and sinc()
        // make static so we can reuse it ... it's a bit costly.
        static ImagePtr cimage0 = getCoeffImage<Pixel>(0, 0, getRadius());

        // as long as we're asked for the same radius, we don't have to recompute cimage0
        if (::fabs(last_radius - getRadius()) > std::numeric_limits<double>::epsilon()) {
            cimage0 = getCoeffImage<Pixel>(0, 0, getRadius());
            last_radius = getRadius();
        }
        cimage0->markPersistent();
        
        // shift to center the aperture on the object being measured
        ImagePtr cimage = afwMath::offsetImage(*cimage0, xcen, ycen);
        
        // pass the image and cimage into the wfluxFunctor to do the sum
        FootprintWeightFlux<MaskedImageT, Image> wfluxFunctor(mimage, cimage);
        afwDetection::Footprint foot(afwImage::BBox(afwImage::PointI(cimage->getX0(), cimage->getY0()),
                                                    cimage->getWidth(), cimage->getHeight()), imageBBox);
        wfluxFunctor.apply(foot);
        flux = wfluxFunctor.getSum();
        fluxErr = wfluxFunctor.getSumVar();
    }
    return boost::make_shared<SincPhotometry>(flux, fluxErr);
}

//
// Explicit instantiations
//
// \cond
#define INSTANTIATE(T) \
    template lsst::afw::image::Image<T>::Ptr getCoeffImage<T>(double const, double const, double const)
    
/*
 * Declare the existence of a "SINC" algorithm to MeasurePhotometry
 */
#define MAKE_PHOTOMETRYS(TYPE)                                          \
    MeasurePhotometry<afwImage::Exposure<TYPE> >::declare("SINC", \
        &SincPhotometry::doMeasure<afwImage::Exposure<TYPE> >, \
        &SincPhotometry::doConfigure \
    )

namespace {
    volatile bool isInstance[] = {
        MAKE_PHOTOMETRYS(float)
#if 0
        ,MAKE_PHOTOMETRYS(double)
#endif
    };
}

INSTANTIATE(float);
#if 0
INSTANTIATE(double);
#endif
    
// \endcond

}}}
