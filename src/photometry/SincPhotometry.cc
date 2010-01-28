// -*- LSST-C++ -*-
#include <numeric>
#include <cmath>
#include <functional>
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image.h"
#include "lsst/meas/algorithms/Measure.h"
#include "lsst/afw/math/Integrate.h"

#include "lsst/meas/algorithms/Photometry.h"
#include "lsst/meas/algorithms/detail/SincPhotometry.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace detection = lsst::afw::detection;
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
template<typename ImageT>
class SincMeasurePhotometry : public MeasurePhotometry<ImageT> {
public:
    typedef MeasurePhotometry<ImageT> MeasurePropertyBase;

    using MeasurePhotometry<ImageT>::getRadius;

    SincMeasurePhotometry(typename ImageT::ConstPtr image) : MeasurePhotometry<ImageT>(image) {}
private:
    Photometry doApply(ImageT const& image, double xcen, double ycen,
                       PSF const* psf, double background) const;
};

/************************************************************************************************************/
namespace {

// sinc function
template<typename T>
inline T sinc(T const x) {
    return (x != 0.0) ? (std::sin(x) / x) : 1.0;
}

// define a circular aperture function object g_i, cos-tapered?
// (perhaps this can be in SincPhotometry.cc)
template<typename CoordT>            
class CircularAperture {
public:
    
    CircularAperture(CoordT const radius, CoordT const taperwidth):
        _radius(radius), _taperwidth(taperwidth) {
        _period = 2.0 * taperwidth;
    }
    
    
    // replace the sinusoid taper with a band-limited
    CoordT operator() (CoordT const x, CoordT const y) const {
        CoordT const xyrad = std::sqrt(x*x + y*y);
        if ( xyrad <= _radius ) {
            return 1.0;
        } else if (xyrad > (_radius) && xyrad < (_radius + _taperwidth) ) {
            return 0.5*(1.0 + std::cos(  (2.0*M_PI)/_period * (xyrad - _radius)) );
        } else {
            return 0.0;
        }
    }
    
private:
    CoordT _radius;
    CoordT _taperwidth;
    CoordT _period;
};
    
template<typename IntegrandT>
class SincAperture : public std::binary_function<IntegrandT, IntegrandT, IntegrandT> {
public:
    
    SincAperture(CircularAperture<IntegrandT> const &ap,
                 double const xcen,   // aperture center x
                 double const ycen,   // aperture center y
                 int const ix,        // sinc center x
                 int const iy         // sinc center y
                )
        : _ap(ap),
          _xcen(xcen), _ycen(ycen) {
        _pixelCenterCorrection = 0.0;
        _ix = ix + _pixelCenterCorrection;
        _iy = iy + _pixelCenterCorrection;
        _xtaper = 10.0;
        _ytaper = 10.0;
    }
    
    IntegrandT operator() (IntegrandT const x, IntegrandT const y) const {
        double const fTransConvention = 1.0*M_PI;
        //double const x = r * std::cos(t);
        //double const y = r * std::sin(t);
        double const xx = fTransConvention*(x - _ix);
        double const yy = fTransConvention*(y - _iy);
        double const fx = 0.5*(1.0 + std::cos(xx/_xtaper)) * sinc<double>(xx);
        double const fy = 0.5*(1.0 + std::cos(yy/_ytaper)) * sinc<double>(yy);
        return (1.0 + _ap(x - _xcen, y - _ycen)*fx*fy);
        //return sinc(fTransConvention*(x - _ix)) * sinc(fTransConvention*(_y - _iy));
    }
    
private: 
    CircularAperture<IntegrandT> const &_ap;
    double _xcen, _ycen;
    double _ix, _iy;
    double _xtaper, _ytaper;
    double _pixelCenterCorrection;
};
    
    
/***********************************************************************************************************/
/**
 * Accumulate sum(x) and sum(x**2)
 */
template<typename T>
struct getSum2 {
    getSum2() : sum(0.0), sum2(0.0) {}
    
    getSum2& operator+(T x) {
        sum += x;
        sum2 += x*x;
        
        return *this;
    }
    
    double sum;                         // \sum_i(x_i)
    double sum2;                        // \sum_i(x_i^2)
};


template <typename MaskedImageT, typename WeightImageT>
class FootprintWeightFlux : public detection::FootprintFunctor<MaskedImageT> {
public:
    FootprintWeightFlux(MaskedImageT const& mimage, ///< The image the source lives in
                        typename WeightImageT::Ptr wimage    ///< The weight image
                       ) : detection::FootprintFunctor<MaskedImageT>(mimage),
                           _wimage(wimage),
                           _sum(0), _x0(0), _y0(0) {}
    
    /// @brief Reset everything for a new Footprint
    void reset() {}        
    void reset(detection::Footprint const& foot) {
        _sum = 0.0;

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
        typename WeightImageT::Pixel wval = (*_wimage)(x - _x0, y - _y0);
        _sum += wval*ival;
    }

    /// Return the Footprint's flux
    double getSum() const { return _sum; }

private:
    typename WeightImageT::Ptr const& _wimage;        // The weight image
    double _sum;                                      // our desired sum
    int _x0, _y0;                                     // the origin of the current Footprint
};

    
} // end of anonymous namespace




template<typename PixelT>
typename afwImage::Image<PixelT>::Ptr getCoeffImage(
                                                 double const xcen0,
                                                 double const ycen0,
                                                 double const radius
                                                ) {
    // @todo this should be in a .paf file with radius
    double const taperwidth = 2.0;
    double const bufferWidth = 10.0;
    
    PixelT initweight = 0.0;
    double const xdwidth = 2.0*(radius + taperwidth + bufferWidth);
    double const ydwidth = 2.0*(radius + taperwidth + bufferWidth);
    int const xwidth = static_cast<int>(xdwidth);
    int const ywidth = static_cast<int>(ydwidth);
    double ip;
    double const xcen = static_cast<double>(xwidth/2) + std::modf(xcen0, &ip);
    double const ycen = static_cast<double>(ywidth/2) + std::modf(ycen0, &ip);
    
    // create an image to hold the coefficient image
    typename afwImage::Image<PixelT>::Ptr cimage =
        typename afwImage::Image<PixelT>::Ptr(new afwImage::Image<PixelT>(xwidth, ywidth, initweight));
    
    // create the aperture function object
    CircularAperture<double> ap(radius, taperwidth);
    
    // ################################################################################
    // integrate over the aperture
    PixelT normalizationSum = 0.0;
    //double const epsilon = 0.01;
    double const limit = radius + taperwidth;
    double const x1 = xcen - limit;
    double const x2 = xcen + limit;
    double const y1 = ycen - limit;
    double const y2 = ycen + limit;
    for (int iY = 0; iY != cimage->getHeight(); ++iY) {
        int iX = 0;
        typename afwImage::Image<PixelT>::x_iterator end = cimage->row_end(iY);
        for (typename afwImage::Image<PixelT>::x_iterator ptr = cimage->row_begin(iY); ptr != end; ++ptr) {
            SincAperture<double> sincAp(ap, xcen, ycen, iX, iY);
            PixelT integral = afwMath::integrate2d(sincAp, x1, x2, y1, y2, 1.0e-8);
            
            // we actually integrated 1+function and now must subtract the excess volume
            double const dx = iX - xcen;
            double const dy = iY - ycen;
            if ( std::sqrt(dx*dx + dy*dy) > xwidth/2) {
                *ptr = 0.0;
            } else {
                *ptr = integral - (x2 - x1)*(y2 - y1);
                normalizationSum += integral;
            }
            ++iX;
        }
    }
    
    // normalize
    PixelT const normalizationFactor = 1.0; ///normalizationSum; //M_PI*radius*radius/normalizationSum;
    for (int iY = 0; iY != cimage->getHeight(); ++iY) {
        int iX = 0;
        typename afwImage::Image<PixelT>::x_iterator end = cimage->row_end(iY);
        for (typename afwImage::Image<PixelT>::x_iterator ptr = cimage->row_begin(iY);
             ptr != end; ++ptr) {
            *ptr *= normalizationFactor;
            ++iX;
        }
    }
    //cimage->writeFits("cimage.fits");
    return cimage;
}



    
    
/**
 * @brief Given an image and a pixel position, return a Photometry using a naive 3x3 weighted moment
 */
template<typename MaskedImageT>
Photometry SincMeasurePhotometry<MaskedImageT>::doApply(MaskedImageT const& img,    ///< The Image 
                                                  double xcen,            ///< object's column position
                                                  double ycen,            ///< object's row position
                                                  PSF const *psf,         ///< image's PSF
                                                  double                  ///< image's background level
                                                 ) const {

    typedef typename MaskedImageT::Image::Pixel Pixel;
    typedef typename afwImage::Image<Pixel> Image;
    typedef typename afwImage::Image<Pixel>::Ptr ImagePtr;
    
    Photometry photometry;              // The photometry to return
    
    int const ixcen = afwImage::positionToIndex(xcen);
    int const iycen = afwImage::positionToIndex(ycen);

    afwImage::BBox imageBBox(afwImage::PointI(img.getX0(), img.getY0()),
                          img.getWidth(), img.getHeight()); // BBox for data image

    static double last_radius = getRadius();

    /* ********************************************************** */
    // Aperture photometry
    {
        // make the coeff image
        // compute c_i as double integral over aperture def g_i(), and sinc()
        static ImagePtr cimage0 = getCoeffImage<Pixel>(0, 0, getRadius());
        
        if (::fabs(last_radius - getRadius()) > std::numeric_limits<double>::epsilon()) {
            cimage0 = getCoeffImage<Pixel>(0, 0, getRadius());
            last_radius = getRadius();
        }
        
        // shift it by the appropriate fractional pixel
        double dummy;
        double const dxpix = std::modf(xcen, &dummy);
        double const dypix = std::modf(ycen, &dummy);
        ImagePtr cimage_tmp = afwMath::offsetImage(*cimage0, dxpix - 0.5, dypix - 0.5);
        
        int const border = 5;
        typename afwImage::BBox coeffBBox(afwImage::PointI(cimage_tmp->getX0() + border, cimage_tmp->getY0() +
                                                     border), cimage_tmp->getWidth() - 2*border,
                                       cimage_tmp->getHeight() - 2*border);
        
        ImagePtr cimage = ImagePtr(new Image(*cimage_tmp, coeffBBox, true));
        
        // pass the image and cimage into the wfluxFunctor to do the sum
        FootprintWeightFlux<MaskedImageT, typename MaskedImageT::Image> wfluxFunctor(img, cimage);
        detection::Footprint foot(afwImage::BBox(afwImage::PointI(cimage->getX0(), cimage->getY0()),
                                              cimage->getWidth(), cimage->getHeight()), imageBBox);
        foot.shift(ixcen - cimage->getWidth()/2, iycen - cimage->getHeight()/2);
        wfluxFunctor.apply(foot);
        
        // if the cimage is correctly made, the sum should be 1.0 ... keep for debugging and then remove this
        getSum2<Pixel> csum;
        csum = std::accumulate(cimage->begin(true), cimage->end(true), csum);
        
        photometry.setApFlux( wfluxFunctor.getSum() );
    }

    /* ***************************************************************** */
    // Weighted aperture photometry, using a PSF weight --- i.e. a PSF flux
    if (psf) {
        PSF::Image::Ptr wimage = psf->getImage(xcen, ycen);
        
        FootprintWeightFlux<MaskedImageT, PSF::Image> wfluxFunctor(img, wimage);
        // Build a rectangular Footprint corresponding to wimage
        detection::Footprint foot(afwImage::BBox(afwImage::PointI(0, 0), psf->getWidth(),
                                                 psf->getHeight()), imageBBox);
        foot.shift(ixcen - psf->getWidth()/2, iycen - psf->getHeight()/2);
        
        wfluxFunctor.apply(foot);
        
        getSum2<PSF::Pixel> sum;
        sum = std::accumulate(wimage->begin(true), wimage->end(true), sum);
        
        photometry.setPsfFlux( wfluxFunctor.getSum()/sum.sum2 );
    } else {
        photometry.setPsfFlux(std::numeric_limits<double>::quiet_NaN());
    }
    
    return photometry;
}

//
// Explicit instantiations
//
// We need to make an instance here so as to register it with MeasurePhotometry
//
// \cond
#define MAKE_PHOTOMETRYS(IMAGE_T)                                       \
    registerMe<SincMeasurePhotometry, afwImage::MaskedImage<IMAGE_T> >("SINC")

namespace {
    volatile bool isInstance[] = {
        MAKE_PHOTOMETRYS(float)
#if 0
        ,MAKE_PHOTOMETRYS(double)
#endif
    };
}

// \endcond

}}}
