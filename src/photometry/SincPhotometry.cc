// -*- LSST-C++ -*-
/**
 * @file SincPhotometry.cc
 *
 * @brief Measure aperture photometry in a more sophisticated way than NaivePhotometry
 * @author Steve Bickerton
 * @ingroup meas/algorithms
 *
 */
#include <numeric>
#include <cmath>
#include <functional>
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image.h"
#include "lsst/meas/algorithms/Measure.h"
#include "lsst/afw/math/Integrate.h"

#include "SincPhotometry.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace detection = lsst::afw::detection;
namespace image = lsst::afw::image;
namespace math = lsst::afw::math;

namespace lsst { namespace meas { namespace algorithms {

/**
 * @brief the (unique) instance of measureSincPhotometry
 */
template<typename ImageT> measureSincPhotometry<ImageT>* measureSincPhotometry<ImageT>::_instance = 0;


namespace {

    // sinc function
    template<typename T>
    T sinc(T const x) {
        return (x != 0.0) ? (std::sin(x) / x) : 1.0;
    }

    // define a circular aperture function object g_i, cos-tapered?
    // (perhaps this can be in SincPhotometry.cc)
    template<typename CoordT>            
    class CircularAperture {
    public:
        
        CircularAperture(CoordT const radius, CoordT const taperwidth):
            _radius(radius), _taperwidth(taperwidth)
            {
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
    class SincAperture : public std::binary_function<IntegrandT,IntegrandT,IntegrandT> {
    public:
        
        SincAperture(CircularAperture<IntegrandT> &ap,
                     double const xcen,   // aperture center x
                     double const ycen,   // aperture center y
                     int const ix,        // sinc center x
                     int const iy         // sinc center y
                    )
            : _ap(ap),
              _xcen(xcen), _ycen(ycen)
            {
            
            _pixel_center_correction = 0.0;
            _ix = ix + _pixel_center_correction;
            _iy = iy + _pixel_center_correction;
            _xtaper = 10.0;
            _ytaper = 10.0;
        }

        IntegrandT operator() (IntegrandT const x, IntegrandT const y) const {
            double const FTconvention = 1.0*M_PI;
            //double const x = r * std::cos(t);
            //double const y = r * std::sin(t);
            double const xx = FTconvention*(x - _ix);
            double const yy = FTconvention*(y - _iy);
            double const fx = 0.5*(1.0 + std::cos(xx/_xtaper)) * sinc<double>(xx);
            double const fy = 0.5*(1.0 + std::cos(yy/_ytaper)) * sinc<double>(yy);
            return (1.0 + _ap(x - _xcen, y - _ycen)*fx*fy);
            //return sinc(FTconvention*(x - _ix)) * sinc(FTconvention*(_y - _iy));
        }
        
    private: 
        CircularAperture<IntegrandT> &_ap;
        double _xcen, _ycen;
        double _ix, _iy;
        double _xtaper, _ytaper;
        double _pixel_center_correction;
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

}
            
    template<typename PixelT>
    typename image::Image<PixelT>::Ptr getCoeffImage(
                                                     double const xcen0,
                                                     double const ycen0,
                                                     double const radius
                                                    ) {
        // @todo this should be in a .paf file with radius
        double const taperwidth = 2.0;
        double const buffer_width = 10.0;

        PixelT initweight = 0.0;
        double const xdwidth = 2.0*(radius + taperwidth + buffer_width);
        double const ydwidth = 2.0*(radius + taperwidth + buffer_width);
        int const xwidth = static_cast<int>(xdwidth);
        int const ywidth = static_cast<int>(ydwidth);
        double ip;
        double const xcen = static_cast<double>(xwidth/2) + std::modf(xcen0, &ip);
        double const ycen = static_cast<double>(ywidth/2) + std::modf(ycen0, &ip);
        
        // create an image to hold the coefficient image
        typename image::Image<PixelT>::Ptr cimage =
            typename image::Image<PixelT>::Ptr(new image::Image<PixelT>(xwidth, ywidth, initweight));

        // create the aperture function object
        CircularAperture<double> ap(radius, taperwidth);

        // ################################################################################
        // integrate over the aperture
        PixelT normalization_sum = 0.0;
        //double const epsilon = 0.01;
        double const limit = radius + taperwidth;
        double const x1 = xcen - limit;
        double const x2 = xcen + limit;
        double const y1 = ycen - limit;
        double const y2 = ycen + limit;
        for (int i_y = 0; i_y != cimage->getHeight(); ++i_y) {
            int i_x = 0;
            for (typename image::Image<PixelT>::x_iterator ptr = cimage->row_begin(i_y),
                     end = cimage->row_end(i_y);
                 ptr != end; ++ptr, ++i_x) {
                SincAperture<double> sincAp(ap, xcen, ycen, i_x, i_y);
                PixelT integral = math::integrate2d(sincAp, x1, x2, y1, y2, 1.0e-8);
                
                // we actually integrated 1+function and now must subtract the excess volume
                double const dx = i_x - xcen;
                double const dy = i_y - ycen;
                if ( std::sqrt(dx*dx + dy*dy) > xwidth/2) {
                    *ptr = 0.0;
                } else {
                    *ptr = integral - (x2 - x1)*(y2 - y1);
                    normalization_sum += integral;
                }
            }
        }

        // normalize
        PixelT const normalization_factor = 1.0; ///normalization_sum; //M_PI*radius*radius/normalization_sum;
        for (int i_y = 0; i_y != cimage->getHeight(); ++i_y) {
            int i_x = 0;
            typename image::Image<PixelT>::x_iterator end = cimage->row_end(i_y);
            for (typename image::Image<PixelT>::x_iterator ptr = cimage->row_begin(i_y);
                 ptr != end; ++ptr) {
                *ptr *= normalization_factor;
                ++i_x;
            }
        }
        //cimage->writeFits("cimage.fits");
        return cimage;
    }
    

namespace {
template <typename MaskedImageT, typename WeightImageT>
class FootprintWeightFlux : public detection::FootprintFunctor<MaskedImageT> {
public:
    FootprintWeightFlux(MaskedImageT const& mimage, ///< The image the source lives in
                        typename WeightImageT::Ptr wimage    ///< The weight image
                       ) : detection::FootprintFunctor<MaskedImageT>(mimage),
                           _wimage(wimage),
                           _sum(0), _x0(0), _y0(0)
        {}
    
    /// @brief Reset everything for a new Footprint
    void reset(detection::Footprint const& foot) {
        _sum = 0.0;

        lsst::afw::image::BBox const& bbox(foot.getBBox());
        _x0 = bbox.getX0();
        _y0 = bbox.getY0();

        if (bbox.getDimensions() != _wimage->getDimensions()) {
            throw LSST_EXCEPT(lsst::pex::exceptions::LengthErrorException,
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

}            

            
/**
 * @brief Given an image and a pixel position, return a Photometry 
 */
template<typename MaskedImageT>
Photometry measureSincPhotometry<MaskedImageT>::doApply(MaskedImageT const& mimage, ///< The MaskedImage
                                                        double xcen,          ///< object's column position
                                                        double ycen,          ///< object's row position
                                                        PSF const *psf,       ///< mimage's PSF
                                                        double background     ///< mimage's background level
                                                       ) const {

    typedef typename MaskedImageT::Image::Pixel Pixel;
    typedef typename image::Image<Pixel>::Image Image;
    typedef typename image::Image<Pixel>::Ptr   ImagePtr;
    
    Photometry photometry;              // The photometry to return
    
    int const ixcen = image::positionToIndex(xcen);
    int const iycen = image::positionToIndex(ycen);

    image::BBox imageBBox(image::PointI(mimage.getX0(), mimage.getY0()),
                          mimage.getWidth(), mimage.getHeight()); // BBox for data image

    static double last_radius = this->_radius;
    
    // Aperture photometry
    {
        // make the coeff image
        // compute c_i as double integral over aperture def g_i(), and sinc()
        static typename MaskedImageT::ImagePtr cimage0 = getCoeffImage<Pixel>(0, 0, this->_radius);

        if ( last_radius != this->_radius ) {
            cimage0 = getCoeffImage<Pixel>(0, 0, this->_radius);
            last_radius = this->_radius;
        }
        
        // shift it by the appropriate fractional pixel
        double dummy;
        double const dxpix = std::modf(xcen, &dummy);
        double const dypix = std::modf(ycen, &dummy);
        ImagePtr cimage_tmp = math::offsetImage(*cimage0, dxpix - 0.5, dypix - 0.5);
        
        int const border = 5;
        typename image::BBox coeffBBox(image::PointI(cimage_tmp->getX0() + border, cimage_tmp->getY0() +
                                                     border), cimage_tmp->getWidth() - 2*border,
                                       cimage_tmp->getHeight() - 2*border);
        
        ImagePtr cimage = ImagePtr(new Image(*cimage_tmp, coeffBBox, true));
        
        cimage->writeFits("cimage.fits");
        cimage0->writeFits("cimage0.fits");
        
        // pass the image and cimage into the wfluxFunctor to do the sum
        FootprintWeightFlux<MaskedImageT,typename MaskedImageT::Image> wfluxFunctor(mimage, cimage);
        detection::Footprint foot(image::BBox(image::PointI(cimage->getX0(), cimage->getY0()),
                                              cimage->getWidth(), cimage->getHeight()), imageBBox);
        foot.shift(ixcen - cimage->getWidth()/2, iycen - cimage->getHeight()/2);
        wfluxFunctor.apply(foot);

        // if the cimage is correctly made, the sum should be 1.0 ... keep for debugging and then remove this
        getSum2<Pixel> csum;
        csum = std::accumulate(cimage->begin(true), cimage->end(true), csum);

        photometry.setApFlux( wfluxFunctor.getSum() );
    }

    // Weighted aperture photometry, using a PSF weight --- i.e. a PSF flux
    {


        PSF::Image::Ptr wimage = psf->getImage(xcen, ycen);
        
        FootprintWeightFlux<MaskedImageT, PSF::Image> wfluxFunctor(mimage, wimage);
        // Build a rectangular Footprint corresponding to wimage
        detection::Footprint foot(image::BBox(image::PointI(0, 0), psf->getWidth(),
                                  psf->getHeight()), imageBBox);
        foot.shift(ixcen - psf->getWidth()/2, iycen - psf->getHeight()/2);

        wfluxFunctor.apply(foot);

        getSum2<PSF::Pixel> sum;
        sum = std::accumulate(wimage->begin(true), wimage->end(true), sum);

        photometry.setPsfFlux( wfluxFunctor.getSum()/sum.sum2 );

        //photometry.setPsfFlux(1.0);
    }
    
    return photometry;
    
}

// Explicit instantiations
//
// We need to make an instance here so as to register it with measurePhotometry
//
// \cond
#define MAKE_PHOTOMETRYFINDERS(IMAGE_T)                                 \
            namespace {                                                 \
                measurePhotometry<lsst::afw::image::MaskedImage<IMAGE_T> >* foo = \
                    measureSincPhotometry<lsst::afw::image::MaskedImage<IMAGE_T> >::getInstance(0); \
            }
            
//MAKE_PHOTOMETRYFINDERS(double)
MAKE_PHOTOMETRYFINDERS(float)

// \endcond

}}}
