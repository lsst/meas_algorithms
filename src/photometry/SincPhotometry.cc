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
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image.h"
#include "lsst/meas/algorithms/Measure.h"
#include "lsst/afw/math/Quadrature.h"

#include "SincPhotometry.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace detection = lsst::afw::detection;
namespace image = lsst::afw::image;
namespace math = lsst::afw::math;

namespace lsst { namespace meas { namespace algorithms {

/**
 * @brief the (unique) instance of measureNaivePhotometry
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
    template<typename T>            
    class CircularAperture {
    public:
        
        CircularAperture(T const radius, T const taperwidth):
            _radius(radius), _taperwidth(taperwidth)
            {
                _period = 4.0 * taperwidth;
            }
        

        // replace the sinusoid taper with a band-limited
        T operator() (T const x, T const y) {
            T xyrad = std::sqrt(x*x + y*y);
            if ( xyrad <= _radius - _taperwidth ) {
                return 1.0;
            } else if (xyrad > (_radius - _taperwidth) && xyrad <= (_radius + _taperwidth) ) {
                return 0.5*(1.0 - std::sin(  (2.0*M_PI)/_period * (xyrad - _radius)));
            } else {
                return 0.0;
            }
        }
        
    private:
        T _radius;
        T _taperwidth;
        T _period;
    };

    template<typename IntegrandT>
    class SincAperture : public math::IntegrandBase {
    public:
        
        SincAperture(CircularAperture<IntegrandT> &ap,
                     double const xcen, double const ycen,
                     int const ix, int const iy)
            : _ap(ap),
              _xcen(xcen), _ycen(ycen),
              _ix(ix), _iy(iy)
            {}

        IntegrandT getY() { return math::IntegrandBase::_y; }
        
        IntegrandT operator() (IntegrandT const x) {
            double const FTconvention = 1.0*M_PI;
            return (1.0 + _ap(x - _xcen, _y - _ycen) * sinc(FTconvention*(x - _ix)) * sinc(FTconvention*(_y - _iy)));
            //return sinc(FTconvention*(x - _ix)) * sinc(FTconvention*(_y - _iy));
        }
        
    private:
        CircularAperture<IntegrandT> &_ap;
        double _xcen, _ycen;
        double _ix, _iy;
        using math::IntegrandBase::_y;
    };
    

/************************************************************************************************************/
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
        double const taperwidth = 0.0;

        PixelT initweight = 0.0;
        int const xwidth = 2*static_cast<int>(radius + 2*taperwidth + 1);
        int const ywidth = 2*static_cast<int>(radius + 2*taperwidth + 1);
        double ip;
        double const xcen = static_cast<double>(xwidth/2) + std::modf(xcen0, &ip);
        double const ycen = static_cast<double>(ywidth/2) + std::modf(ycen0, &ip);
        
        // make the image
        typename image::Image<PixelT>::Ptr cimage =
            typename image::Image<PixelT>::Ptr(new image::Image<PixelT>(xwidth, ywidth, initweight));

        // create the aperture function object
        CircularAperture<double> ap(radius, taperwidth);

        // ################################################################################
        // integrate over the aperture
        for (int i_y = 0; i_y != cimage->getHeight(); ++i_y) {
            int i_x = 0;
            for (typename image::Image<PixelT>::x_iterator ptr = cimage->row_begin(i_y), end = cimage->row_end(i_y);
                 ptr != end; ++ptr, ++i_x) {
                SincAperture<double> sincAp(ap, xcen, ycen, i_x, i_y);
                PixelT integral = math::romb2D(sincAp, 0, xwidth - 1, 0, ywidth - 1);

                // we actually integrated 1+function and now must subtract the excess volume
                *ptr = integral - (xwidth-1.0)*(ywidth-1.0);
                //*ptr = ap(i_x - xcen,i_y-ycen);
                //SincAperture<double> sincAp(ap, xcen, ycen, xcen, xcen);
                //sincAp.setY(i_y);
                //*ptr = sincAp(i_x);
            }
        }
        
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
                              (boost::format("Footprint at %d,%d -- %d,%d is wrong size for %d x %d weight image") %
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
Photometry measureSincPhotometry<MaskedImageT>::doApply(MaskedImageT const& mimage, ///< The MaskedImage wherein dwells the object
                                                        double xcen,          ///< object's column position
                                                        double ycen,          ///< object's row position
                                                        PSF const* psf,       ///< mimage's PSF
                                                        double background     ///< mimage's background level
                                                       ) const {
    Photometry photometry;              // The photometry to return
    
    int const ixcen = image::positionToIndex(xcen);
    int const iycen = image::positionToIndex(ycen);

    image::BBox imageBBox(image::PointI(mimage.getX0(), mimage.getY0()),
                          mimage.getWidth(), mimage.getHeight()); // BBox for data image

    // Aperture photometry
    {
        // make the coeff image
        // compute c_i as double integral over aperture def g_i(), and sinc()
        typename MaskedImageT::ImagePtr cimage =
            getCoeffImage<typename MaskedImageT::Image::Pixel>(xcen, ycen, this->_radius);

        cimage->writeFits("cimage.fits");
        
        // pass the image and cimage into the wfluxFunctor to do the sum
        FootprintWeightFlux<MaskedImageT,typename MaskedImageT::Image> wfluxFunctor(mimage, cimage);
        detection::Footprint foot(image::BBox(image::PointI(0, 0), cimage->getWidth(), cimage->getHeight()), imageBBox);
        foot.shift(ixcen - cimage->getWidth()/2, iycen - cimage->getHeight()/2);
        wfluxFunctor.apply(foot);
        
        // if the cimage is correctly made, the sum should be 1.0 ... keep for debugging and then remove this
        getSum2<typename MaskedImageT::Image::Pixel> sum;
        sum = std::accumulate(cimage->begin(true), cimage->end(true), sum);
        
        photometry.setApFlux( wfluxFunctor.getSum()/sum.sum2 );
    }

    // Weighted aperture photometry, using a PSF weight --- i.e. a PSF flux
    {
        PSF::ImageT::Ptr wimage = psf->getImage(xcen, ycen);
        
        FootprintWeightFlux<MaskedImageT, PSF::ImageT> wfluxFunctor(mimage, wimage);
        // Build a rectangular Footprint corresponding to wimage
        detection::Footprint foot(image::BBox(image::PointI(0, 0), psf->getWidth(), psf->getHeight()), imageBBox);
        foot.shift(ixcen - psf->getWidth()/2, iycen - psf->getHeight()/2);

        wfluxFunctor.apply(foot);

        getSum2<PSF::PixelT> sum;
        sum = std::accumulate(wimage->begin(true), wimage->end(true), sum);

        photometry.setPsfFlux( wfluxFunctor.getSum()/sum.sum2 );
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
