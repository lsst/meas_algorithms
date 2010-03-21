// -*- LSST-C++ -*-
#include <limits>
#include <numeric>
#include "Eigen/LU"
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image.h"
#include "lsst/meas/algorithms/Measure.h"

#include "lsst/meas/algorithms/Photometry.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace detection = lsst::afw::detection;
namespace afwImage = lsst::afw::image;

namespace lsst {
namespace meas {
namespace algorithms {

namespace {
/**
 * @brief A class that knows how to calculate photometrys as a simple unweighted first moment
 * of the 3x3 region around a pixel
 */
template<typename ImageT>
class NaiveMeasurePhotometry : public MeasurePhotometry<ImageT> {
public:
    typedef MeasurePhotometry<ImageT> MeasurePropertyBase;

    using MeasurePhotometry<ImageT>::getRadius;

    NaiveMeasurePhotometry(typename ImageT::ConstPtr image) : MeasurePhotometry<ImageT>(image) {}
private:
    Photometry doApply(ImageT const& image, double xcen, double ycen,
                       PSF const* psf, double background) const;
};

namespace {
            
template <typename MaskedImageT>
class FootprintFlux : public detection::FootprintFunctor<MaskedImageT> {
public:
    explicit FootprintFlux(MaskedImageT const& mimage ///< The image the source lives in
                 ) : detection::FootprintFunctor<MaskedImageT>(mimage),
                     _sum(0) {}

    /// @brief Reset everything for a new Footprint
    void reset() {
        _sum = 0.0;
    }
    void reset(detection::Footprint const&) {}        

    /// @brief method called for each pixel by apply()
    void operator()(typename MaskedImageT::xy_locator loc, ///< locator pointing at the pixel
                    int,                                   ///< column-position of pixel
                    int                                    ///< row-position of pixel
                   ) {
        typename MaskedImageT::Image::Pixel val = loc.image(0, 0);
        _sum += val;
    }

    /// Return the Footprint's flux
    double getSum() const { return _sum; }

private:
    double _sum;
};


template <typename MaskedImageT, typename WeightImageT>
class FootprintWeightFlux : public detection::FootprintFunctor<MaskedImageT> {
public:
    FootprintWeightFlux(MaskedImageT const& mimage,          ///< The image the source lives in
                        typename WeightImageT::Ptr wimage    ///< The weight image
                       ) : detection::FootprintFunctor<MaskedImageT>(mimage),
                           _wimage(wimage),
                           _sum(0), _x0(0), _y0(0) {}
    
    /// @brief Reset everything for a new Footprint
    void reset(detection::Footprint const& foot) {
        _sum = 0.0;
        
        afwImage::BBox const& bbox(foot.getBBox());
        _x0 = bbox.getX0();
        _y0 = bbox.getY0();

        if (bbox.getDimensions() != _wimage->getDimensions()) {
            throw LSST_EXCEPT(pexExceptions::LengthErrorException,
                              (boost::format("Footprint at %d,%d -- %d,%d is wrong size for "
                                             "%d x %d weight image") %
                               bbox.getX0() % bbox.getY0() % bbox.getX1() % bbox.getY1() %
                               _wimage->getWidth() % _wimage->getHeight()).str());
        }
    }
    void reset() {}
    
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

            
/*****************************************************************************************************/
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
    
} // end of anonymous namespace

    

    
/**
 * @brief Given an image and a pixel position, return a Photometry using a naive 3x3 weighted moment
 */
template<typename MaskedImageT>
Photometry NaiveMeasurePhotometry<MaskedImageT>::doApply(MaskedImageT const& img,   ///< The Image 
                                                   double xcen,    ///< object's column position
                                                   double ycen,    ///< object's row position
                                                   PSF const *psf, ///< image's PSF
                                                   double          ///< image's background level
                                                  ) const {

    Photometry photometry;              // The photometry to return
    
    int const ixcen = afwImage::positionToIndex(xcen);
    int const iycen = afwImage::positionToIndex(ycen);
    
    afwImage::BBox imageBBox(afwImage::PointI(img.getX0(), img.getY0()),
                          img.getWidth(), img.getHeight()); // BBox for data image

    /* ******************************************************* */
    // Aperture photometry
    {
        FootprintFlux<MaskedImageT> fluxFunctor(img);
        
        detection::Footprint const foot(afwImage::BCircle(afwImage::PointI(ixcen, iycen), getRadius()),
                                        imageBBox);
        fluxFunctor.apply(foot);
        photometry.setApFlux( fluxFunctor.getSum() );
    }


    /* ******************************************************** */
    // Weighted aperture photometry, using a PSF weight --- i.e. a PSF flux
    if (psf) {
        
        PSF::Image::Ptr wimage = psf->getImage(xcen, ycen);
        
        FootprintWeightFlux<MaskedImageT, PSF::Image> wfluxFunctor(img, wimage);
        
        // Build a rectangular Footprint corresponding to wimage
        detection::Footprint foot(afwImage::BBox(afwImage::PointI(0, 0),
                                              psf->getWidth(), psf->getHeight()), imageBBox);
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
    registerMe<NaiveMeasurePhotometry, afwImage::MaskedImage<IMAGE_T> >("NAIVE")
    
volatile bool isInstance[] = {
    MAKE_PHOTOMETRYS(float)
#if 0
    ,MAKE_PHOTOMETRYS(double)
#endif
};

// \endcond

}}}}
