// -*- LSST-C++ -*-
/**
 * @file NaivePhotometry.cc
 *
 * @brief Measure adaptive photometry in a simple way.
 * @author Steve Bickerton
 * @ingroup meas/algorithms
 *
 */
#include <limits>
#include "Eigen/LU"
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image.h"
#include "lsst/meas/algorithms/Measure.h"

#include "NaivePhotometry.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace detection = lsst::afw::detection;
namespace image = lsst::afw::image;

namespace lsst { namespace meas { namespace algorithms {

/**
 * @brief the (unique) instance of measureNaivePhotometry
 */
template<typename ImageT> measureNaivePhotometry<ImageT>* measureNaivePhotometry<ImageT>::_instance = 0;

namespace {
            
template <typename MaskedImageT>
class FootprintFlux : public detection::FootprintFunctor<MaskedImageT> {
public:
    FootprintFlux(MaskedImageT const& mimage                    ///< The image the source lives in
                 ) : detection::FootprintFunctor<MaskedImageT>(mimage),
                     _sum(0)
        {}

    /// @brief Reset everything for a new Footprint
    void reset() {
        _sum = 0.0;
    }

    /// @brief method called for each pixel by apply()
    void operator()(typename MaskedImageT::xy_locator loc, ///< locator pointing at the pixel
                    int x,                                 ///< column-position of pixel
                    int y                                  ///< row-position of pixel
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
        typename WeightImageT::Pixel wval = _wimage->at(x - _x0, y - _y0)[0];
        _sum += wval * ival;
    }

    /// Return the Footprint's flux
    double getSum() const { return _sum; }

private:
    typename WeightImageT::Ptr const& _wimage;        // The weight image
    double _sum;                                      // our desired sum
    int _x0, _y0;                                     // the origin of the current Footprint
};

}            

            
/************************************************************************************************************/
/**
 * @brief Given an image and a pixel position, return a Photometry 
 */
template<typename MaskedImageT>
Photometry measureNaivePhotometry<MaskedImageT>::doApply(MaskedImageT const& mimage, ///< The MaskedImage wherein dwells the object
                                              double xcen,          ///< object's column position
                                              double ycen,          ///< object's row position
                                              PSF const*,           ///< mimage's PSF
                                              double background     ///< mimage's background level
                                             ) const {
    xcen -= mimage.getX0();              // work in image Pixel coordinates
    ycen -= mimage.getY0();
    
    Photometry photometry;                         // The photometry to return


    // Aperture photometry
    FootprintFlux<MaskedImageT> fluxFunctor(mimage);
    image::PointI center(image::positionToIndex(xcen), image::positionToIndex(ycen));
    detection::Footprint const foot(image::BCircle(center, this->_radius));
    fluxFunctor.apply(foot);
    photometry.setApFlux( fluxFunctor.getSum() );

    // Weighted aperture photometry
    typedef typename MaskedImageT::Image WeightImageT;
    typename WeightImageT::Ptr wimage(new WeightImageT(foot.getBBox().getDimensions()));
    *wimage = 1;
    FootprintWeightFlux<MaskedImageT, WeightImageT> wfluxFunctor(mimage, wimage);
    wfluxFunctor.apply(foot);
    photometry.setPsfFlux( wfluxFunctor.getSum() );

    return photometry;
}

//
// Explicit instantiations
//
// We need to make an instance here so as to register it with measurePhotometry
//
// \cond
#define MAKE_PHOTOMETRYFINDERS(IMAGE_T)                                 \
            namespace {                                                 \
                measurePhotometry<lsst::afw::image::MaskedImage<IMAGE_T> >* foo = \
                    measureNaivePhotometry<lsst::afw::image::MaskedImage<IMAGE_T> >::getInstance(0); \
            }
            
MAKE_PHOTOMETRYFINDERS(float)


// \endcond

}}}
