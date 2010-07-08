#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/meas/algorithms/Centroid.h"
#include "all.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;

namespace lsst { namespace meas { namespace algorithms {

namespace {
/**
 * @brief A class that knows how to calculate centroids using a 2-D Gaussian fitter from Dave Monet
 */
template<typename ImageT>
class GaussianMeasureCentroid : public MeasureCentroid<ImageT> {
public:
    typedef MeasureCentroid<ImageT> MeasurePropertyBase;
    
    GaussianMeasureCentroid(typename ImageT::ConstPtr image) : MeasureCentroid<ImageT>(image) {}
private:
    Centroid doApply(ImageT const& image, int x, int y, lsst::afw::detection::Psf const*, double) const;
};

/**
 * @brief Given an image and a pixel position, return a Centroid using a naive 3x3 weighted moment
 */
template<typename ImageT>
Centroid GaussianMeasureCentroid<ImageT>::doApply(ImageT const& image, ///< The Image wherein dwells the object
                                                  int x,               ///< object's column position
                                                  int y,               ///< object's row position
                                                  lsst::afw::detection::Psf const*,          ///< image's PSF
                                                  double ///< image's background level
                                                 ) const
{
    x -= image.getX0();                 // work in image Pixel coordinates
    y -= image.getY0();

#if 0
    if (sum == 0.0) {
        throw LSST_EXCEPT(pexExceptions::RuntimeErrorException,
                          (boost::format("Object at (%d, %d) has no counts") % x % y).str());
    }
#endif

    FittedModel fit = twodg(image, x, y); // here's the fitter

    return Centroid(lsst::afw::image::indexToPosition(image.getX0()) + fit.params[FittedModel::X0],
                    lsst::afw::image::indexToPosition(image.getY0()) + fit.params[FittedModel::Y0]);
}
//}
//
// Explicit instantiations
//
// We need to make an instance here so as to register it with MeasureCentroid
//
// \cond
#define MAKE_CENTROIDERS(IMAGE_T) \
    registerMe<GaussianMeasureCentroid, lsst::afw::image::Image<IMAGE_T> >("GAUSSIAN")
                
    volatile bool isInstance[] = {
        MAKE_CENTROIDERS(int),
        MAKE_CENTROIDERS(float)
    };
// \endcond

}}}}
