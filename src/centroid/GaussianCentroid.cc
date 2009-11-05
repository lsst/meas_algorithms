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
    static bool registerMe(std::string const& name);
protected:
    friend class MeasureCentroidFactory<GaussianMeasureCentroid>;
    GaussianMeasureCentroid() : MeasureCentroid<ImageT>() {}
private:
    Centroid doApply(ImageT const& image, int x, int y, PSF const*, double) const;
};

/**
 * Register the factory that builds GaussianMeasureCentroid
 *
 * \note This function returns bool so that it can be used in an initialisation at file scope to do the actual
 * registration
 */
template<typename ImageT>
bool GaussianMeasureCentroid<ImageT>::registerMe(std::string const& name) {
    static bool _registered = false;

    if (!_registered) {
        MeasureCentroidFactory<GaussianMeasureCentroid> *factory =
            new MeasureCentroidFactory<GaussianMeasureCentroid>();
        factory->markPersistent();

        GaussianMeasureCentroid::declare(name, factory);
        _registered = true;
    }

    return true;
}

/**
 * @brief Given an image and a pixel position, return a Centroid using a naive 3x3 weighted moment
 */
template<typename ImageT>
Centroid GaussianMeasureCentroid<ImageT>::doApply(ImageT const& image, ///< The Image wherein dwells the object
                                          int x,               ///< object's column position
                                          int y,               ///< object's row position
                                          PSF const*,          ///< image's PSF
                                          double ///< image's background level
                                         ) const {
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

//
// Explicit instantiations
//
// We need to make an instance here so as to register it with MeasureCentroid
//
// \cond
#define MAKE_CENTROIDERS(IMAGE_T) \
    bool isInstance = GaussianMeasureCentroid<lsst::afw::image::Image<IMAGE_T> >::registerMe("GAUSSIAN");
                
MAKE_CENTROIDERS(float)

// \endcond

}}}}
