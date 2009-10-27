#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/meas/algorithms/Centroid.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;

namespace lsst { namespace meas { namespace algorithms {

namespace {
/**
 * @brief A class that knows how to calculate centroids
 */
template<typename ImageT>
class SillyMeasureCentroid : public MeasureCentroid<ImageT> {
public:
    static bool registerMe(std::string const& name);
protected:
    friend class MeasureCentroidFactory<SillyMeasureCentroid>;
    SillyMeasureCentroid() : MeasureCentroid<ImageT>() {}
private:
    Centroid doApply(ImageT const& image, int x, int y, PSF const* psf, double background) const;
};

/**
 * Register the factory that builds SillyMeasureCentroid
 *
 * \note This function returns bool so that it can be used in an initialisation at file scope to do the actual
 * registration
 */
template<typename ImageT>
bool SillyMeasureCentroid<ImageT>::registerMe(std::string const& name) {
    static bool _registered = false;

    if (!_registered) {
        MeasureCentroidFactory<SillyMeasureCentroid> *factory =
            new MeasureCentroidFactory<SillyMeasureCentroid>();
        factory->markPersistent();

        SillyMeasureCentroid::declare(name, factory);
        _registered = true;
    }

    return true;
}

/**
 * @brief Given an image and a pixel position, return a Centroid offset by (1, 1) from initial position
 */
template<typename ImageT>
Centroid SillyMeasureCentroid<ImageT>::doApply(ImageT const& image, ///< The Image wherein dwells the object
                                          int x,               ///< object's column position
                                          int y,               ///< object's row position
                                          PSF const*,          ///< image's PSF
                                          double background    ///< image's background level
                                         ) const {

    return Centroid(lsst::afw::image::indexToPosition(x) + 1.0,
                    lsst::afw::image::indexToPosition(y) + 1.0);
}

//
// Explicit instantiations
//
// We need to make an instance here so as to register it with MeasureCentroid
//
// \cond
#define MAKE_CENTROIDERS(IMAGE_T) \
    bool isInstance = SillyMeasureCentroid<lsst::afw::image::Image<IMAGE_T> >::registerMe("SILLY");
                
MAKE_CENTROIDERS(float)

// \endcond

}}}}
