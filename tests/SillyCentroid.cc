#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/meas/algorithms/Centroid.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;

namespace lsst { namespace meas { namespace algorithms {

namespace {
/**
 * @brief A class that knows how to calculate centroids as a simple unweighted first moment
 * of the 3x3 region around a pixel
 */
template<typename ImageT>
class NaiveMeasureCentroid : public MeasureCentroid<ImageT> {
public:
    static bool registerMe(std::string const& name);
protected:
    friend class MeasureCentroidFactory<NaiveMeasureCentroid>;
    NaiveMeasureCentroid() : MeasureCentroid<ImageT>() {}
private:
    Centroid doApply(ImageT const& image, int x, int y, PSF const* psf, double background) const;
};

/**
 * Register the factory that builds NaiveMeasureCentroid
 *
 * \note This function returns bool so that it can be used in an initialisation at file scope to do the actual
 * registration
 */
template<typename ImageT>
bool NaiveMeasureCentroid<ImageT>::registerMe(std::string const& name) {
    static bool _registered = false;

    if (!_registered) {
        MeasureCentroidFactory<NaiveMeasureCentroid> *factory =
            new MeasureCentroidFactory<NaiveMeasureCentroid>();
        factory->markPersistent();

        NaiveMeasureCentroid::declare(name, factory);
        _registered = true;
    }

    return true;
}

/**
 * @brief Given an image and a pixel position, return a Centroid offset by (1, 1) from initial position
 */
template<typename ImageT>
Centroid NaiveMeasureCentroid<ImageT>::doApply(ImageT const& image, ///< The Image wherein dwells the object
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
    bool isInstance = NaiveMeasureCentroid<lsst::afw::image::Image<IMAGE_T> >::registerMe("SILLY");
                
MAKE_CENTROIDERS(float)

// \endcond

}}}}
