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
        MeasureCentroidFactory<NaiveMeasureCentroid> *factory = new MeasureCentroidFactory<NaiveMeasureCentroid>();
        factory->markPersistent();

        NaiveMeasureCentroid::declare(name, factory);
        _registered = true;
    }

    return true;
}

/**
 * @brief Given an image and a pixel position, return a Centroid using a naive 3x3 weighted moment
 */
template<typename ImageT>
Centroid NaiveMeasureCentroid<ImageT>::doApply(ImageT const& image, ///< The Image wherein dwells the object
                                          int x,               ///< object's column position
                                          int y,               ///< object's row position
                                          PSF const*,          ///< image's PSF
                                          double background    ///< image's background level
                                         ) const {
    x -= image.getX0();                 // work in image Pixel coordinates
    y -= image.getY0();

    typename ImageT::xy_locator im = image.xy_at(x, y);

    double const sum =
        (im(-1,  1) + im( 0,  1) + im( 1,  1) +
         im(-1,  0) + im( 0,  0) + im( 1,  0) +
         im(-1, -1) + im( 0, -1) + im( 1, -1)) - 9*background;

    if (sum == 0.0) {
        throw LSST_EXCEPT(pexExceptions::RuntimeErrorException,
                          (boost::format("Object at (%d, %d) has no counts") % x % y).str());
    }

    double const sum_x =
        -im(-1,  1) + im( 1,  1) +
        -im(-1,  0) + im( 1,  0) +
        -im(-1, -1) + im( 1, -1);
    double const sum_y =
        (im(-1,  1) + im( 0,  1) + im( 1,  1)) -
        (im(-1, -1) + im( 0, -1) + im( 1, -1));

    return Centroid(lsst::afw::image::indexToPosition(x + image.getX0()) + sum_x/sum,
                    lsst::afw::image::indexToPosition(y + image.getY0()) + sum_y/sum);
}

//
// Explicit instantiations
//
// We need to make an instance here so as to register it with MeasureCentroid
//
// \cond
#define MAKE_CENTROIDERS(IMAGE_T) \
    bool b = NaiveMeasureCentroid<lsst::afw::image::Image<IMAGE_T> >::registerMe("NAIVE");
                
MAKE_CENTROIDERS(float)

// \endcond

}}}}
