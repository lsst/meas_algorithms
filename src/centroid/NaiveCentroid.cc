// -*- LSST-C++ -*-
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/meas/algorithms/Measure.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace afwDetection = lsst::afw::detection;
namespace afwImage = lsst::afw::image;

namespace lsst {
namespace meas {
namespace algorithms {

namespace {

/**
 * @brief A class that knows how to calculate centroids as a simple unweighted first moment
 * of the 3x3 region around a pixel
 */
class NaiveAstrometry : public afwDetection::Astrometry
{
public:
    typedef boost::shared_ptr<NaiveAstrometry> Ptr;
    typedef boost::shared_ptr<NaiveAstrometry const> ConstPtr;

    /// Ctor
    NaiveAstrometry(double x, double xErr, double y, double yErr)
    {
        init();                         // This allocates space for fields added by defineSchema
        set<X>(x);                      // ... if you don't, these set calls will fail an assertion
        set<X_ERR>(xErr);               // the type of the value must match the schema
        set<Y>(y);
        set<Y_ERR>(yErr);
    }

    /// Add desired fields to the schema
    virtual void defineSchema(afwDetection::Schema::Ptr schema ///< our schema; == _mySchema
                     ) {
        Astrometry::defineSchema(schema);
    }

    static bool doConfigure(lsst::pex::policy::Policy const& policy);

    template<typename MaskedImageT>
    static Astrometry::Ptr doMeasure(typename MaskedImageT::ConstPtr im, afwDetection::Peak const&);
};

/**
 * @brief Given an image and a pixel position, return a Centroid using a naive 3x3 weighted moment
 */
/**
 * Process the image; calculate values
 */
template<typename MaskedImageT>
afwDetection::Astrometry::Ptr NaiveAstrometry::doMeasure(typename MaskedImageT::ConstPtr image,
                                                         afwDetection::Peak const& peak)
{
    int x = static_cast<int>(peak.getIx() + 0.5);
    int y = static_cast<int>(peak.getIy() + 0.5);

    x -= image->getX0();                 // work in image Pixel coordinates
    y -= image->getY0();

    typename MaskedImageT::Image::xy_locator im = image->getImage()->xy_at(x, y);

    float const background = 0.0;

    double const sum =
        (im(-1,  1) + im( 0,  1) + im( 1,  1) +
         im(-1,  0) + im( 0,  0) + im( 1,  0) +
         im(-1, -1) + im( 0, -1) + im( 1, -1)) - 9*background;

    if (sum == 0.0) {
        throw LSST_EXCEPT(pexExceptions::RuntimeErrorException,
                          (boost::format("Object at (%d, %d) has no counts") %
                           peak.getIx() % peak.getIy()).str());
    }

    double const sum_x =
        -im(-1,  1) + im( 1,  1) +
        -im(-1,  0) + im( 1,  0) +
        -im(-1, -1) + im( 1, -1);
    double const sum_y =
        (im(-1,  1) + im( 0,  1) + im( 1,  1)) -
        (im(-1, -1) + im( 0, -1) + im( 1, -1));

    return boost::make_shared<NaiveAstrometry>(
        lsst::afw::image::indexToPosition(x + image->getX0()) + sum_x/sum, 0.0,
        lsst::afw::image::indexToPosition(y + image->getY0()) + sum_y/sum, 0.0);
}

//
// Explicit instantiations
//
// We need to make an instance here so as to register it with MeasureCentroid
//
// \cond
/**
 * Declare the existence of a "NAIVE" algorithm
 */
#define INSTANTIATE(TYPE) \
    NewMeasureAstrometry<afwImage::MaskedImage<TYPE> >::declare("NAIVE", \
        &NaiveAstrometry::doMeasure<afwImage::MaskedImage<TYPE> > \
        )

volatile bool isInstance[] = {
    INSTANTIATE(int),
    INSTANTIATE(float),
    INSTANTIATE(double)
};

// \endcond

}}}}
