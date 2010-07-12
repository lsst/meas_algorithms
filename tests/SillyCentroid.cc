// -*- LSST-C++ -*-
/**
 * @file
 */
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
 * @brief A class that knows how to calculate centroids by guessing the wrong answer
 */
class SillyAstrometry : public afwDetection::Astrometry
{
public:
    typedef boost::shared_ptr<SillyAstrometry> Ptr;
    typedef boost::shared_ptr<SillyAstrometry const> ConstPtr;

    /// Ctor
    SillyAstrometry(double x, double xErr, double y, double yErr)
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

    template<typename MaskedImageT>
    static Astrometry::Ptr doMeasure(typename MaskedImageT::ConstPtr im, afwDetection::Peak const&);
};

/**
 * @brief Given an image and a pixel position, return a Centroid offset by (1, 1) from initial position
 */
template<typename MaskedImageT>
afwDetection::Astrometry::Ptr SillyAstrometry::doMeasure(typename MaskedImageT::ConstPtr image,
                                                         afwDetection::Peak const& peak)
{
    double const posErr = std::numeric_limits<double>::quiet_NaN();
    return boost::make_shared<SillyAstrometry>(peak.getFx() + 1.0, posErr,
                                               peak.getFy() + 1.0, posErr);
}

/*
 * Declare the existence of a "SILLY" algorithm to MeasureAstrometry
 *
 * \cond
 */
#define INSTANTIATE(TYPE) \
    NewMeasureAstrometry<afwImage::MaskedImage<TYPE> >::declare("SILLY", \
        &SillyAstrometry::doMeasure<afwImage::MaskedImage<TYPE> > \
        )

volatile bool isInstance[] = {
    INSTANTIATE(float),
    INSTANTIATE(double)
};

// \endcond

}}}}
