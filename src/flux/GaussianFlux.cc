// -*- LSST-C++ -*-
#include <numeric>
#include <cmath>
#include <functional>
#include "boost/make_shared.hpp"
#include "boost/tuple/tuple.hpp"
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/geom/Angle.h"
#include "lsst/afw/image.h"
#include "lsst/afw/math/Integrate.h"
#include "lsst/afw/coord/Coord.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/meas/algorithms/detail/SdssShape.h"
#include "lsst/meas/algorithms/FluxControl.h"
#include "lsst/meas/algorithms/GaussianFluxControl.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace afwDet = lsst::afw::detection;
namespace afwGeom = lsst::afw::geom;
namespace afwImage = lsst::afw::image;
namespace afwMath = lsst::afw::math;
namespace afwCoord = lsst::afw::coord;

namespace lsst {
namespace meas {
namespace algorithms {

namespace {

/**
 * @brief A class that knows how to calculate fluxes using the GAUSSIAN photometry algorithm
 * @ingroup meas/algorithms
 */
class GaussianFlux : public FluxAlgorithm {
public:

    GaussianFlux(GaussianFluxControl const & ctrl, afw::table::Schema & schema) :
        FluxAlgorithm(
            ctrl, schema,
            "linear fit to an elliptical Gaussian with shape parameters set by adaptive moments"
        )
    {
        if (ctrl.fixed) {
            _centroidKey = schema[ctrl.centroid];
            _shapeKey = schema[ctrl.shape];
            _shapeFlagKey = schema[ctrl.shape + ctrl.shapeFlag];
        }
    }

private:

    /// Record measurement and measure PSF correction
    ///
    /// Common code for both _apply and _applyForced.
    /// Most arguments are the same as for those.
    template <typename PixelT>
    void _measurement(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const& exposure,
        afw::geom::Point2D const& center,
        std::pair<double, double> const& measurement ///< Measurement (flux, error) to record
        ) const;

    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    template <typename PixelT>
    void _applyForced(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center,
        afw::table::SourceRecord const & reference,
        afw::geom::AffineTransform const & refToMeas
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(GaussianFlux);

    afw::table::Centroid::MeasKey _centroidKey;
    afw::table::Shape::MeasKey _shapeKey;
    afw::table::Key<afw::table::Flag> _shapeFlagKey;
};

/************************************************************************************************************/
    
template<typename ImageT>
std::pair<double, double>
getGaussianFlux(
    ImageT const& mimage,           // the data to process
    double background,               // background level
    double xcen, double ycen,         // centre of object
    double shiftmax,                  // max allowed centroid shift
    int maxIter=detail::SDSS_SHAPE_MAX_ITER, ///< Maximum number of iterations
    float tol1=detail::SDSS_SHAPE_TOL1, ///< Convergence tolerance for e1,e2
    float tol2=detail::SDSS_SHAPE_TOL2, ///< Convergence tolerance for FWHM
    PTR(detail::SdssShapeImpl) shape=PTR(detail::SdssShapeImpl)() // SDSS shape measurement
) {
    double flux = std::numeric_limits<double>::quiet_NaN();
    double fluxErr = std::numeric_limits<double>::quiet_NaN();

    if (!shape) {
        shape = boost::make_shared<detail::SdssShapeImpl>();
    }

    if (!detail::getAdaptiveMoments(mimage, background, xcen, ycen, shiftmax, shape.get(),
                                    maxIter, tol1, tol2)) {
        ;                               // Should set a flag here
    } else {
        double const scale = shape->getFluxScale();
        flux = scale*shape->getI0();
        fluxErr = scale*shape->getI0Err();
    }

    return std::make_pair(flux, fluxErr);
}

/************************************************************************************************************/
/**
 * Calculate the desired gaussian flux
 */

template <typename PixelT>
void GaussianFlux::_measurement(
    afw::table::SourceRecord& source,
    afw::image::Exposure<PixelT> const& exposure,
    afw::geom::Point2D const& center,
    std::pair<double, double> const& measurement
) const {
    source.set(getKeys().meas, measurement.first);
    source.set(getKeys().err, measurement.second);
    source.set(getKeys().flag, false);
}


template <typename PixelT>
void GaussianFlux::_apply(
    afw::table::SourceRecord & source, 
    afw::image::Exposure<PixelT> const& exposure,
    afw::geom::Point2D const & center
) const {
    source.set(getKeys().flag, true); // say we've failed so that's the result if we throw
    typename afw::image::Exposure<PixelT>::MaskedImageT const& mimage = exposure.getMaskedImage();

    double const xcen = center.getX() - mimage.getX0(); ///< column position in image pixel coords
    double const ycen = center.getY() - mimage.getY0(); ///< row position

    GaussianFluxControl const & ctrl = static_cast<GaussianFluxControl const &>(getControl());

    std::pair<double, double> result;
    if (ctrl.fixed) {
        // Fixed aperture, defined by SDSS shape measurement made elsewhere
        if (source.get(_shapeFlagKey)) {
            throw LSST_EXCEPT(pexExceptions::RuntimeErrorException, "Shape measurement failed");
        }
        detail::SdssShapeImpl sdss(source.get(_centroidKey), source.get(_shapeKey));
        result = detail::getFixedMomentsFlux(mimage, ctrl.background, xcen, ycen, sdss);
    } else {
        // FIXME: propagate SDSS shape measurement flags.
        /*
         * Find the object's adaptive-moments.  N.b. it would be better to use the SdssShape measurement
         * as this code repeats the work of that measurement
         */
        result = getGaussianFlux(mimage, ctrl.background, xcen, ycen, ctrl.shiftmax, ctrl.maxIter,
                                 ctrl.tol1, ctrl.tol2);
    }

    _measurement(source, exposure, center, result);
}

template <typename PixelT>
void GaussianFlux::_applyForced(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center,
    afw::table::SourceRecord const & reference,
    afw::geom::AffineTransform const & refToMeas
    ) const
{
    source.set(getKeys().flag, true); // bad unless we get all the way to success at the end
    GaussianFluxControl const& ctrl = static_cast<GaussianFluxControl const &>(this->getControl());
    int const x0 = exposure.getX0(), y0 = exposure.getY0();
    // Fixed aperture, defined by SDSS shape measurement on the reference
    afw::geom::ellipses::Quadrupole const& refShape =
        reference.get(reference.getSchema().find<afw::table::Shape::MeasTag>(ctrl.shape).key);
    detail::SdssShapeImpl sdss(center, refShape);
    std::pair<double, double> const& result =
        detail::getFixedMomentsFlux(exposure.getMaskedImage(), ctrl.background,
                                    center.getX() - x0, center.getY() - y0, sdss);
    _measurement(source, exposure, center, result);
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(GaussianFlux);

} // anonymous namespace

PTR(AlgorithmControl) GaussianFluxControl::_clone() const {
    return boost::make_shared<GaussianFluxControl>(*this);
}

PTR(Algorithm) GaussianFluxControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const &
) const {
    return boost::make_shared<GaussianFlux>(*this, boost::ref(schema));
}

}}}
