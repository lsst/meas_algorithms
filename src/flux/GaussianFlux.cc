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
#include "lsst/meas/algorithms/ScaledFlux.h"
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
class GaussianFlux : public FluxAlgorithm, public ScaledFlux {
public:

    GaussianFlux(GaussianFluxControl const & ctrl, afw::table::Schema & schema) :
        FluxAlgorithm(
            ctrl, schema,
            "linear fit to an elliptical Gaussian with shape parameters set by adaptive moments"        
        ),
        _fluxCorrectionKeys(ctrl.name, schema)
    {
        if (ctrl.fixed) {
            _centroidKey = schema[ctrl.centroid];
            _shapeKey = schema[ctrl.shape];
        }
    }

    virtual afw::table::KeyTuple<afw::table::Flux> getFluxKeys(int n=0) const {
        return FluxAlgorithm::getKeys();
    }

    virtual ScaledFlux::KeyTuple getFluxCorrectionKeys(int n=0) const {
        return _fluxCorrectionKeys;
    }

private:
    
    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(GaussianFlux);

    ScaledFlux::KeyTuple _fluxCorrectionKeys;
    afw::table::Centroid::MeasKey _centroidKey;
    afw::table::Shape::MeasKey _shapeKey;
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




/*
 * Apply the algorithm to the PSF model
 */
double getPsfFactor(afwDet::Psf const & psf, afw::geom::Point2D const & center, double shiftmax,
                    int maxIter=detail::SDSS_SHAPE_MAX_ITER, float tol1=detail::SDSS_SHAPE_TOL1,
                    float tol2=detail::SDSS_SHAPE_TOL2) {

    typedef afwDet::Psf::Image PsfImageT;
    PTR(PsfImageT) psfImage; // the image of the PSF
    PTR(PsfImageT) psfImageNoPad;   // Unpadded image of PSF
    
    int const pad = 5;
    try {
        psfImageNoPad = psf.computeImage(center);
        
        psfImage = PTR(PsfImageT)(
            new PsfImageT(psfImageNoPad->getDimensions() + afwGeom::Extent2I(2*pad))
            );
        afwGeom::BoxI middleBBox(afwGeom::Point2I(pad, pad), psfImageNoPad->getDimensions());
        
        PTR(PsfImageT) middle(new PsfImageT(*psfImage, middleBBox, afwImage::LOCAL));
        *middle <<= *psfImageNoPad;
    } catch (lsst::pex::exceptions::Exception & e) {
        LSST_EXCEPT_ADD(e, (boost::format("Computing PSF at (%.3f, %.3f)")
                            % center.getX() % center.getY()).str());
        throw e;
    }
    // Estimate the GaussianFlux for the Psf
    double const psfXCen = 0.5*(psfImage->getWidth() - 1); // Center of (21x21) image is (10.0, 10.0)
    double const psfYCen = 0.5*(psfImage->getHeight() - 1);
    std::pair<double, double> const result = getGaussianFlux(*psfImage, 0.0, psfXCen, psfYCen, shiftmax,
                                                             maxIter, tol1, tol2);
    return result.first;
}

/************************************************************************************************************/
/**
 * Calculate the desired gaussian flux
 */

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

    source.set(getKeys().meas, result.first);
    source.set(getKeys().err, result.second);
    source.set(getKeys().flag, false); // If PSF factor fails, we'll set this flag in the correction stage,
                                       // but we clear it now in case we don't care about corrections.

    source.set(_fluxCorrectionKeys.psfFactorFlag, true);
    if (!exposure.hasPsf()) {
        throw LSST_EXCEPT(pexExceptions::RuntimeErrorException, "No PSF provided for Gaussian photometry");
    }
    double psfFactor = getPsfFactor(*exposure.getPsf(), center, ctrl.shiftmax,
                                    ctrl.maxIter, ctrl.tol1, ctrl.tol2);
    source.set(_fluxCorrectionKeys.psfFactor, psfFactor);
    source.set(_fluxCorrectionKeys.psfFactorFlag, false);
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
