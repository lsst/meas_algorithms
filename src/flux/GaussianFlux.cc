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

/**
 * @brief A class that knows how to calculate fluxes using the GAUSSIAN photometry algorithm
 * @ingroup meas/algorithms
 */
template<typename ExposureT>
class GaussianFlux : public Algorithm<ExposureT>
{
public:
    typedef Algorithm<ExposureT> AlgorithmT;

    GaussianFlux(GaussianFluxControl const & ctrl, afw::table::Schema & schema) :
        AlgorithmT(ctrl), _fixed(ctrl.fixed), _shiftmax(ctrl.shiftmax), _background(ctrl.background),
        _keys(
            addFluxFields(
                schema,
                ctrl.name,
                "linear fit to an elliptical Gaussian with shape parameters set by adaptive moments"
            )
        )
    {
        if (_fixed) {
            _centroidKey = schema[ctrl.centroid];
            _shapeKey = schema[ctrl.shape];
        }
    }

    virtual void apply(afw::table::SourceRecord &, ExposurePatch<ExposureT> const&) const;

private:
    bool _fixed;
    double _shiftmax;
    double _background;
    afw::table::KeyTuple<afw::table::Flux> _keys;
    afw::table::Key< afw::table::Point<double> > _centroidKey;
    afw::table::Key< afw::table::Moments<double> > _shapeKey;
};

/************************************************************************************************************/

namespace {
    
template<typename ImageT>
std::pair<double, double>
getGaussianFlux(
        ImageT const& mimage,           // the data to process
        double background,               // background level
        double xcen, double ycen,         // centre of object
        double shiftmax,                  // max allowed centroid shift
        PTR(detail::SdssShapeImpl) shape=PTR(detail::SdssShapeImpl)() // SDSS shape measurement
               )
{
    double flux = std::numeric_limits<double>::quiet_NaN();
    double fluxErr = std::numeric_limits<double>::quiet_NaN();

    if (!shape) {
        shape = boost::make_shared<detail::SdssShapeImpl>();
    }

    if (!detail::getAdaptiveMoments(mimage, background, xcen, ycen, shiftmax, shape.get())) {
        ;                               // Should set a flag here
    } else {
        double const scale = shape->getFluxScale();
        flux = scale*shape->getI0();
        fluxErr = scale*shape->getI0Err();
    }

    return std::make_pair(flux, fluxErr);
}


/*
 * Calculate aperture correction
 *
 * The multiplier returned will correct the measured flux for an object so that if it's a PSF we'll
 * get the aperture corrected psf flux
 */
double getApertureCorrection(afwDet::Psf::ConstPtr psf, double xcen, double ycen, double shiftmax)
{
    if (!psf) {
        throw LSST_EXCEPT(pexExceptions::RuntimeErrorException, "No PSF provided for Gaussian photometry");
    }

    typedef afwDet::Psf::Image PsfImageT;
    PsfImageT::Ptr psfImage; // the image of the PSF
    PsfImageT::Ptr psfImageNoPad;   // Unpadded image of PSF
    
    int const pad = 5;
    try {
        psfImageNoPad = psf->computeImage(afwGeom::PointD(xcen, ycen));
        
        psfImage = PsfImageT::Ptr(
            new PsfImageT(psfImageNoPad->getDimensions() + afwGeom::Extent2I(2*pad))
            );
        afwGeom::BoxI middleBBox(afwGeom::Point2I(pad, pad),
                                 psfImageNoPad->getDimensions());
        
        PsfImageT::Ptr middle(new PsfImageT(*psfImage, middleBBox, afwImage::LOCAL));
        *middle <<= *psfImageNoPad;
        psfImage->setXY0(0, 0);     // SHOULD NOT BE NEEDED; psfXCen should be 0.0 and fix getGaussianFlux
    } catch (lsst::pex::exceptions::Exception & e) {
        LSST_EXCEPT_ADD(e, (boost::format("Computing PSF at (%.3f, %.3f)") % xcen % ycen).str());
        throw e;
    }
    // Estimate the GaussianFlux for the Psf
    double const psfXCen = 0.5*(psfImage->getWidth() - 1); // Center of (21x21) image is (10.0, 10.0)
    double const psfYCen = 0.5*(psfImage->getHeight() - 1);
    std::pair<double, double> const result = getGaussianFlux(*psfImage, 0.0, psfXCen, psfYCen, shiftmax);
    double const psfGaussianFlux = result.first;
#if 0
    double const psfGaussianFluxErr = result.second; // NaN -- no variance in the psfImage
#endif

    // Need to correct to the PSF flux
    double psfFlux = std::accumulate(psfImageNoPad->begin(), psfImageNoPad->end(), 0.0);
    return psfFlux/psfGaussianFlux;
}
                  
} // anonymous namespace

/************************************************************************************************************/
/**
 * Calculate the desired gaussian flux
 */
template<typename ExposureT>
void GaussianFlux<ExposureT>::apply(
    afw::table::SourceRecord & source,
    ExposurePatch<ExposureT> const & patch
) const {
    CONST_PTR(ExposureT) exposure = patch.getExposure();
    typename ExposureT::MaskedImageT const& mimage = exposure->getMaskedImage();

    double const xcen = patch.getCenter().getX() - mimage.getX0(); ///< column position in image pixel coords
    double const ycen = patch.getCenter().getY() - mimage.getY0(); ///< row position

    std::pair<double, double> result;
    if (_fixed) {
        // Fixed aperture, defined by SDSS shape measurement made elsewhere
        afw::geom::AffineTransform const& transform = patch.fromStandard();
        detail::SdssShapeImpl sdss(source.get(_centroidKey), source.get(_shapeKey));
        result = detail::getFixedMomentsFlux(mimage, _background, xcen, ycen, *sdss.transform(transform));
    } else {
        /*
         * Find the object's adaptive-moments.  N.b. it would be better to use the SdssShape measurement
         * as this code repeats the work of that measurement
         */
        result = getGaussianFlux(mimage, _background, xcen, ycen, _shiftmax);
    }
    double flux = result.first;
    double fluxErr = result.second;

    /*
     * Correct the measured flux for our object so that if it's a PSF we'll
     * get the aperture corrected psf flux
     */
    double correction = getApertureCorrection(exposure->getPsf(), xcen, ycen, _shiftmax);
    flux *=    correction;
    fluxErr *= correction;

    source.set(_keys.meas, flux);
    source.set(_keys.err, fluxErr);
    source.set(_keys.flag, true);
}

LSST_ALGORITHM_CONTROL_PRIVATE_IMPL(GaussianFluxControl, GaussianFlux)

}}}
