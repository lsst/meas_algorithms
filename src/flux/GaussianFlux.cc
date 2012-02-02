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
        ),
        _badApCorrKey(
            schema.addField<afw::table::Flag>(
                ctrl.name + ".flags.badapcorr",
                "the GaussianFlux algorithm's built-in aperture correction failed"
            )
        )
    {
        if (ctrl.fixed) {
            _centroidKey = schema[ctrl.centroid];
            _shapeKey = schema[ctrl.shape];
        }
    }

private:
    
    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(GaussianFlux);

    afw::table::Key<afw::table::Flag> _badApCorrKey;
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
        result = getGaussianFlux(mimage, ctrl.background, xcen, ycen, ctrl.shiftmax);
    }
    source.set(getKeys().meas, result.first);
    source.set(getKeys().err, result.second);
    
    /*
     * Correct the measured flux for our object so that if it's a PSF we'll
     * get the aperture corrected psf flux
     */
    source.set(_badApCorrKey, true);
    double correction = getApertureCorrection(exposure.getPsf(), xcen, ycen, ctrl.shiftmax);
    source[getKeys().meas] *=    correction;
    source[getKeys().err] *= correction;
    source.set(_badApCorrKey, false);

    source.set(getKeys().flag, false);
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(GaussianFlux);
                  
} // anonymous namespace

PTR(AlgorithmControl) GaussianFluxControl::_clone() const {
    return boost::make_shared<GaussianFluxControl>(*this);
}

PTR(Algorithm) GaussianFluxControl::_makeAlgorithm(afw::table::Schema & schema) const {
    return boost::make_shared<GaussianFlux>(*this, boost::ref(schema));
}

}}}
