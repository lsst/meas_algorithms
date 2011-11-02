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
#include "lsst/meas/algorithms/Measure.h"

#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/detection/Photometry.h"
#include "lsst/meas/algorithms/detail/SdssShape.h"
#include "lsst/meas/algorithms/Photometry.h"

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
class GaussianPhotometer : public Algorithm<afwDet::Photometry, ExposureT>
{
public:
    typedef Algorithm<afwDet::Photometry, ExposureT> AlgorithmT;
    typedef boost::shared_ptr<GaussianPhotometer> Ptr;
    typedef boost::shared_ptr<GaussianPhotometer const> ConstPtr;

    /// Ctor
    GaussianPhotometer(bool fixed=false, double shiftmax=10.0, double background=0.0) :
        AlgorithmT(), _fixed(fixed), _shiftmax(shiftmax), _background(background) {}

    virtual std::string getName() const { return "GAUSSIAN"; }

    virtual PTR(AlgorithmT) clone() const {
        return boost::make_shared<GaussianPhotometer<ExposureT> >(_fixed, _shiftmax, _background);
    }

    virtual void configure(lsst::pex::policy::Policy const& policy) {
        if (policy.isBool("fixed")) {
            _fixed = policy.getBool("fixed");
        }
        if (policy.isDouble("background")) {
            _background = policy.getDouble("background");
        } 
        if (policy.isDouble("shiftmax")) {
            _shiftmax = policy.getDouble("shiftmax");
        }
    }

    virtual PTR(afwDet::Photometry) measureSingle(afwDet::Source const&, afwDet::Source const&,
                                                  ExposurePatch<ExposureT> const&) const;

private:
    double _shiftmax;
    double _background;
    bool _fixed;
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
afwDet::Photometry::Ptr GaussianPhotometer<ExposureT>::measureSingle(
    afwDet::Source const& target,
    afwDet::Source const& source,
    ExposurePatch<ExposureT> const& patch
    ) const
{
    CONST_PTR(ExposureT) exposure = patch.getExposure();
    typename ExposureT::MaskedImageT const& mimage = exposure->getMaskedImage();

    double const xcen = patch.getCenter().getX() - mimage.getX0(); ///< column position in image pixel coords
    double const ycen = patch.getCenter().getY() - mimage.getY0(); ///< row position

    std::pair<double, double> result;
    if (_fixed) {
        // Fixed aperture, defined by SDSS shape measurement made elsewhere
        PTR(afwDet::Shape) shape = source.getShape()->find("SDSS");
        if (!shape) {
            throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException,
                              "No SDSS shape found in source.");
        }
        afwGeom::AffineTransform const& transform = patch.fromStandard();
        detail::SdssShapeImpl sdss(*shape);
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

    return boost::make_shared<afwDet::Photometry>(flux, fluxErr);
}

// Declare the existence of a "GAUSSIAN" algorithm to MeasurePhotometry
LSST_DECLARE_ALGORITHM(GaussianPhotometer, afwDet::Photometry);

}}}
