// -*- LSST-C++ -*-
#include <numeric>
#include <cmath>
#include <functional>
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/image.h"
#include "lsst/afw/math/Integrate.h"
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
    GaussianPhotometer(double apRadius=9.0, double shiftmax=10.0, double background=0.0) : 
        AlgorithmT(), _apRadius(apRadius), _shiftmax(shiftmax), _background(background) {}

    virtual std::string getName() const { return "GAUSSIAN"; }

    virtual PTR(AlgorithmT) clone() const {
        return boost::make_shared<GaussianPhotometer<ExposureT> >(_apRadius, _shiftmax, _background);
    }

    virtual PTR(afwDet::Photometry) measureNull(void) const {
        const double NaN = std::numeric_limits<double>::quiet_NaN();
        return boost::make_shared<afwDet::Photometry>(NaN, NaN);
    }

    virtual void configure(lsst::pex::policy::Policy const& policy) {
        if (policy.isDouble("apRadius")) {
            _apRadius = policy.getDouble("apRadius");
        } 
        if (policy.isDouble("background")) {
            _background = policy.getDouble("background");
        } 
        if (policy.isDouble("shiftmax")) {
            _shiftmax = policy.getDouble("shiftmax");
        } 
    }

    virtual PTR(afwDet::Photometry) measureOne(ExposurePatch<ExposureT> const&, afwDet::Source const&) const;

private:
    double _apRadius;
    double _shiftmax;
    double _background;
};

/************************************************************************************************************/

namespace {
    
template<typename ImageT>
std::pair<double, double>
getGaussianFlux(
        ImageT const& mimage,           // the data to process
        double background,               // background level
        double xcen, double ycen,         // centre of object
        double shiftmax                  // max allowed centroid shift
               )
{
    double flux = std::numeric_limits<double>::quiet_NaN();
    double fluxErr = std::numeric_limits<double>::quiet_NaN();

    detail::SdssShapeImpl shapeImpl;

    if (!detail::getAdaptiveMoments(mimage, background, xcen, ycen, shiftmax, &shapeImpl)) {
        ;                               // Should set a flag here
    } else {
        /*
         * The shape is an ellipse that's axis-aligned in (u, v) [<uv> = 0] after rotation by theta:
         * <x^2> + <y^2> = <u^2> + <v^2>
         * <x^2> - <y^2> = cos(2 theta)*(<u^2> - <v^2>)
         * 2*<xy>        = sin(2 theta)*(<u^2> - <v^2>)
         */
        double const Mxx = shapeImpl.getIxx(); // <x^2>
        double const Mxy = shapeImpl.getIxy(); // <xy>
        double const Myy = shapeImpl.getIyy(); // <y^2>
        
        double const Muu_p_Mvv = Mxx + Myy;                             // <u^2> + <v^2>
        double const Muu_m_Mvv = ::sqrt(::pow(Mxx - Myy, 2) + 4*::pow(Mxy, 2)); // <u^2> - <v^2>
        double const Muu = 0.5*(Muu_p_Mvv + Muu_m_Mvv);
        double const Mvv = 0.5*(Muu_p_Mvv - Muu_m_Mvv);
        
        double const scale = 2*M_PI*::sqrt(Muu*Mvv);
        flux = scale*shapeImpl.getI0();
        fluxErr = scale*shapeImpl.getI0Err();
    }

    return std::make_pair(flux, fluxErr);
}
}

/************************************************************************************************************/
/**
 * Calculate the desired gaussian flux
 */
template<typename ExposureT>
afwDet::Photometry::Ptr GaussianPhotometer<ExposureT>::measureOne(ExposurePatch<ExposureT> const& patch,
                                                                  afwDet::Source const& source) const
{
    typedef typename ExposureT::MaskedImageT MaskedImageT;
    typedef typename MaskedImageT::Image Image;
    typedef typename Image::Pixel Pixel;
    typedef typename Image::Ptr ImagePtr;

    CONST_PTR(ExposureT) exposure = patch.getExposure();
    CONST_PTR(afwDet::Peak) peak = patch.getPeak();

    MaskedImageT const& mimage = exposure->getMaskedImage();
    
    double const xcen = peak->getFx() - mimage.getX0(); ///< object's column position in image pixel coords
    double const ycen = peak->getFy() - mimage.getY0();  ///< object's row position
    /*
     * Find the object's adaptive-moments.  N.b. it would be better to use the SdssShape measurement
     * as this code repeats the work of that measurement
     */
    double flux, fluxErr;
    {
        std::pair<double, double> flux_fluxErr = getGaussianFlux(mimage, _background, xcen, ycen, _shiftmax);
        flux = flux_fluxErr.first;
        fluxErr = flux_fluxErr.second;
    }

    /*
     * Calculate aperture correction
     */
    afwDet::Psf::ConstPtr psf = exposure->getPsf();
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
    std::pair<double, double> psf_flux_fluxErr =
        getGaussianFlux(*psfImage, 0.0, psfXCen, psfYCen, _shiftmax);
    double const psfGaussianFlux = psf_flux_fluxErr.first;
#if 0
    double const psfGaussianFluxErr = psf_flux_fluxErr.second; // NaN -- no variance in the psfImage
#endif

    // Need to correct to the PSF flux
    double psfFlux = std::accumulate(psfImageNoPad->begin(true), psfImageNoPad->end(true), 0.0);

    /*
     * Correct the measured flux for our object so that if it's a PSF we'll
     * get the aperture corrected psf flux
     */
    double correction = psfFlux/psfGaussianFlux;
    flux *=    correction;
    fluxErr *= correction;

    return boost::make_shared<afwDet::Photometry>(flux, fluxErr);
}

// Declare the existence of a "GAUSSIAN" algorithm to MeasurePhotometry
DECLARE_ALGORITHM(GaussianPhotometer, afwDet::Photometry);

}}}
