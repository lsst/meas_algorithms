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
namespace afwDetection = lsst::afw::detection;
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
class GaussianPhotometry : public afwDetection::Photometry
{
public:
    typedef boost::shared_ptr<GaussianPhotometry> Ptr;
    typedef boost::shared_ptr<GaussianPhotometry const> ConstPtr;

    /// Ctor
    GaussianPhotometry(double flux, double fluxErr=std::numeric_limits<double>::quiet_NaN()) :
        afwDetection::Photometry(flux, fluxErr) {}

    /// Add desired fields to the schema
    virtual void defineSchema(afwDetection::Schema::Ptr schema ///< our schema; == _mySchema
                     ) {
        Photometry::defineSchema(schema);
    }

    static bool doConfigure(lsst::pex::policy::Policy const& policy)
    {
        if (policy.isDouble("apRadius")) {
            _apRadius = policy.getDouble("apRadius");
        } 
        if (policy.isDouble("background")) {
            _background = policy.getDouble("background");
        } 
        if (policy.isDouble("shiftmax")) {
            _shiftmax = policy.getDouble("shiftmax");
        } 
        
        return true;
    }
    template<typename ImageT>
    static Photometry::Ptr doMeasure(CONST_PTR(ImageT),
                                     CONST_PTR(afwDetection::Peak),
                                     CONST_PTR(afwDetection::Source)
                                    );

private:
    static double _apRadius;
    static double _background;
    static double _shiftmax;

    GaussianPhotometry(void) : afwDetection::Photometry() { }
    LSST_SERIALIZE_PARENT(afwDetection::Photometry)
};

LSST_REGISTER_SERIALIZER(GaussianPhotometry)

double GaussianPhotometry::_apRadius = 9.0;   // Radius (in pixels) for PSF aperture correction
double GaussianPhotometry::_background = 0.0; // the frame's background level
double GaussianPhotometry::_shiftmax = 10;    // Max allowed centroid shift

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
afwDetection::Photometry::Ptr GaussianPhotometry::doMeasure(CONST_PTR(ExposureT) exposure,
                                                            CONST_PTR(afwDetection::Peak) peak,
                                                            CONST_PTR(afwDetection::Source)
                                                           )
{
    double flux = std::numeric_limits<double>::quiet_NaN();
    double fluxErr = std::numeric_limits<double>::quiet_NaN();

    if (!peak) {
        return boost::make_shared<GaussianPhotometry>(flux, fluxErr);
    }

    typedef typename ExposureT::MaskedImageT MaskedImageT;
    typedef typename MaskedImageT::Image Image;
    typedef typename Image::Pixel Pixel;
    typedef typename Image::Ptr ImagePtr;

    MaskedImageT const& mimage = exposure->getMaskedImage();
    
    double const xcen = peak->getFx() - mimage.getX0(); ///< object's column position in image pixel coords
    double const ycen = peak->getFy() - mimage.getY0();  ///< object's row position
    /*
     * Find the object's adaptive-moments.  N.b. it would be better to use the SdssShape measurement
     * as this code repeats the work of that measurement
     */
    {
        std::pair<double, double> flux_fluxErr = getGaussianFlux(mimage, _background, xcen, ycen, _shiftmax);
        flux = flux_fluxErr.first;
        fluxErr = flux_fluxErr.second;
    }
    /*
     * Calculate aperture correction
     */
    afwDetection::Psf::ConstPtr psf = exposure->getPsf();
    if (psf) {
        afwDetection::Psf::Image::Ptr psfImage; // the image of the PSF
        
        try {
#if 0
            psfImage = psf->computeImage(afwGeom::PointD(xcen, ycen));
#else  // Pad the PSF image
            {
                typedef afwDetection::Psf::Image PsfImageT;
                PsfImageT::Ptr psfImage0 = psf->computeImage(afwGeom::PointD(xcen, ycen));

                int const pad = 5;
                psfImage = PsfImageT::Ptr(
                    new PsfImageT(psfImage0->getDimensions() + afwGeom::Extent2I(2*pad))
                );
                afwGeom::BoxI middleBBox(afwGeom::Point2I(pad, pad),
                                          psfImage0->getDimensions());

                PsfImageT::Ptr middle(new PsfImageT(*psfImage, middleBBox, afwImage::LOCAL));
                *middle <<= *psfImage0;
            }
#endif
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
        /*
         * And the flux out to the Canonical radius used for aperture corrections
         */
        double psfApFlux = 0.0;            // The measured Psf flux for the psf model
        try {
            afwImage::MaskedImage<double> psfMimage(psfImage);
            psfApFlux =
                photometry::calculateSincApertureFlux(psfMimage, psfXCen, psfYCen, 0.0, _apRadius).first;
        } catch (lsst::pex::exceptions::Exception & e) {
            std::cerr << (boost::format("Failed to measure %dx%d image at (%.1f, %.1f)")
                          % psfImage->getWidth() % psfImage->getHeight()
                          % xcen % ycen).str() << std::endl;
            psfApFlux = std::accumulate(psfImage->begin(true), psfImage->end(true), psfApFlux);
        }
        /*
         * Correct the measured flux for our object so that if it's a PSF we'll
         * get the aperture corrected psf flux
         */
        double const psfGaussianFlux = psf_flux_fluxErr.first;
#if 0
        double const psfGaussianFluxErr = psf_flux_fluxErr.second; // NaN -- no variance in the psfImage
#endif

        flux *=    psfApFlux/psfGaussianFlux;
        fluxErr *= psfApFlux/psfGaussianFlux;
    }

    return boost::make_shared<GaussianPhotometry>(flux, fluxErr);
}

/*
 * Declare the existence of a "GAUSSIAN" algorithm to MeasurePhotometry
 *
 * \cond
 */
#define MAKE_PHOTOMETRYS(TYPE)                                          \
    MeasurePhotometry<afwImage::Exposure<TYPE> >::declare("GAUSSIAN", \
        &GaussianPhotometry::doMeasure<afwImage::Exposure<TYPE> >, \
        &GaussianPhotometry::doConfigure \
    )

namespace {
    volatile bool isInstance[] = {
        MAKE_PHOTOMETRYS(float)
#if 0
        ,MAKE_PHOTOMETRYS(double)
#endif
    };
}
    
// \endcond

}}}
