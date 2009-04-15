/**
 * @file
 */
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/meas/algorithms/Centroid.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace afwImage = lsst::afw::image;

namespace lsst { namespace meas { namespace algorithms { namespace{
/**
 * @brief A class that knows how to calculate centroids using the SDSS centroiding algorithm
 */
template<typename ImageT>
class SdssMeasureCentroid : public MeasureCentroid<ImageT> {
public:
    static bool registerMe(std::string const& name);
protected:
    friend class MeasureCentroidFactory<SdssMeasureCentroid>;
    SdssMeasureCentroid() : MeasureCentroid<ImageT>() {}
private:
    
    Centroid doApply(ImageT const& image, int x, int y, PSF const*, double background) const;
};

/**
 * Register the factory that builds SdssMeasureCentroid
 *
 * \note This function returns bool so that it can be used in an initialisation at file scope to do the actual
 * registration
 */
template<typename ImageT>
bool SdssMeasureCentroid<ImageT>::registerMe(std::string const& name) {
    static bool _registered = false;

    if (!_registered) {
        MeasureCentroidFactory<SdssMeasureCentroid> *factory = new MeasureCentroidFactory<SdssMeasureCentroid>();
        factory->markPersistent();

        SdssMeasureCentroid::declare(name, factory);
        _registered = true;
    }

    return true;
}

/************************************************************************************************************/

float const AMPAST4 = 1.33;           // amplitude of `4th order' corr compared to theory

/* 
 * Do the Gaussian quartic interpolation for the position
 * of the maximum for three equally spaced values vm,v0,vp, assumed to be at
 * abscissae -1,0,1; the answer is returned as *cen
 *
 * Return 0 is all is well, otherwise 1
 */
static int
    inter4(float vm, float v0, float vp, float *cen)
{
    float const sp = v0 - vp;
    float const sm = v0 - vm;
    float const d2 = sp + sm;
    float const s = 0.5*(vp - vm);

    if(d2 <= 0.0f || v0 <= 0.0f) {
        return(1);
    }
    
    *cen = s/d2*(1. + AMPAST4*sp*sm/(d2*v0));

    return fabs(*cen) < 1 ? 0 : 1;
}

/*****************************************************************************/
/*
 * Calculate error in centroid
 */
float astrom_errors(float gain,		// CCD's gain
                    float esky,		// noise-equivalent sky, including background variance
                    float A,                // abs(peak value in raw image)
                    float tau2,		// Object is N(0,tau2)
                    float As,               // abs(peak value in smoothed image)
                    float s,                // slope across central pixel
                    float d,                // curvature at central pixel
                    float sigma,		// width of smoothing filter
                    int quartic_bad)        // was quartic estimate bad?
{
    const float k = quartic_bad ? 0 : AMPAST4; /* quartic correction coeff */
    const float sigma2 = sigma*sigma;	/* == sigma^2 */
    float sVar, dVar;			/* variances of s and d */
    float xVar;				/* variance of centroid, x */

    if(As == 0.0 || d == 0.0) {
        return(1e3);
    }

    if(sigma == 0) {			/* no smoothing; no covariance */
        sVar = esky/gain/2;                 /* sky */
        dVar = 6*esky/gain;

        sVar += 0.5*(A/gain)*exp(-1/(2*tau2));
        dVar += (A/gain)*(4*exp(-1/(2*tau2)) + 2*exp(-1/(2*tau2)));
    } else {				/* smoothed */
        sVar = esky/gain/(8*M_PI*sigma2)*(1 - exp(-1/sigma2));
        dVar = esky/gain/(2*M_PI*sigma2)*
            (3 - 4*exp(-1/(4*sigma2)) + exp(-1/sigma2));

        sVar += (A/gain)/(12*M_PI*sigma2)*(exp(-1/(3*sigma2)) - exp(-1/sigma2));
        dVar += (A/gain)/(3*M_PI*sigma2)*
            (2 - 3*exp(-1/(3*sigma2)) + exp(-1/sigma2));
    }

    xVar = sVar*pow(1/d + k/(4*As)*(1 - 12*s*s/(d*d)), 2) +
        dVar*pow(s/(d*d) - k/(4*As)*8*s*s/(d*d*d), 2);

    return(xVar >= 0 ? sqrt(xVar) : NAN);
}

/**
 * @brief Given an image and a pixel position, return a Centroid using the SDSS algorithm
 */
template<typename ImageT>
Centroid SdssMeasureCentroid<ImageT>::doApply(ImageT const& image, ///< The Image wherein dwells the object
                                         int x,               ///< object's column position
                                         int y,               ///< object's row position
                                         PSF const* psf,      ///< image's PSF (NULL if image is already smoothed)
                                         double background    ///< image's background level
                                        ) const {
    x -= image.getX0();                 // work in image Pixel coordinates
    y -= image.getY0();
    /*
     * If a PSF is provided, smooth the object with that PSF
     */
    typename ImageT::xy_locator im;                    // locator for the (possible smoothed) image
    typename ImageT::template ImageTypeFactory<>::type tmp(3, 3); // a (small piece of the) smoothed image, if needed

    if (psf == NULL) {                  // image is presumably already smoothed
        im = image.xy_at(x, y);
    } else {
        int const kWidth = psf->getKernel()->getWidth();
        int const kHeight = psf->getKernel()->getHeight();

        afwImage::BBox bbox(afwImage::PointI(x - 2 - kWidth/2, y - 2 - kHeight/2), 3 + kWidth + 1, 3 + kHeight + 1);

        ImageT subImage = ImageT(image, bbox);     // image to smooth, a shallow copy
        ImageT smoothedImage = ImageT(image, bbox, true); // image to smooth into, a deep copy.  Forgets [XY]0
        psf->convolve(smoothedImage, subImage);

        tmp <<= ImageT(smoothedImage, afwImage::BBox(afwImage::PointI(1 + kWidth/2, 1 + kHeight/2), 3, 3));
        im = tmp.xy_at(1, 1);
    }
    /*
     * find a first quadratic estimate
     */
    double const d2x = 2*im(0, 0) - im(-1,  0) - im(1, 0);
    double const d2y = 2*im(0, 0) - im( 0, -1) - im(0, 1);
    double const sx =  0.5*(im(1, 0) - im(-1,  0));
    double const sy =  0.5*(im(0, 1) - im( 0, -1));

    if(d2x == 0.0 || d2y == 0.0) {
        throw LSST_EXCEPT(pexExceptions::RuntimeErrorException,
                          (boost::format("Object at (%d, %d) has a vanishing 2nd derivative") % x % y).str());
    }
    if(d2x < 0.0 || d2y < 0.0) {
        throw LSST_EXCEPT(pexExceptions::RuntimeErrorException,
                          (boost::format("Object at (%d, %d) is not a maximum: d2I/dx2, d2I/dy2 = %g %f")
                           % x % y % d2x % d2y).str());
    }

    double const dx0 = sx/d2x;
    double const dy0 = sy/d2y;          // first guess
   
    double vpk = im(0, 0) + 0.5*(sx*dx0 + sy*dy0); // height of peak in image
    int sign_peak = 1;
    if(vpk < 0) {
        sign_peak = -1;
        vpk = -vpk;
    }
/*
 * now evaluate maxima on stripes
 */
    float m0x = 0, m1x = 0, m2x = 0;
    float m0y = 0, m1y = 0, m2y = 0;
    
    int quartic_bad = 0;
    quartic_bad += inter4(im(-1, -1), im( 0, -1), im( 1, -1), &m0x);
    quartic_bad += inter4(im(-1,  0), im( 0,  0), im( 1,  0), &m1x);
    quartic_bad += inter4(im(-1,  1), im( 0,  1), im( 1,  1), &m2x);
   
    quartic_bad += inter4(im(-1, -1), im(-1,  0), im(-1,  1), &m0y);
    quartic_bad += inter4(im( 0, -1), im( 0,  0), im( 0,  1), &m1y);
    quartic_bad += inter4(im( 1, -1), im( 1,  0), im( 1,  1), &m2y);

    double xc, yc;                              // position of maximum
    double sigma_x2, sigma_y2;                  // widths^2 in x and y

    if(quartic_bad) {			// >= 1 quartic interpolator is bad
        xc = dx0;
        yc = dy0;
        sigma_x2 = vpk/d2x;		// widths^2 in x
        sigma_y2 = vpk/d2y;		//             and y
   } else {
        double const smx = 0.5*(m2x-m0x);
        double const smy = 0.5*(m2y-m0y);
        double const dm2x = m1x - 0.5*(m0x+m2x);
        double const dm2y = m1y - 0.5*(m0y+m2y);      
        double const dx = m1x + dy0*(smx - dy0*dm2x); // first quartic approx
        double const dy = m1y + dx0*(smy - dx0*dm2y);
        double const dx4 = m1x + dy*(smx - dy*dm2x);	// second quartic approx
        double const dy4 = m1y + dx*(smy - dx*dm2y);
      
        xc = dx4;
        yc = dy4;
        sigma_x2 = vpk/d2x - (1 + 6*dx0*dx0)/4; // widths^2 in x
        sigma_y2 = vpk/d2y - (1 + 6*dy0*dy0)/4; //             and y
    }
    /*
     * Now for the errors.
     */
    float const gain = 1; float const dark_variance = 0; // XXX
    float const sigma = 0;		// sigma of Gaussian to smooth image with (if < 0, assume that the
                                        // region is already smoothed)

    float esky = background + gain*dark_variance; // noise-equivalent sky, including background variance
      
    float tau_x2 = sigma_x2;            // width^2 of _un_ smoothed object
    float tau_y2 = sigma_y2; 
    tau_x2 -= sigma*sigma;		// correct for smoothing
    tau_y2 -= sigma*sigma;

    if(tau_x2 <= sigma*sigma) {         // problem; sigma_x2 must be bad
        tau_x2 = sigma*sigma;
    }
    if(tau_y2 <= sigma*sigma) {         // sigma_y2 must be bad
        tau_y2 = sigma*sigma;
    }

    float const A = vpk*sqrt((sigma_x2/tau_x2)*(sigma_y2/tau_y2)); // peak of Unsmoothed object

    double const dxc = astrom_errors(gain, esky, A, tau_x2, vpk, sx, d2x, fabs(sigma), quartic_bad);
    double const dyc = astrom_errors(gain, esky, A, tau_y2, vpk, sy, d2y, fabs(sigma), quartic_bad);

    return Centroid(Centroid::xyAndError(afwImage::indexToPosition(x + image.getX0()) + xc, dxc),
                    Centroid::xyAndError(afwImage::indexToPosition(y + image.getY0()) + yc, dyc));
}

//
// Explicit instantiations
//
// We need to call registerMe here to register SdssMeasureCentroid; doRegister returns bool solely to allow us
// to assign it to a file-static variable, and thus register at programme startup
//
// \cond
#define MAKE_CENTROIDERS(IMAGE_T) \
    bool b = SdssMeasureCentroid<afwImage::Image<IMAGE_T> >::registerMe("SDSS");
                
MAKE_CENTROIDERS(float)

// \endcond

}}}}
