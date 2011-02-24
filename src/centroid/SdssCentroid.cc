// -*- LSST-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
/**
 * @file
 */
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/math/ConvolveImage.h"
#include "lsst/meas/algorithms/Measure.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace afwDetection = lsst::afw::detection;
namespace afwGeom = lsst::afw::geom;
namespace afwImage = lsst::afw::image;
namespace afwMath = lsst::afw::math;

namespace lsst {
namespace meas {
namespace algorithms {
    
namespace{
/**
 * @brief A class that knows how to calculate centroids using the SDSS centroiding algorithm
 */
class SdssAstrometry : public afwDetection::Astrometry
{
public:
    typedef boost::shared_ptr<SdssAstrometry> Ptr;
    typedef boost::shared_ptr<SdssAstrometry const> ConstPtr;

    /// Ctor
    SdssAstrometry(double x, double xErr, double y, double yErr) : afwDetection::Astrometry(x, xErr, y, yErr) {}

    /// Add desired fields to the schema
    virtual void defineSchema(afwDetection::Schema::Ptr schema ///< our schema; == _mySchema
                     ) {
        Astrometry::defineSchema(schema);
    }

    static bool doConfigure(lsst::pex::policy::Policy const& policy) { return true; }

    template<typename ExposureT>
    static Astrometry::Ptr doMeasure(typename ExposureT::ConstPtr im, afwDetection::Peak const*);

private:
    SdssAstrometry(void) : afwDetection::Astrometry() { }
    LSST_SERIALIZE_PARENT(afwDetection::Astrometry)
};

LSST_REGISTER_SERIALIZER(SdssAstrometry)


/************************************************************************************************************/

float const AMPAST4 = 1.33;           // amplitude of `4th order' corr compared to theory

/* 
 * Do the Gaussian quartic interpolation for the position
 * of the maximum for three equally spaced values vm,v0,vp, assumed to be at
 * abscissae -1,0,1; the answer is returned as *cen
 *
 * Return 0 is all is well, otherwise 1
 */
static int inter4(float vm, float v0, float vp, float *cen) {
    float const sp = v0 - vp;
    float const sm = v0 - vm;
    float const d2 = sp + sm;
    float const s = 0.5*(vp - vm);

    if (d2 <= 0.0f || v0 <= 0.0f) {
        return(1);
    }
    
    *cen = s/d2*(1.0 + AMPAST4*sp*sm/(d2*v0));
    
    return fabs(*cen) < 1 ? 0 : 1;
}

/*****************************************************************************/
/*
 * Calculate error in centroid
 */
float astrom_errors(float skyVar,       // variance of pixels at the sky level
                    float sourceVar,    // variance in peak due to excess counts over sky
                    float A,            // abs(peak value in raw image)
                    float tau2,         // Object is N(0,tau2)
                    float As,           // abs(peak value in smoothed image)
                    float s,            // slope across central pixel
                    float d,            // curvature at central pixel
                    float sigma,        // width of smoothing filter
                    int quarticBad) {   // was quartic estimate bad?

    float const k = quarticBad ? 0 : AMPAST4; /* quartic correction coeff */
    float const sigma2 = sigma*sigma;   /* == sigma^2 */
    float sVar, dVar;                   /* variances of s and d */
    float xVar;                         /* variance of centroid, x */
    
    if (fabs(As) < std::numeric_limits<float>::min() ||
        fabs(d)  < std::numeric_limits<float>::min()) {
        return(1e3);
    }

    if (sigma <= 0) {                    /* no smoothing; no covariance */
        sVar = 0.5*skyVar;               /* due to sky */
        dVar = 6*skyVar;

        sVar += 0.5*sourceVar*exp(-1/(2*tau2));
        dVar += sourceVar*(4*exp(-1/(2*tau2)) + 2*exp(-1/(2*tau2)));
    } else {                            /* smoothed */
        sVar = skyVar/(8*M_PI*sigma2)*(1 - exp(-1/sigma2));
        dVar = skyVar/(2*M_PI*sigma2)*(3 - 4*exp(-1/(4*sigma2)) + exp(-1/sigma2));

        sVar += sourceVar/(12*M_PI*sigma2)*(exp(-1/(3*sigma2)) - exp(-1/sigma2));
        dVar += sourceVar/(3*M_PI*sigma2)*(2 - 3*exp(-1/(3*sigma2)) + exp(-1/sigma2));
    }

    xVar = sVar*pow(1/d + k/(4*As)*(1 - 12*s*s/(d*d)), 2) +
        dVar*pow(s/(d*d) - k/(4*As)*8*s*s/(d*d*d), 2);

    return(xVar >= 0 ? sqrt(xVar) : NAN);
}

/**
 * @brief Given an image and a pixel position, return a Centroid using the SDSS algorithm
 */
template<typename ExposureT>
afwDetection::Astrometry::Ptr SdssAstrometry::doMeasure(typename ExposureT::ConstPtr exposure,
                                                        afwDetection::Peak const* peak)
{
    if (!peak) {
        double const pos = std::numeric_limits<double>::quiet_NaN();
        double const posErr = std::numeric_limits<double>::quiet_NaN();
        return boost::make_shared<SdssAstrometry>(pos, posErr, pos, posErr);
    }

    typedef typename ExposureT::MaskedImageT MaskedImageT;
    typedef typename MaskedImageT::Image ImageT;
    MaskedImageT const& mimage = exposure->getMaskedImage();
    ImageT const& image = *mimage.getImage();
    lsst::afw::detection::Psf::ConstPtr psf = exposure->getPsf();

    int x = static_cast<int>(peak->getIx() + 0.5);
    int y = static_cast<int>(peak->getIy() + 0.5);

    x -= image.getX0();                 // work in image Pixel coordinates
    y -= image.getY0();

    if (x < 0 || x >= image.getWidth() || y < 0 || y >= image.getHeight()) {
         throw LSST_EXCEPT(lsst::pex::exceptions::LengthErrorException,
                           (boost::format("Object at (%d, %d) is off the frame") % x % y).str());
    }
    /*
     * If a PSF is provided, smooth the object with that PSF
     */
    typename ImageT::xy_locator im;                               // locator for the (possible smoothed) image
    typename ImageT::template ImageTypeFactory<>::type tmp(3, 3); // a (small piece of the) smoothed image

    if (psf == NULL) {                  // image is presumably already smoothed
        im = image.xy_at(x, y);
    } else {
        afwMath::Kernel::ConstPtr kernel = psf->getLocalKernel(afwGeom::makePointD(x, y));
        int const kWidth = kernel->getWidth();
        int const kHeight = kernel->getHeight();

        afwImage::BBox bbox(afwImage::PointI(x - 2 - kWidth/2, y - 2 - kHeight/2),
                            3 + kWidth + 1, 3 + kHeight + 1);
        
        ImageT subImage = ImageT(image, bbox);     // image to smooth, a shallow copy
        ImageT smoothedImage = ImageT(image, bbox, true); // image to smooth into, a deep copy.  Forgets [XY]0

        afwMath::convolve(smoothedImage, subImage, *kernel, afwMath::ConvolutionControl());
        
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

    if (d2x == 0.0 || d2y == 0.0) {
        throw LSST_EXCEPT(pexExceptions::RuntimeErrorException,
                          (boost::format("Object at (%d, %d) has a vanishing 2nd derivative") % x % y).str());
    }
    if (d2x < 0.0 || d2y < 0.0) {
        throw LSST_EXCEPT(pexExceptions::RuntimeErrorException,
                          (boost::format("Object at (%d, %d) is not a maximum: d2I/dx2, d2I/dy2 = %g %f")
                           % x % y % d2x % d2y).str());
    }

    double const dx0 = sx/d2x;
    double const dy0 = sy/d2y;          // first guess

    if (fabs(dx0) > 10.0 || fabs(dy0) > 10.0) {
        throw LSST_EXCEPT(pexExceptions::RuntimeErrorException,
                          (boost::format("Object at (%d, %d) has an almost vanishing 2nd derivative:"
                                         " sx, d2x, sy, d2y = %f %f %f %f")
                           % x % y % sx % d2x % sy % d2y).str());
    }

    double vpk = im(0, 0) + 0.5*(sx*dx0 + sy*dy0); // height of peak in image
    if (vpk < 0) {
        vpk = -vpk;
    }
/*
 * now evaluate maxima on stripes
 */
    float m0x = 0, m1x = 0, m2x = 0;
    float m0y = 0, m1y = 0, m2y = 0;
    
    int quarticBad = 0;
    quarticBad += inter4(im(-1, -1), im( 0, -1), im( 1, -1), &m0x);
    quarticBad += inter4(im(-1,  0), im( 0,  0), im( 1,  0), &m1x);
    quarticBad += inter4(im(-1,  1), im( 0,  1), im( 1,  1), &m2x);
   
    quarticBad += inter4(im(-1, -1), im(-1,  0), im(-1,  1), &m0y);
    quarticBad += inter4(im( 0, -1), im( 0,  0), im( 0,  1), &m1y);
    quarticBad += inter4(im( 1, -1), im( 1,  0), im( 1,  1), &m2y);

    double xc, yc;                              // position of maximum
    double sigmaX2, sigmaY2;                    // widths^2 in x and y

    if (quarticBad) {                   // >= 1 quartic interpolator is bad
        xc = dx0;
        yc = dy0;
        sigmaX2 = vpk/d2x;             // widths^2 in x
        sigmaY2 = vpk/d2y;             //             and y
   } else {
        double const smx = 0.5*(m2x - m0x);
        double const smy = 0.5*(m2y - m0y);
        double const dm2x = m1x - 0.5*(m0x + m2x);
        double const dm2y = m1y - 0.5*(m0y + m2y);      
        double const dx = m1x + dy0*(smx - dy0*dm2x); // first quartic approx
        double const dy = m1y + dx0*(smy - dx0*dm2y);
        double const dx4 = m1x + dy*(smx - dy*dm2x);    // second quartic approx
        double const dy4 = m1y + dx*(smy - dx*dm2y);
      
        xc = dx4;
        yc = dy4;
        sigmaX2 = vpk/d2x - (1 + 6*dx0*dx0)/4; // widths^2 in x
        sigmaY2 = vpk/d2y - (1 + 6*dy0*dy0)/4; //             and y
    }
    /*
     * Now for the errors.
     */
    float const sigma = 0;              // sigma of Gaussian to smooth image with (if < 0, assume that the
                                        // region is already smoothed)

    float tauX2 = sigmaX2;            // width^2 of _un_ smoothed object
    float tauY2 = sigmaY2; 
    tauX2 -= sigma*sigma;              // correct for smoothing
    tauY2 -= sigma*sigma;

    if (tauX2 <= sigma*sigma) {         // problem; sigmaX2 must be bad
        tauX2 = sigma*sigma;
    }
    if (tauY2 <= sigma*sigma) {         // sigmaY2 must be bad
        tauY2 = sigma*sigma;
    }

    float const skyVar = (mimage.xy_at(x - 1, y - 1).variance() + 
                          mimage.xy_at(x    , y - 1).variance() + 
                          mimage.xy_at(x + 1, y - 1).variance() + 
                          mimage.xy_at(x - 1, y    ).variance() + 
                          mimage.xy_at(x    , y    ).variance() + 
                          mimage.xy_at(x + 1, y    ).variance() + 
                          mimage.xy_at(x - 1, y + 1).variance() + 
                          mimage.xy_at(x    , y + 1).variance() + 
                          mimage.xy_at(x + 1, y + 1).variance()
                         )/8.0;         // Variance in sky, including background variance
    float const A = vpk*sqrt((sigmaX2/tauX2)*(sigmaY2/tauY2)); // peak of Unsmoothed object
    float const sourceVar = mimage.xy_at(x, y).variance(); // extra variance of peak due to its photons

    double const dxc = astrom_errors(skyVar, sourceVar, A, tauX2, vpk, sx, d2x, fabs(sigma), quarticBad);
    double const dyc = astrom_errors(skyVar, sourceVar, A, tauY2, vpk, sy, d2y, fabs(sigma), quarticBad);

    double const xCenter = afwImage::indexToPosition(x + image.getX0()) + xc;
    double const yCenter = afwImage::indexToPosition(y + image.getY0()) + yc;

    return boost::make_shared<SdssAstrometry>(xCenter, dxc, yCenter, dyc);
}

/*
 * Declare the existence of a "SDSS" algorithm to MeasureAstrometry
 *
 * \cond
 */
#define INSTANTIATE(TYPE) \
    MeasureAstrometry<afwImage::Exposure<TYPE> >::declare("SDSS", \
        &SdssAstrometry::doMeasure<afwImage::Exposure<TYPE> >,       \
        &SdssAstrometry::doConfigure \
        )

volatile bool isInstance[] = {
    INSTANTIATE(int),
    INSTANTIATE(float),
    INSTANTIATE(double)
};

// \endcond

}}}}
