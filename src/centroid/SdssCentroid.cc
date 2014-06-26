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
#include "lsst/meas/algorithms/DoubleGaussianPsf.h"
#include "lsst/afw/math/ConvolveImage.h"
#include "lsst/afw/math/offsetImage.h"
#include "lsst/afw/geom/Angle.h"
#include "lsst/meas/algorithms/PSF.h"
#include "lsst/meas/algorithms/CentroidControl.h"
#include "lsst/utils/ieee.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace afwGeom = lsst::afw::geom;
namespace afwImage = lsst::afw::image;
namespace afwMath = lsst::afw::math;

namespace lsst {
namespace meas {
namespace algorithms {
    
namespace {

/**
 * @brief A class that knows how to calculate centroids using the SDSS centroiding algorithm
 */
class SdssCentroid : public CentroidAlgorithm {
public:

    SdssCentroid(SdssCentroidControl const & ctrl, afw::table::Schema & schema) :
        CentroidAlgorithm(ctrl, schema, "SDSS-algorithm centroid measurement")
    {}

private:
    
    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(SdssCentroid);

};

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
        sVar = skyVar/(8*afwGeom::PI*sigma2)*(1 - exp(-1/sigma2));
        dVar = skyVar/(2*afwGeom::PI*sigma2)*(3 - 4*exp(-1/(4*sigma2)) + exp(-1/sigma2));

        sVar += sourceVar/(12*afwGeom::PI*sigma2)*(exp(-1/(3*sigma2)) - exp(-1/sigma2));
        dVar += sourceVar/(3*afwGeom::PI*sigma2)*(2 - 3*exp(-1/(3*sigma2)) + exp(-1/sigma2));
    }

    xVar = sVar*pow(1/d + k/(4*As)*(1 - 12*s*s/(d*d)), 2) +
        dVar*pow(s/(d*d) - k/(4*As)*8*s*s/(d*d*d), 2);

    return(xVar >= 0 ? sqrt(xVar) : NAN);
}

/************************************************************************************************************/
/*
 * Estimate the position of an object, assuming we know that it's approximately the size of the PSF 
 */
namespace {

#if 1
template<typename ImageXy_locatorT, typename VarImageXy_locatorT>
void doMeasureCentroidImpl(double *xCenter, // output; x-position of object
                       double *dxc,     // output; error in xCenter
                       double *yCenter, // output; y-position of object
                       double *dyc,     // output; error in yCenter
                       double *sizeX2, double *sizeY2, // output; object widths^2 in x and y directions
                       ImageXy_locatorT im, // Locator for the pixel values
                       VarImageXy_locatorT vim, // Locator for the image containing the variance
                       double smoothingSigma // Gaussian sigma of already-applied smoothing filter
                      )
{
    /*
     * find a first quadratic estimate
     */
    double const d2x = 2*im(0, 0) - im(-1,  0) - im(1, 0);
    double const d2y = 2*im(0, 0) - im( 0, -1) - im(0, 1);
    double const sx =  0.5*(im(1, 0) - im(-1,  0));
    double const sy =  0.5*(im(0, 1) - im( 0, -1));

    if (d2x == 0.0 || d2y == 0.0) {
        throw LSST_EXCEPT(pexExceptions::RuntimeError, "Object has a vanishing 2nd derivative");
    }
    if (d2x < 0.0 || d2y < 0.0) {
        throw LSST_EXCEPT(pexExceptions::RuntimeError,
                          (boost::format("Object is not at a maximum: d2I/dx2, d2I/dy2 = %g %g")
                           % d2x % d2y).str());
    }

    double const dx0 = sx/d2x;
    double const dy0 = sy/d2y;          // first guess

    if (fabs(dx0) > 10.0 || fabs(dy0) > 10.0) {
        throw LSST_EXCEPT(pexExceptions::RuntimeError,
                          (boost::format("Object has an almost vanishing 2nd derivative:"
                                         " sx, d2x, sy, d2y = %f %f %f %f")
                           % sx % d2x % sy % d2y).str());
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

    double xc, yc;                      // position of maximum
    double sigmaX2, sigmaY2;            // widths^2 in x and y of smoothed object

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
    float tauX2 = sigmaX2;              // width^2 of _un_ smoothed object
    float tauY2 = sigmaY2; 
    tauX2 -= smoothingSigma*smoothingSigma;              // correct for smoothing
    tauY2 -= smoothingSigma*smoothingSigma;

    if (tauX2 <= smoothingSigma*smoothingSigma) {         // problem; sigmaX2 must be bad
        tauX2 = smoothingSigma*smoothingSigma;
    }
    if (tauY2 <= smoothingSigma*smoothingSigma) {         // sigmaY2 must be bad
        tauY2 = smoothingSigma*smoothingSigma;
    }

    float const skyVar = (vim(-1, -1) + vim( 0, -1) + vim( 1, -1) + 
                          vim(-1,  0)               + vim( 1,  0) + 
                          vim(-1,  1) + vim( 0,  1) + vim( 1,  1))/8.0; // Variance in sky
    float const sourceVar = vim(0, 0);                         // extra variance of peak due to its photons
    float const A = vpk*sqrt((sigmaX2/tauX2)*(sigmaY2/tauY2)); // peak of Unsmoothed object

    *xCenter = xc;
    *yCenter = yc;

    *dxc = astrom_errors(skyVar, sourceVar, A, tauX2, vpk, sx, d2x, fabs(smoothingSigma), quarticBad);
    *dyc = astrom_errors(skyVar, sourceVar, A, tauY2, vpk, sy, d2y, fabs(smoothingSigma), quarticBad);

    *sizeX2 = tauX2;                    // return the estimates of the (object size)^2
    *sizeY2 = tauY2;
}
#endif

template<typename MaskedImageXy_locatorT>
void doMeasureCentroidImpl(double *xCenter, // output; x-position of object
                       double *dxc,     // output; error in xCenter
                       double *yCenter, // output; y-position of object
                       double *dyc,     // output; error in yCenter
                       double *sizeX2, double *sizeY2, // output; object widths^2 in x and y directions
                       double *peakVal,                // output; peak of object
                       MaskedImageXy_locatorT mim, // Locator for the pixel values
                       double smoothingSigma // Gaussian sigma of already-applied smoothing filter
                      )
{
    /*
     * find a first quadratic estimate
     */
    double const d2x = 2*mim.image(0, 0) - mim.image(-1,  0) - mim.image(1, 0);
    double const d2y = 2*mim.image(0, 0) - mim.image( 0, -1) - mim.image(0, 1);
    double const sx =  0.5*(mim.image(1, 0) - mim.image(-1,  0));
    double const sy =  0.5*(mim.image(0, 1) - mim.image( 0, -1));

    if (d2x == 0.0 || d2y == 0.0) {
        throw LSST_EXCEPT(pexExceptions::RuntimeError, "Object has a vanishing 2nd derivative");
    }
    if (d2x < 0.0 || d2y < 0.0) {
        throw LSST_EXCEPT(pexExceptions::RuntimeError,
                          (boost::format("Object is not at a maximum: d2I/dx2, d2I/dy2 = %g %g")
                           % d2x % d2y).str());
    }

    double const dx0 = sx/d2x;
    double const dy0 = sy/d2y;          // first guess

    if (fabs(dx0) > 10.0 || fabs(dy0) > 10.0) {
        throw LSST_EXCEPT(pexExceptions::RuntimeError,
                          (boost::format("Object has an almost vanishing 2nd derivative:"
                                         " sx, d2x, sy, d2y = %f %f %f %f")
                           % sx % d2x % sy % d2y).str());
    }

    double vpk = mim.image(0, 0) + 0.5*(sx*dx0 + sy*dy0); // height of peak in image
    if (vpk < 0) {
        vpk = -vpk;
    }
/*
 * now evaluate maxima on stripes
 */
    float m0x = 0, m1x = 0, m2x = 0;
    float m0y = 0, m1y = 0, m2y = 0;
    
    int quarticBad = 0;
    quarticBad += inter4(mim.image(-1, -1), mim.image( 0, -1), mim.image( 1, -1), &m0x);
    quarticBad += inter4(mim.image(-1,  0), mim.image( 0,  0), mim.image( 1,  0), &m1x);
    quarticBad += inter4(mim.image(-1,  1), mim.image( 0,  1), mim.image( 1,  1), &m2x);
   
    quarticBad += inter4(mim.image(-1, -1), mim.image(-1,  0), mim.image(-1,  1), &m0y);
    quarticBad += inter4(mim.image( 0, -1), mim.image( 0,  0), mim.image( 0,  1), &m1y);
    quarticBad += inter4(mim.image( 1, -1), mim.image( 1,  0), mim.image( 1,  1), &m2y);

    double xc, yc;                      // position of maximum
    double sigmaX2, sigmaY2;            // widths^2 in x and y of smoothed object

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
    float tauX2 = sigmaX2;              // width^2 of _un_ smoothed object
    float tauY2 = sigmaY2; 
    tauX2 -= smoothingSigma*smoothingSigma;              // correct for smoothing
    tauY2 -= smoothingSigma*smoothingSigma;

    if (tauX2 <= smoothingSigma*smoothingSigma) {         // problem; sigmaX2 must be bad
        tauX2 = smoothingSigma*smoothingSigma;
    }
    if (tauY2 <= smoothingSigma*smoothingSigma) {         // sigmaY2 must be bad
        tauY2 = smoothingSigma*smoothingSigma;
    }

    float const skyVar = (mim.variance(-1, -1) + mim.variance( 0, -1) + mim.variance( 1, -1) + 
                          mim.variance(-1,  0)                        + mim.variance( 1,  0) + 
                          mim.variance(-1,  1) + mim.variance( 0,  1) + mim.variance( 1,  1)
                         )/8.0; // Variance in sky
    float const sourceVar = mim.variance(0, 0);                // extra variance of peak due to its photons
    float const A = vpk*sqrt((sigmaX2/tauX2)*(sigmaY2/tauY2)); // peak of Unsmoothed object

    *xCenter = xc;
    *yCenter = yc;

    *dxc = astrom_errors(skyVar, sourceVar, A, tauX2, vpk, sx, d2x, fabs(smoothingSigma), quarticBad);
    *dyc = astrom_errors(skyVar, sourceVar, A, tauY2, vpk, sy, d2y, fabs(smoothingSigma), quarticBad);

    *sizeX2 = tauX2;                    // return the estimates of the (object size)^2
    *sizeY2 = tauY2;

    *peakVal = vpk;
}    

template<typename MaskedImageT>
std::pair<MaskedImageT, double>
smoothAndBinImage(CONST_PTR(lsst::afw::detection::Psf) psf,
            int const x, const int y,
            MaskedImageT const& mimage,
            int binX, int binY)
{
    afwGeom::Point2D const center(x + mimage.getX0(), y + mimage.getY0());
    afwGeom::ellipses::Quadrupole const& shape = psf->computeShape(center);
    double const smoothingSigma = shape.getDeterminantRadius();
#if 0
    double const nEffective = psf->computeEffectiveArea(); // not implemented yet (#2821)
#else
    double const nEffective = 4*M_PI*smoothingSigma*smoothingSigma; // correct for a Gaussian
#endif

    afwMath::Kernel::ConstPtr kernel = psf->getLocalKernel(center);
    int const kWidth = kernel->getWidth();
    int const kHeight = kernel->getHeight();

    afwGeom::BoxI bbox(afwGeom::Point2I(x - binX*(2 + kWidth/2), y - binY*(2 + kHeight/2)),
                       afwGeom::ExtentI(binX*(3 + kWidth + 1), binY*(3 + kHeight + 1)));
        
    // image to smooth, a shallow copy
    MaskedImageT subImage = MaskedImageT(mimage, bbox, afwImage::LOCAL);
    PTR(MaskedImageT) binnedImage = afwMath::binImage(subImage, binX, binY, lsst::afw::math::MEAN);
    binnedImage->setXY0(subImage.getXY0());
    // image to smooth into, a deep copy.  
    MaskedImageT smoothedImage = MaskedImageT(*binnedImage, true);
    assert(smoothedImage.getWidth()/2  == kWidth/2  + 2); // assumed by the code that uses smoothedImage
    assert(smoothedImage.getHeight()/2 == kHeight/2 + 2);

    afwMath::convolve(smoothedImage, *binnedImage, *kernel, afwMath::ConvolutionControl());
    *smoothedImage.getVariance() *= binX*binY*nEffective; // We want the per-pixel variance, so undo the
                                                    // effects of binning and smoothing

    return std::make_pair(smoothedImage, smoothingSigma);
}

}

/************************************************************************************************************/
/**
 * @brief Given an image and a pixel position, return a Centroid using the SDSS algorithm
 */

template<typename PixelT>
void SdssCentroid::_apply(
    afw::table::SourceRecord & source, 
    afw::image::Exposure<PixelT> const& exposure,
    afw::geom::Point2D const & center
) const {
    source.set(getKeys().meas, center); // better than NaN
    source.set(getKeys().flag, true); // say we've failed so that's the result if we throw
    typedef typename afw::image::Exposure<PixelT>::MaskedImageT MaskedImageT;
    typedef typename MaskedImageT::Image ImageT;
    typedef typename MaskedImageT::Variance VarianceT;

    MaskedImageT const& mimage = exposure.getMaskedImage();
    ImageT const& image = *mimage.getImage();
    CONST_PTR(lsst::afw::detection::Psf) psf = exposure.getPsf();

    int const x = static_cast<int>(center.getX() + 0.5) - image.getX0(); // in image Pixel coords
    int const y = static_cast<int>(center.getY() + 0.5) - image.getY0();

    if (x < 0 || x >= image.getWidth() || y < 0 || y >= image.getHeight()) {
         throw LSST_EXCEPT(lsst::pex::exceptions::LengthError,
                           (boost::format("Object at (%d, %d) is off the frame") % x % y).str());
    }
    /*
     * If a PSF is provided, smooth the object with that PSF
     */
    if (!psf) {                  // image is presumably already smoothed
        // FIXME: the above logic is probably bad; this option should probably be a config parameter
        psf.reset(new DoubleGaussianPsf(11, 11, 0.01));
    }
    
    SdssCentroidControl const & ctrl = static_cast<SdssCentroidControl const &>(getControl());

    int binX = 1;
    int binY = 1;
    double xc=0., yc=0., dxc=0., dyc=0.;            // estimated centre and error therein
    for(int binsize = 1; binsize <= ctrl.binmax; binsize *= 2) {
        std::pair<MaskedImageT, double> result = smoothAndBinImage(psf, x, y, mimage, binX, binY);
        MaskedImageT const smoothedImage = result.first;
        double const smoothingSigma = result.second;

        typename MaskedImageT::xy_locator mim = smoothedImage.xy_at(smoothedImage.getWidth()/2,
                                                                    smoothedImage.getHeight()/2);

        try {
            double sizeX2, sizeY2;      // object widths^2 in x and y directions
            double peakVal;             // peak intensity in image

            doMeasureCentroidImpl(&xc, &dxc, &yc, &dyc, &sizeX2, &sizeY2, &peakVal, mim, smoothingSigma);

            if(binsize > 1) {
                // dilate from the lower left corner of central pixel
                xc = (xc + 0.5)*binX - 0.5;
                dxc *= binX;
                sizeX2 *= binX*binX;

                yc = (yc + 0.5)*binY - 0.5;
                dyc *= binY;
                sizeY2 *= binY*binY;
            }
      
            xc += x;                    // xc, yc are measured relative to pixel (x, y)
            yc += y;

            double const fac = ctrl.wfac*(1 + smoothingSigma*smoothingSigma);
            double const facX2 = fac*binX*binX;
            double const facY2 = fac*binY*binY;

            if (sizeX2 < facX2 && ::pow(xc - x, 2) < facX2 &&
                sizeY2 < facY2 && ::pow(yc - y, 2) < facY2) {
                if (binsize > 1 || ctrl.peakMin < 0.0 || peakVal > ctrl.peakMin) {
                    break;
                }
            }

            if (sizeX2 >= facX2 || ::pow(xc - x, 2) >= facX2) {
                binX *= 2;
            }
            if (sizeY2 >= facY2 || ::pow(yc - y, 2) >= facY2) {
                binY *= 2;
            }
        } catch(pexExceptions::Exception &e) {
            LSST_EXCEPT_ADD(e, (boost::format("Object %d at (%d, %d)")
                                % source.getId() % x % y).str());
            throw e;
        }
    }

    source.set(getKeys().flag, false);
    source.set(
        getKeys().meas,
        afw::geom::Point2D(
            afwImage::indexToPosition(xc + image.getX0()),
            afwImage::indexToPosition(yc + image.getY0())
        )
    );
    // FIXME: should include off-diagonal term in covariance
    source.set(getKeys().err(0, 0), dxc*dxc);
    source.set(getKeys().err(1, 1), dyc*dyc);
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(SdssCentroid);

} // anonymous

PTR(AlgorithmControl) SdssCentroidControl::_clone() const {
    return boost::make_shared<SdssCentroidControl>(*this);
}

PTR(Algorithm) SdssCentroidControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const &
) const {
    return boost::make_shared<SdssCentroid>(*this, boost::ref(schema));
}

}}} // namespace lsst::meas::algorithms
