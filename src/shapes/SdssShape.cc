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
 * Measure adaptive moments.
 *
 * Originally provided by Phil Fischer, based on code from Tim McKay's group.  Error calculations by Dave
 * Johnston.  Major reworking by RHL for SDSS, and now a major rewrite for LSST
 */
#include "Eigen/LU"
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/geom/Angle.h"
#include "lsst/meas/algorithms/Measure.h"
#include "lsst/meas/algorithms/detail/SdssShape.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace afwDetection = lsst::afw::detection;
namespace afwImage = lsst::afw::image;
namespace afwGeom = lsst::afw::geom;

namespace lsst {
namespace meas {
namespace algorithms {

namespace {
    int const MAXIT = 100;              // \todo from Policy XXX
#if 0
    double const TOL1 = 0.001;             // \todo from Policy XXX
    double const TOL2 = 0.01;              // \todo from Policy XXX
#else                                   // testing
    double const TOL1 = 0.00001;           // \todo from Policy XXX
    double const TOL2 = 0.0001;            // \todo from Policy XXX
#endif

/*****************************************************************************/
/*
 * Error analysis, courtesy of David Johnston, University of Chicago
 */
/*
 * This function takes the 4 Gaussian parameters A, sigmaXXW and the
 * sky variance and fills in the Fisher matrix from the least squares fit.
 *
 * Following "Numerical Recipes in C" section 15.5, it ignores the 2nd
 * derivative parts and so the fisher matrix is just a function of these
 * best fit model parameters. The components are calculated analytically.
 */
detail::SdssShapeImpl::Matrix4
calc_fisher(detail::SdssShapeImpl const& shape, // the Shape that we want the the Fisher matrix for
            float bkgd_var              // background variance level for object
           )
{
    float const A = shape.getI0();     // amplitude
    float const sigma11W = shape.getIxx();
    float const sigma12W = shape.getIxy();
    float const sigma22W = shape.getIyy();
    
    double const D = sigma11W*sigma22W - sigma12W*sigma12W;
   
    if (D <= std::numeric_limits<double>::epsilon()) {
        throw LSST_EXCEPT(lsst::pex::exceptions::DomainErrorException,
                          "Determinant is too small calculating Fisher matrix");
    }
/*
 * a normalization factor
 */
    if (bkgd_var <= 0.0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::DomainErrorException,
                          (boost::format("Background variance must be positive (saw %g)") % bkgd_var).str());
    }
    double const F = afwGeom::PI*sqrt(D)/bkgd_var;
/*
 * Calculate the 10 independent elements of the 4x4 Fisher matrix 
 */
    detail::SdssShapeImpl::Matrix4 fisher;

    double fac = F*A/(4.0*D);
    fisher(0, 0) =  F;
    fisher(0, 1) =  fac*sigma22W;
    fisher(1, 0) =  fisher(0, 1);
    fisher(0, 2) =  fac*sigma11W;                      
    fisher(2, 0) =  fisher(0, 2);
    fisher(0, 3) = -fac*2*sigma12W;    
    fisher(3, 0) =  fisher(0, 3);
    
    fac = 3.0*F*A*A/(16.0*D*D);
    fisher(1, 1) =  fac*sigma22W*sigma22W;
    fisher(2, 2) =  fac*sigma11W*sigma11W;
    fisher(3, 3) =  fac*4.0*(sigma12W*sigma12W + D/3.0);
    
    fisher(1, 2) =  fisher(3, 3)/4.0;
    fisher(2, 1) =  fisher(1, 2);
    fisher(1, 3) =  fac*(-2*sigma22W*sigma12W);
    fisher(3, 1) =  fisher(1, 3);
    fisher(2, 3) =  fac*(-2*sigma11W*sigma12W);
    fisher(3, 2) =  fisher(2, 3);
    
    return fisher;
}
//
// Here's a class to allow us to get the Image and variance from an Image or MaskedImage
//
template<typename ImageT>               // general case
struct ImageAdaptor {
    typedef ImageT Image;

    Image const& getImage(ImageT const& image) const {
        return image;
    }

    double getVariance(ImageT const&, int, int) {
        return std::numeric_limits<double>::quiet_NaN();
    }
};
    
template<typename T>                    // specialise to a MaskedImage
struct ImageAdaptor<afwImage::MaskedImage<T> > {
    typedef typename afwImage::MaskedImage<T>::Image Image;

    Image const& getImage(afwImage::MaskedImage<T> const& mimage) const {
        return *mimage.getImage();
    }

    double getVariance(afwImage::MaskedImage<T> const& mimage, int ix, int iy) {
        return mimage.at(ix, iy).variance();
    }
};

}

/************************************************************************************************************/

namespace detail {
/**
 * Workhorse for adaptive moments
 */
template<typename ImageT>
bool
getAdaptiveMoments(ImageT const& mimage, ///< the data to process
            double bkgd,                 ///< background level
            double xcen,                 ///< x-centre of object
            double ycen,                 ///< y-centre of object
            double shiftmax,             ///< max allowed centroid shift
            detail::SdssShapeImpl *shape ///< a place to store desired data
           )
{
    float ampW = 0;                     // amplitude of best-fit Gaussian
    double sum;                         // sum of intensity*weight
    double sumx, sumy;                  // sum ((int)[xy])*intensity*weight
    double sumxx, sumxy, sumyy;         // sum {x^2,xy,y^2}*intensity*weight
    double sums4;                       // sum intensity*weight*exponent^2
    float const xcen0 = xcen;           // initial centre
    float const ycen0 = ycen;           //                of object

    double sigma11W = 1.5;              // quadratic moments of the
    double sigma12W = 0.0;              //     weighting fcn;
    double sigma22W = 1.5;              //               xx, xy, and yy

    double w11 = -1, w12 = -1, w22 = -1;        // current weights for moments; always set when iter == 0
    float e1_old = 1e6, e2_old = 1e6;           // old values of shape parameters e1 and e2
    float sigma11_ow_old = 1e6;                 // previous version of sigma11_ow
    
    typename ImageAdaptor<ImageT>::Image const &image = ImageAdaptor<ImageT>().getImage(mimage);

    if (lsst::utils::isnan(xcen) || lsst::utils::isnan(ycen)) {
        // Can't do anything
        shape->setFlags(shape->getFlags() | Flags::SHAPE_UNWEIGHTED_BAD);
        return false;
    }

    bool interpflag = false;            // interpolate finer than a pixel?
    lsst::afw::geom::BoxI bbox;
    int iter = 0;                       // iteration number
    for (; iter < MAXIT; iter++) {
        bbox = set_amom_bbox(image.getWidth(), image.getHeight(),
                             xcen, ycen, sigma11W, sigma12W, sigma22W);

        double const detW = sigma11W*sigma22W - sigma12W*sigma12W; // determinant of sigmaXXW matrix
        if (lsst::utils::isnan(detW) ||
            detW < std::numeric_limits<float>::epsilon()) { // a suitably small number
            shape->setFlags(shape->getFlags() | Flags::SHAPE_UNWEIGHTED);
            break;
        }

#if 0                                   // this form was numerically unstable on my G4 powerbook
        assert(detW >= 0.0);
#else
        assert(sigma11W*sigma22W >= sigma12W*sigma12W - std::numeric_limits<float>::epsilon());
#endif

        {
            const double ow11 = w11;    // old
            const double ow12 = w12;    //     values
            const double ow22 = w22;    //            of w11, w12, w22

            w11 = sigma22W/detW;
            w12 = -sigma12W/detW;
            w22 = sigma11W/detW;

            float const xinterp = 0.25; // I.e. 0.5*0.5
            if (sigma11W < xinterp || sigma22W < xinterp || detW < xinterp*xinterp) {
                if (!interpflag) {
                    interpflag = true;       // N.b.: stays set for this object
                    if (iter > 0) {
                        sigma11_ow_old = 1.e6; // force at least one more iteration
                        w11 = ow11;
                        w12 = ow12;
                        w22 = ow22;
                        iter--;         // we didn't update wXX
                    }
                }
            }
        }

        if (calcmom(image, xcen, ycen, bbox, bkgd, interpflag, w11, w12, w22,
                    &sum, &sumx, &sumy, &sumxx, &sumxy, &sumyy, &sums4) < 0) {
            shape->setFlags(shape->getFlags() | Flags::SHAPE_UNWEIGHTED);
            break;
        }

        ampW = sum/(afwGeom::PI*sqrt(detW));
#if 0
/*
 * Find new centre
 *
 * This is only needed if we update the centre; if we use the input position we've already done the work
 */
        xcen = sumx/sum;
        ycen = sumy/sum;
#endif
        shape->setX(sumx/sum); // update centroid.  N.b. we're not setting errors here
        shape->setY(sumy/sum);

        if (fabs(shape->getX() - xcen0) > shiftmax || fabs(shape->getY() - ycen0) > shiftmax) {
            shape->setFlags(shape->getFlags() | Flags::SHAPE_SHIFT);
        }
/*
 * OK, we have the centre. Proceed to find the second moments.
 */
        float const sigma11_ow = sumxx/sum; // quadratic moments of
        float const sigma22_ow = sumyy/sum; //          weight*object
        float const sigma12_ow = sumxy/sum; //                 xx, xy, and yy 

        if (sigma11_ow <= 0 || sigma22_ow <= 0) {
            shape->setFlags(shape->getFlags() | Flags::SHAPE_UNWEIGHTED);
            break;
        }

        float const d = sigma11_ow + sigma22_ow; // current values of shape parameters
        float const e1 = (sigma11_ow - sigma22_ow)/d;
        float const e2 = 2.0*sigma12_ow/d;
/*
 * Did we converge?
 */
        if (iter > 0 &&
           fabs(e1 - e1_old) < TOL1 && fabs(e2 - e2_old) < TOL1 &&
           fabs(sigma11_ow/sigma11_ow_old - 1.0) < TOL2 ) {
            break;                              // yes; we converged
        }

        e1_old = e1;
        e2_old = e2;
        sigma11_ow_old = sigma11_ow;
/*
 * Didn't converge, calculate new values for weighting function
 *
 * The product of two Gaussians is a Gaussian:
 * <x^2 exp(-a x^2 - 2bxy - cy^2) exp(-Ax^2 - 2Bxy - Cy^2)> = 
 *                            <x^2 exp(-(a + A) x^2 - 2(b + B)xy - (c + C)y^2)>
 * i.e. the inverses of the covariances matrices add.
 *
 * We know sigmaXX_ow and sigmaXXW, the covariances of the weighted object
 * and of the weights themselves.  We can estimate the object's covariance as
 *   sigmaXX_ow^-1 - sigmaXXW^-1
 * and, as we want to find a set of weights with the _same_ covariance as the
 * object we take this to be the an estimate of our correct weights.
 *
 * N.b. This assumes that the object is roughly Gaussian.
 * Consider the object:
 *   O == delta(x + p) + delta(x - p)
 * the covariance of the weighted object is equal to that of the unweighted
 * object, and this prescription fails badly.  If we detect this, we set
 * the Flags::SHAPE_UNWEIGHTED bit, and calculate the UNweighted moments
 * instead.
 */
        {
            double det_n;                       // determinant of nXX matrix
            double det_ow;                      // determinant of sigmaXX_ow matrix
            float n11, n12, n22;                // elements of inverse of next guess at weighting function
            float ow11, ow12, ow22;             // elements of inverse of sigmaXX_ow

            det_ow = sigma11_ow*sigma22_ow - sigma12_ow*sigma12_ow;

            if (det_ow <= 0) {
                shape->setFlags(shape->getFlags() | Flags::SHAPE_UNWEIGHTED);
                break;
            }
         
            ow11 =  sigma22_ow/det_ow;
            ow12 = -sigma12_ow/det_ow;
            ow22 =  sigma11_ow/det_ow;

            n11 = ow11 - w11;
            n12 = ow12 - w12;
            n22 = ow22 - w22;
            det_n = n11*n22 - n12*n12;

            if (det_n <= 0) {            // product-of-Gaussians assumption failed
                shape->setFlags(shape->getFlags() | Flags::SHAPE_UNWEIGHTED);
                break;
            }
      
            sigma11W = n22/det_n;
            sigma12W = -n12/det_n;
            sigma22W = n11/det_n;
        }

        if (sigma11W <= 0 || sigma22W <= 0) {
            shape->setFlags(shape->getFlags() | Flags::SHAPE_UNWEIGHTED);
            break;
        }
    }

    if (iter == MAXIT) {
        shape->setFlags(shape->getFlags() | Flags::SHAPE_MAXITER | Flags::SHAPE_UNWEIGHTED);
    }

    if (sumxx + sumyy == 0.0) {
        shape->setFlags(shape->getFlags() | Flags::SHAPE_UNWEIGHTED);
    }
/*
 * Problems; try calculating the un-weighted moments
 */
    if (shape->getFlags() & Flags::SHAPE_UNWEIGHTED) {
        w11 = w22 = w12 = 0;
        if (calcmom(image, xcen, ycen, bbox, bkgd, interpflag, w11, w12, w22,
                   &sum, &sumx, &sumy, &sumxx, &sumxy, &sumyy, NULL) < 0 || sum <= 0) {
            shape->setFlags((shape->getFlags() & ~Flags::SHAPE_UNWEIGHTED) | Flags::SHAPE_UNWEIGHTED_BAD);

            if (sum > 0) {
                shape->setIxx(1/12.0);      // a single pixel
                shape->setIxy(0.0);
                shape->setIyy(1/12.0);
            }
            
            return false;
        }

        sigma11W = sumxx/sum;          // estimate of object moments
        sigma12W = sumxy/sum;          //   usually, object == weight
        sigma22W = sumyy/sum;          //      at this point
    }

    shape->setI0(ampW);
    shape->setIxx(sigma11W);
    shape->setIxy(sigma12W);
    shape->setIyy(sigma22W);
    shape->setIxy4(sums4/sum);

    if (shape->getIxx() + shape->getIyy() != 0.0) {
        int const ix = lsst::afw::image::positionToIndex(xcen);
        int const iy = lsst::afw::image::positionToIndex(ycen);
        
        if (ix >= 0 && ix < mimage.getWidth() && iy >= 0 && iy < mimage.getHeight()) {
            float const bkgd_var =
                ImageAdaptor<ImageT>().getVariance(mimage, ix, iy); // XXX Overestimate as it includes object

            if (bkgd_var > 0.0) {                                   // NaN is not > 0.0
                if (!(shape->getFlags() & Flags::SHAPE_UNWEIGHTED)) {
                    detail::SdssShapeImpl::Matrix4 fisher = calc_fisher(*shape, bkgd_var); // Fisher matrix 
                    shape->setCovar(fisher.inverse());
                }
            }
        }
    }

    return true;
}
}

namespace {

/**
 * @brief A class that knows how to calculate the SDSS adaptive moment shape measurements
 */
class SdssShape : public afwDetection::Shape
{
public:
    typedef boost::shared_ptr<SdssShape> Ptr;
    typedef boost::shared_ptr<SdssShape const> ConstPtr;

    /// Ctor
    SdssShape(double x, double xErr, double y, double yErr,
              double ixx, double ixxErr, double ixy, double ixyErr, double iyy, double iyyErr) :
        afwDetection::Shape(x, xErr, y, yErr, ixx, ixxErr, ixy, ixyErr, iyy, iyyErr) {}

    /// Add desired fields to the schema
    virtual void defineSchema(afwDetection::Schema::Ptr schema ///< our schema; == _mySchema
                     ) {
        Shape::defineSchema(schema);
    }

    template<typename ExposureT>
    static Shape::Ptr doMeasure(CONST_PTR(ExposureT) exposure,
                                CONST_PTR(afwDetection::Peak) peak,
                                CONST_PTR(afwDetection::Source)
                               );

    static bool doConfigure(lsst::pex::policy::Policy const& policy)
    {
        if (policy.isDouble("background")) {
            _background = policy.getDouble("background");
        } 
        
        return true;
    }
private:
    static double _background;
    SdssShape(void) : afwDetection::Shape() { }
    LSST_SERIALIZE_PARENT(afwDetection::Shape)
};

LSST_REGISTER_SERIALIZER(SdssShape)

double SdssShape::_background = 0.0;    // the frame's background level

/************************************************************************************************************/
/*
 * Decide on the bounding box for the region to examine while calculating
 * the adaptive moments
 */
lsst::afw::geom::BoxI set_amom_bbox(int width, int height, // size of region
                                     float xcen, float ycen,        // centre of object
                                     double sigma11_w,              // quadratic moments of the
                                     double ,                       //         weighting function
                                     double sigma22_w,              //                    xx, xy, and yy
                                     float maxRad = 1000              // Maximum radius of area to use
                                    )
{
    float rad = 4*sqrt(((sigma11_w > sigma22_w) ? sigma11_w : sigma22_w));
        
    if (rad > maxRad) {
        rad = maxRad;
    }
        
    int ix0 = static_cast<int>(xcen - rad - 0.5);
    ix0 = (ix0 < 0) ? 0 : ix0;
    int iy0 = static_cast<int>(ycen - rad - 0.5);
    iy0 = (iy0 < 0) ? 0 : iy0;
    lsst::afw::geom::Point2I llc(ix0, iy0); // Desired lower left corner
        
    int ix1 = static_cast<int>(xcen + rad + 0.5);
    if (ix1 >= width) {
        ix1 = width - 1;
    }
    int iy1 = static_cast<int>(ycen + rad + 0.5);
    if (iy1 >= height) {
        iy1 = height - 1;
    }
    lsst::afw::geom::Point2I urc(ix1, iy1); // Desired upper right corner
        
    return lsst::afw::geom::BoxI(llc, urc);
}   

/*****************************************************************************/
/*
 * Calculate weighted moments of an object up to 2nd order
 */
template<typename ImageT>
static int
calcmom(ImageT const& image,            // the image data
        float xcen, float ycen,         // centre of object
        lsst::afw::geom::BoxI bbox,    // bounding box to consider
        float bkgd,                     // data's background level
        bool interpflag,                // interpolate within pixels?
        double w11, double w12, double w22, // weights
        double *psum, double *psumx, double *psumy, // sum w*I, sum [xy]*w*I
        double *psumxx, double *psumxy, double *psumyy, // sum [xy]^2*w*I
        double *psums4                                  // sum w*I*weight^2 or NULL
       )
{
    
    float tmod, ymod;
    float X, Y;                          // sub-pixel interpolated [xy]
    float weight;
    float tmp;
    double sum, sumx, sumy, sumxx, sumyy, sumxy, sums4;
#define RECALC_W 0                      // estimate sigmaXX_w within BBox?
#if RECALC_W
    double wsum, wsumxx, wsumxy, wsumyy;

    wsum = wsumxx = wsumxy = wsumyy = 0;
#endif

    assert(w11 >= 0);                   // i.e. it was set
    if (fabs(w11) > 1e6 || fabs(w12) > 1e6 || fabs(w22) > 1e6) {
        return(-1);
    }

    sum = sumx = sumy = sumxx = sumxy = sumyy = sums4 = 0;

    int const ix0 = bbox.getMinX();       // corners of the box being analyzed
    int const ix1 = bbox.getMaxX();
    int const iy0 = bbox.getMinY();       // corners of the box being analyzed
    int const iy1 = bbox.getMaxY();

    if (ix0 < 0 || ix1 >= image.getWidth() || iy0 < 0 || iy1 >= image.getHeight()) {
        return -1;
    }

    for (int i = iy0; i <= iy1; ++i) {
        typename ImageT::x_iterator ptr = image.x_at(ix0, i);
        float const y = i - ycen;
        float const y2 = y*y;
        float const yl = y - 0.375;
        float const yh = y + 0.375;
        for (int j = ix0; j <= ix1; ++j, ++ptr) {
            float x = j - xcen;
            if (interpflag) {
                float const xl = x - 0.375;
                float const xh = x + 0.375;
               
                float expon = xl*xl*w11 + yl*yl*w22 + 2.0*xl*yl*w12;
                tmp = xh*xh*w11 + yh*yh*w22 + 2.0*xh*yh*w12;
                expon = (expon > tmp) ? expon : tmp;
                tmp = xl*xl*w11 + yh*yh*w22 + 2.0*xl*yh*w12;
                expon = (expon > tmp) ? expon : tmp;
                tmp = xh*xh*w11 + yl*yl*w22 + 2.0*xh*yl*w12;
                expon = (expon > tmp) ? expon : tmp;
               
                if (expon <= 9.0) {
                    tmod = *ptr - bkgd;
                    for (Y = yl; Y <= yh; Y += 0.25) {
                        double const interpY2 = Y*Y;
                        for (X = xl; X <= xh; X += 0.25) {
                            double const interpX2 = X*X;
                            double const interpXy = X*Y;
                            expon = interpX2*w11 + 2*interpXy*w12 + interpY2*w22;
                            weight = exp(-0.5*expon);
                           
                            ymod = tmod*weight;
                            sum += ymod;
                            sumx += ymod*(X + xcen);
                            sumy += ymod*(Y + ycen);
#if RECALC_W
                            wsum += weight;
                           
                            tmp = interpX2*weight;
                            wsumxx += tmp;
                            sumxx += tmod*tmp;
                           
                            tmp = interpXy*weight;
                            wsumxy += tmp;
                            sumxy += tmod*tmp;
                           
                            tmp = interpY2*weight;
                            wsumyy += tmp;
                            sumyy += tmod*tmp;
#else
                            sumxx += interpX2*ymod;
                            sumxy += interpXy*ymod;
                            sumyy += interpY2*ymod;
#endif
                            sums4 += expon*expon*ymod;
                        }
                    }
                }
            } else {
                float x2 = x*x;
                float xy = x*y;
                float expon = x2*w11 + 2*xy*w12 + y2*w22;
               
                if (expon <= 14.0) {
                    weight = exp(-0.5*expon);
                    tmod = *ptr - bkgd;
                    ymod = tmod*weight;
                    sum += ymod;
                    sumx += ymod*j;
                    sumy += ymod*i;
#if RECALC_W
                    wsum += weight;
                   
                    tmp = x2*weight;
                    wsumxx += tmp;
                    sumxx += tmod*tmp;
                   
                    tmp = xy*weight;
                    wsumxy += tmp;
                    sumxy += tmod*tmp;
                   
                    tmp = y2*weight;
                    wsumyy += tmp;
                    sumyy += tmod*tmp;
#else
                    sumxx += x2*ymod;
                    sumxy += xy*ymod;
                    sumyy += y2*ymod;
#endif
                    sums4 += expon*expon*ymod;
                }
            }
        }
    }
   
    *psum = sum;
    *psumx = sumx;
    *psumy = sumy;
    *psumxx = sumxx;
    *psumxy = sumxy;
    *psumyy = sumyy;
    if (psums4 != NULL) {
        *psums4 = sums4;
    }

#if RECALC_W
    if (wsum > 0) {
        double det = w11*w22 - w12*w12;
        wsumxx /= wsum;
        wsumxy /= wsum;
        wsumyy /= wsum;
        printf("%g %g %g  %g %g %g\n", w22/det, -w12/det, w11/det, wsumxx, wsumxy, wsumyy);
    }
#endif

    return((sum > 0 && sumxx > 0 && sumyy > 0) ? 0 : -1);
}

/************************************************************************************************************/
/**
 * @brief Given an image and a pixel position, return a Shape using the SDSS algorithm
 */
template<typename ExposureT>
afwDetection::Shape::Ptr SdssShape::doMeasure(CONST_PTR(ExposureT) exposure,
                                              CONST_PTR(afwDetection::Peak) peak,
                                              CONST_PTR(afwDetection::Source)
                                             )
{
    if (!peak) {
        double const NaN = std::numeric_limits<double>::quiet_NaN();
        return boost::shared_ptr<SdssShape>(new SdssShape(NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN));
    }

    typedef typename ExposureT::MaskedImageT MaskedImageT;
    MaskedImageT const& mimage = exposure->getMaskedImage();

    double xcen = peak->getFx();         // object's column position
    double ycen = peak->getFy();         // object's row position

    xcen -= mimage.getX0();             // work in image Pixel coordinates
    ycen -= mimage.getY0();
    
    float shiftmax = 1;                 // Max allowed centroid shift \todo XXX set shiftmax from Policy
    if (shiftmax < 2) {
        shiftmax = 2;
    } else if (shiftmax > 10) {
        shiftmax = 10;
    }

    detail::SdssShapeImpl shapeImpl;
    (void)detail::getAdaptiveMoments(mimage, _background, xcen, ycen, shiftmax, &shapeImpl);
/*
 * We need to measure the PSF's moments even if we failed on the object
 * N.b. This isn't yet implemented (but the code's available from SDSS)
 */
    double const x = shapeImpl.getX();
    double const xErr = shapeImpl.getXErr();
    double const y = shapeImpl.getY();
    double const yErr = shapeImpl.getYErr();
    double const ixx = shapeImpl.getIxx();
    double const ixy = shapeImpl.getIxy();
    double const iyy = shapeImpl.getIyy();
    double const ixxErr = shapeImpl.getIxxErr();
    double const ixyErr = shapeImpl.getIxyErr();
    double const iyyErr = shapeImpl.getIyyErr();

    /*
     * Can't use boost::make_shared here as it's limited to 9 arguments
     */
    PTR(SdssShape) shape = boost::shared_ptr<SdssShape>(new SdssShape(x, xErr, y, yErr,
                                                      ixx, ixxErr, ixy, ixyErr, iyy, iyyErr));
    shape->set<SHAPE_STATUS, short>(shapeImpl.getFlags());

    return shape;
}

/*
 * Declare the existence of a "SDSS" algorithm to MeasureShape
 *
 * @cond
 */
#define INSTANTIATE(TYPE) \
    MeasureShape<afwImage::Exposure<TYPE> >::declare("SDSS", \
        &SdssShape::doMeasure<afwImage::Exposure<TYPE> >, \
        &SdssShape::doConfigure \
        )

volatile bool isInstance[] = {
    INSTANTIATE(int),
    INSTANTIATE(float),
    INSTANTIATE(double)
};

// \endcond

}

#undef INSTANTIATE
#define INSTANTIATE(TYPE) \
    template bool detail::getAdaptiveMoments<lsst::afw::image::MaskedImage<TYPE> >(afwImage::MaskedImage<TYPE> const&, double, double, double, double, detail::SdssShapeImpl*); \
    template bool detail::getAdaptiveMoments<lsst::afw::image::Image<TYPE> >(afwImage::Image<TYPE> const&, double, double, double, double, detail::SdssShapeImpl*); \

INSTANTIATE(int);
INSTANTIATE(float);
INSTANTIATE(double);

}}}

