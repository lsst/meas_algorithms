/**
 * Measure adaptive moments.
 *
 * Originally provided by Phil Fischer, based on code from Tim McKay's
 * group.  Error calculations by Dave Johnston.  Major reworking by RHL
 * for SDSS, and now a major rewrite for LSST
 */
#include <limits>
#include "Eigen/LU"
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image.h"
#include "lsst/meas/algorithms/Measure.h"

#include "SdssShape.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;

namespace lsst { namespace meas { namespace algorithms {

/**
 * @brief the (unique) instance of SdssmeasureShape
 */
template<typename ImageT> SdssmeasureShape<ImageT>* SdssmeasureShape<ImageT>::_instance = 0;

/************************************************************************************************************/
/*****************************************************************************/
/*
 * Decide on the bounding box for the region to examine while calculating
 * the adaptive moments
 */
namespace {
    lsst::afw::image::BBox set_amom_bbox(int width, int height, // size of region
                                         float xcen, float ycen,	// centre of object
                                         double sigma11_w,		// quadratic moments of the
                                         double sigma12_w,		//         weighting function
                                         double sigma22_w,		//                    xx, xy, and yy
                                         float maxRad=1000              // Maximum radius of area to use
        ) {
        float rad = 4*sqrt(((sigma11_w > sigma22_w) ? sigma11_w : sigma22_w));
        
        assert(sigma12_w == sigma12_w);	// pretend to use this argument
        
        if(rad > maxRad) {
            rad = maxRad;
        }
        
        int ix0 = xcen - rad - 0.5;
        ix0 = (ix0 < 0) ? 0 : ix0;
        int iy0 = ycen - rad - 0.5;
        iy0 = (iy0 < 0) ? 0 : iy0;
        lsst::afw::image::PointI llc(ix0, iy0); // Desired lower left corner
        
        int ix1 = xcen + rad + 0.5;
        if(ix1 >= width) {
            ix1 = width - 1;
        }
        int iy1 = ycen + rad + 0.5;
        if(iy1 >= height) {
            iy1 = height - 1;
        }
        lsst::afw::image::PointI urc(ix1, iy1); // Desired upper right corner
        
        return lsst::afw::image::BBox(llc, urc);
    }   
}

/*****************************************************************************/
/*
 * Calculate weighted moments of an object up to 2nd order
 */
template<typename ImageT>
static int
calcmom(ImageT const& image,		// the image data
	float xcen, float ycen,		// centre of object
        lsst::afw::image::BBox bbox,    // bounding box to consider
	float bkgd,			// data's background level
	bool interpflag,                // interpolate within pixels?
	double w11, double w12, double w22, /* weights */
	double *psum, double *psumx, double *psumy, /* sum w*I, sum [xy]*w*I */
	double *psumxx, double *psumxy, double *psumyy, /* sum [xy]^2*w*I */
	double *psums4)			/* sum w*I*weight^2 or NULL */
{   
    float tmod,ymod;
    float X,Y;				/* sub-pixel interpolated [xy] */
    float weight;
    float tmp;
    double sum,sumx,sumy,sumxx,sumyy,sumxy,sums4;
#define RECALC_W 0			/* estimate sigmaXX_w within BBox? */
#if RECALC_W
    double wsum, wsumxx, wsumxy, wsumyy;

    wsum = wsumxx = wsumxy = wsumyy = 0;
#endif

    assert(w11 >= 0);			/* i.e. it was set */
    if(fabs(w11) > 1e6 || fabs(w12) > 1e6 || fabs(w22) > 1e6) {
        return(-1);
    }

    sum = sumx = sumy = sumxx = sumxy = sumyy = sums4 = 0;

    int const ix0 = bbox.getX0();       /* corners of the box being analyzed */
    int const ix1 = bbox.getX1();
    int const iy0 = bbox.getY0();       /* corners of the box being analyzed */
    int const iy1 = bbox.getY1();
   
    for (int i = iy0; i <= iy1; ++i) {
        typename ImageT::x_iterator ptr = image.x_at(ix0, i);
        float const y = i - ycen;
        float const y2 = y*y;
        float const yl = y - 0.375;
        float const yh = y + 0.375;
        for (int j = ix0;j <= ix1; ++j, ++ptr) {
            float x = j - xcen;
            if(interpflag) {
                float const xl = x - 0.375;
                float const xh = x + 0.375;
               
                float expon = xl*xl*w11 + yl*yl*w22 + 2.*xl*yl*w12;
                tmp = xh*xh*w11 + yh*yh*w22 + 2.*xh*yh*w12;
                expon = (expon > tmp) ? expon : tmp;
                tmp = xl*xl*w11 + yh*yh*w22 + 2.*xl*yh*w12;
                expon = (expon > tmp) ? expon : tmp;
                tmp = xh*xh*w11 + yl*yl*w22 + 2.*xh*yl*w12;
                expon = (expon > tmp) ? expon : tmp;
               
                if(expon <= 9.0) {
                    tmod = *ptr - bkgd;
                    for(Y = yl; Y <= yh;Y += 0.25) {
                        double const interp_y2 = Y*Y;
                        for(X = xl; X <= xh; X += 0.25) {
                            double const interp_x2 = X*X;
                            double const interp_xy = X*Y;
                            expon = interp_x2*w11 + 2*interp_xy*w12 + interp_y2*w22;
                            weight = exp(-0.5*expon);
                           
                            ymod = tmod*weight;
                            sum += ymod;
                            sumx += ymod*(X+xcen);
                            sumy += ymod*(Y+ycen);
#if RECALC_W
                            wsum += weight;
                           
                            tmp = interp_x2*weight;
                            wsumxx += tmp;
                            sumxx += tmod*tmp;
                           
                            tmp = interp_xy*weight;
                            wsumxy += tmp;
                            sumxy += tmod*tmp;
                           
                            tmp = interp_y2*weight;
                            wsumyy += tmp;
                            sumyy += tmod*tmp;
#else
                            sumxx += interp_x2*ymod;
                            sumxy += interp_xy*ymod;
                            sumyy += interp_y2*ymod;
#endif
                            sums4 += expon*expon*ymod;
                        }
                    }
                }
            } else {
                float x2 = x*x;
                float xy = x*y;
                float expon = x2*w11 + 2*xy*w12 + y2*w22;
               
                if(expon <= 14.0) {
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
   
    *psum = sum; *psumx = sumx; *psumy = sumy;
    *psumxx = sumxx; *psumxy = sumxy; *psumyy = sumyy;
    if(psums4 != NULL) { *psums4 = sums4; }

#if RECALC_W
    if(wsum > 0) {
        double det = w11*w22 - w12*w12;
        wsumxx /= wsum; wsumxy /= wsum; wsumyy /= wsum;
        printf("%g %g %g  %g %g %g\n", w22/det, -w12/det, w11/det,
               wsumxx, wsumxy, wsumyy);
    }
#endif

    return((sum > 0 && sumxx > 0 && sumyy > 0) ? 0 : -1);
}

#define MAXIT 100                       // \todo from Policy XXX
#if 0
#define TOL1 0.001                      // \todo from Policy XXX
#define TOL2 0.01                       // \todo from Policy XXX
#else  // testing
#define TOL1 0.00001
#define TOL2 0.0001
#endif

/*
 * Workhorse for adaptive moments
 */
template<typename ImageT>
static bool
get_moments(ImageT const& image,        // the data to process
	    float bkgd,			/* background level */
	    float xcen, float ycen,	/* centre of object */
	    float shiftmax,		/* max allowed centroid shift */
            Shape *shape                // the Shape to fill out
           ) {
    float amp_w = 0;			/* amplitude of best-fit Gaussian */
    double sum;				/* sum of intensity*weight */
    double sumx, sumy;			/* sum ((int)[xy])*intensity*weight */
    double sumxx, sumxy, sumyy;		/* sum {x^2,xy,y^2}*intensity*weight */
    double sums4;			/* sum intensity*weight*exponent^2 */
    const float xcen0 = xcen;		/* initial centre */
    const float ycen0 = ycen;		/*                of object */

    double sigma11_w = 1.5;             // quadratic moments of the
    double sigma12_w = 0.0;             //     weighting fcn;
    double sigma22_w = 1.5;             //               xx, xy, and yy

    double w11 = -1, w12 = -1, w22 = -1;        // current weights for moments; always set when iter == 0
    float e1_old = 1e6, e2_old = 1e6;		// old values of shape parameters e1 and e2
    float sigma11_ow_old = 1e6;                 /* previous version of sigma11_ow */
   
    bool interpflag = false;            //* interpolate finer than a pixel?
    lsst::afw::image::BBox bbox;
    int iter = 0;                       // iteration number
    for(; iter < MAXIT; iter++) {
        bbox = set_amom_bbox(image.getWidth(), image.getHeight(), xcen, ycen, sigma11_w, sigma12_w, sigma22_w);

        double const det_w = sigma11_w*sigma22_w - sigma12_w*sigma12_w; // determinant of sigmaXX_w matrix
#if 0					/* this form was numerically unstable on my G4 powerbook */
        assert(det_w >= 0.0);
#else
        assert(sigma11_w*sigma22_w >= sigma12_w*sigma12_w - std::numeric_limits<float>::epsilon());
#endif
        if(det_w < std::numeric_limits<float>::epsilon()) { // a suitably small number
            shape->setFlags(shape->getFlags() | Flags::SHAPE_UNWEIGHTED);
            break;
        }

        {
            const double ow11 = w11;	/* old */
            const double ow12 = w12;	/*     values */
            const double ow22 = w22;	/*            of w11, w12, w22 */

            w11 = sigma22_w/det_w; w12 = -sigma12_w/det_w; w22 = sigma11_w/det_w;

            float const xinterp = 0.25; // I.e. 0.5*0.5
            if(sigma11_w < xinterp || sigma22_w < xinterp || det_w < xinterp*xinterp) {
                if(!interpflag) {
                    interpflag = true;       /* N.b.: stays set for this object */
                    if(iter > 0) {
                        sigma11_ow_old = 1.e6; /* force at least one more iteration*/
                        w11 = ow11; w12 = ow12; w22 = ow22;
                        iter--;		/* we didn't update wXX */
                    }
                }
            }
        }

        if (calcmom(image, xcen, ycen, bbox, bkgd, interpflag, w11, w12, w22,
                    &sum, &sumx, &sumy, &sumxx, &sumxy, &sumyy, &sums4) < 0) {
            shape->setFlags(shape->getFlags() | Flags::SHAPE_UNWEIGHTED);
            break;
        }

        amp_w = sum/(M_PI*sqrt(det_w));
#if 0
/*
 * Find new centre
 *
 * This is only needed if we update the centre; if we use the input position we've already done the work
 */
        xcen = sumx/sum;
        ycen = sumy/sum;
#endif
        shape->getCentroid().setX(sumx/sum); // update centroid.  N.b. we're not setting errors here
        shape->getCentroid().setY(sumy/sum);

        if(fabs(shape->getCentroid().getX() - xcen0) > shiftmax ||
           fabs(shape->getCentroid().getY() - ycen0) > shiftmax) {
            shape->setFlags(shape->getFlags() | Flags::SHAPE_SHIFT);
        }
/*
 * OK, we have the centre. Proceed to find the second moments.
 */
        float const sigma11_ow = sumxx/sum; // quadratic moments of
        float const sigma22_ow = sumyy/sum; //          weight*object
        float const sigma12_ow = sumxy/sum; //                 xx, xy, and yy 

        if(sigma11_ow <= 0 || sigma22_ow <= 0) {
            shape->setFlags(shape->getFlags() | Flags::SHAPE_UNWEIGHTED);
            break;
        }

        float const d = sigma11_ow + sigma22_ow; // current values of shape parameters
        float const e1 = (sigma11_ow - sigma22_ow)/d;
        float const e2 = 2.*sigma12_ow/d;
/*
 * Did we converge?
 */
        if(iter > 0 &&
           fabs(e1 - e1_old) < TOL1 && fabs(e2 - e2_old) < TOL1 &&
           fabs(sigma11_ow/sigma11_ow_old - 1.) < TOL2 ) {
            break;				/* yes; we converged */
        }

        e1_old = e1; e2_old = e2; sigma11_ow_old = sigma11_ow;
/*
 * Didn't converge, calculate new values for weighting function
 *
 * The product of two Gaussians is a Gaussian:
 * <x^2 exp(-a x^2 - 2bxy - cy^2) exp(-Ax^2 - 2Bxy - Cy^2)> = 
 *                            <x^2 exp(-(a + A) x^2 - 2(b + B)xy - (c + C)y^2)>
 * i.e. the inverses of the covariances matrices add.
 *
 * We know sigmaXX_ow and sigmaXX_w, the covariances of the weighted object
 * and of the weights themselves.  We can estimate the object's covariance as
 *   sigmaXX_ow^-1 - sigmaXX_w^-1
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
            double det_n;			/* determinant of nXX matrix */
            double det_ow;			/* determinant of sigmaXX_ow matrix */
            float n11, n12, n22;		/* elements of inverse of
                                                   next guess at weighting function */
            float ow11, ow12, ow22;	/* elements of inverse of sigmaXX_ow */

            det_ow = sigma11_ow*sigma22_ow - sigma12_ow*sigma12_ow;

            if(det_ow <= 0) {
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

            if(det_n <= 0) {		/* product-of-Gaussians
                                           assumption failed */
                shape->setFlags(shape->getFlags() | Flags::SHAPE_UNWEIGHTED);
                break;
            }
      
            sigma11_w = n22/det_n;
            sigma12_w = -n12/det_n;
            sigma22_w = n11/det_n;
        }

        if(sigma11_w <= 0 || sigma22_w <= 0) {
            shape->setFlags(shape->getFlags() | Flags::SHAPE_UNWEIGHTED);
            break;
        }
    }

    if(iter == MAXIT) {
        shape->setFlags(shape->getFlags() | Flags::SHAPE_MAXITER | Flags::SHAPE_UNWEIGHTED);
    }

    if(sumxx + sumyy == 0.0) {
        shape->setFlags(shape->getFlags() | Flags::SHAPE_UNWEIGHTED);
    }
/*
 * Problems; try calculating the un-weighted moments
 */
    if(shape->getFlags() & Flags::SHAPE_UNWEIGHTED) {
        w11 = w22 = w12 = 0;
        if(calcmom(image, xcen, ycen, bbox, bkgd, interpflag, w11, w12, w22,
                   &sum, &sumx, &sumy, &sumxx, &sumxy, &sumyy, NULL) < 0 || sum <= 0) {
            shape->setFlags((shape->getFlags() & ~Flags::SHAPE_UNWEIGHTED) | Flags::SHAPE_UNWEIGHTED_BAD);
            return false;
        }

        sigma11_w = sumxx/sum;		/* estimate of object moments */
        sigma12_w = sumxy/sum;		/*   usually, object == weight */
        sigma22_w = sumyy/sum;		/*      at this point */
    }

    shape->setM0(amp_w);
    shape->setMxx(sigma11_w);
    shape->setMxy(sigma12_w);
    shape->setMyy(sigma22_w);
    shape->setMxy4(sums4/sum);

    return true;
}

/*****************************************************************************/
/*
 * Error analysis, courtesy of David Johnston, University of Chicago
 */
/*
 * This function takes the 4 Gaussian parameters A, sigmaXX_w and the
 * sky variance and fills in the Fisher matrix from the least squares fit.
 *
 * Following "Numerical Recipes in C" section 15.5, it ignores the 2nd
 * derivative parts and so the fisher matrix is just a function of these
 * best fit model parameters. The components are calculated analytically.
 */
Shape::Matrix4 calc_fisher(Shape *shape,        // the Shape that we want the the Fisher matrix for
                         float bkgd_var         // background variance level for object
                        ) {
    float const A = shape->getM0();           /* amplitude */
    float const sigma11_w = shape->getMxx();
    float const sigma12_w = shape->getMxy();
    float const sigma22_w = shape->getMyy();
    
    double const D = sigma11_w*sigma22_w - sigma12_w*sigma12_w;
   
    if(D <= std::numeric_limits<double>::epsilon()) {
        throw LSST_EXCEPT(lsst::pex::exceptions::DomainErrorException,
                          "Determinant is too small calculating Fisher matrix");
    }
/*
 * a normalization factor
 */
    assert(bkgd_var > 0.0);
    double const F = M_PI*sqrt(D)/bkgd_var;
/*
 * Calculate the 10 independent elements of the 4x4 Fisher matrix 
 */
    Shape::Matrix4 fisher;

    double fac = F*A/(4.0*D);
    fisher(0, 0) =  F;
    fisher(0, 1) =  fac*sigma22_w;
    fisher(1, 0) =  fisher(0, 1);
    fisher(0, 2) =  fac*sigma11_w;			
    fisher(2, 0) =  fisher(0, 2);
    fisher(0, 3) = -fac*2*sigma12_w;	
    fisher(3, 0) =  fisher(0, 3);
    
    fac = 3.0*F*A*A/(16.0*D*D);
    fisher(1, 1) =  fac*sigma22_w*sigma22_w;
    fisher(2, 2) =  fac*sigma11_w*sigma11_w;
    fisher(3, 3) =  fac*4.0*(sigma12_w*sigma12_w + D/3.0);
    
    fisher(1, 2) =  fisher(3, 3)/4.0;
    fisher(2, 1) =  fisher(1, 2);
    fisher(1, 3) =  fac*(-2*sigma22_w*sigma12_w);
    fisher(3, 1) =  fisher(1, 3);
    fisher(2, 3) =  fac*(-2*sigma11_w*sigma12_w);
    fisher(3, 2) =  fisher(2, 3);
    
    return fisher;
}

/************************************************************************************************************/
/**
 * @brief Given an image and a pixel position, return a Shape using the SDSS algorithm
 */
template<typename ImageT>
Shape SdssmeasureShape<ImageT>::doApply(ImageT const& image, ///< The Image wherein dwells the object
                                       double xcen,          ///< object's column position
                                       double ycen,          ///< object's row position
                                       PSF const*,           ///< image's PSF
                                       double background     ///< image's background level
                                      ) const {
    xcen -= image.getX0();              // work in image Pixel coordinates
    ycen -= image.getY0();
    
    float shiftmax = 1;;                /// Max allowed centroid shift \todo XXX set shiftmax from Policy
    if(shiftmax < 2) {
        shiftmax = 2;
    } else if(shiftmax > 10) {
        shiftmax = 10;
    }

#if 1
    float const bkgd_var = 0;            /// \todo XXX set background variance
#else
    float const bkgd_var = sky/gain + dark_variance; // background per-pixel variance
#endif

    Shape shape;                         // The shape to return
    bool status = get_moments(image, background, xcen, ycen, shiftmax, &shape);
/*
 * We need to measure the PSF's moments even if we failed on the object
 * N.b. This isn't yet implemented (but the code's available from SDSS)
 */
    if(!status) {
        return shape;                    // failed
    }

    if(shape.getMxx() + shape.getMyy() != 0.0) {
        if(!shape.getFlags() & Flags::SHAPE_UNWEIGHTED) {
            Shape::Matrix4 fisher = calc_fisher(&shape, bkgd_var); // Fisher matrix 
            shape.setCovar(fisher.inverse());
        }
    }

    return shape;
}

//
// Explicit instantiations
//
// We need to make an instance here so as to register it with measureShape
//
// \cond
#define MAKE_SHAPEFINDERS(IMAGE_T) \
                namespace { \
                    measureShape<lsst::afw::image::Image<IMAGE_T> >* foo = \
                        SdssmeasureShape<lsst::afw::image::Image<IMAGE_T> >::getInstance(); \
                }
                
MAKE_SHAPEFINDERS(float)


// \endcond

}}}
