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
 
#include <algorithm>
#include "Eigen/Core"
#include "Eigen/LU"

#include "all.h"

#define USE_WEIGHT 0                    // zweight is only set, not used.  It isn't even set if this is false
struct Raster {
    int x, y;                           // x, y index of rastered pixel
    double zo;                          // pixel's value
#if USE_WEIGHT
    double zweight;                     // weight for the pixel
#endif
};

struct Fit2d {
    typedef Eigen::Matrix<double, FittedModel::NPARAM, FittedModel::NPARAM> Matrix;
    typedef Eigen::Matrix<double, FittedModel::NPARAM, 1> Vector;

    template<typename PixelT>
    explicit Fit2d(afwImage::Image<PixelT> const& im) : wide(32), rast(*new std::vector<Raster>(wide*wide)) {
        ncols = im.getWidth();
        nrows = im.getHeight();

        dc2zmin = 0.3;
        dc2zmax = 3.5;
        tooclose = 3;
        toosmall = 1.5;
        nitmax = 15;
        flamdok = 1.0e-7;
        ratiomin = -1.0e-7;
        lost = 3.5;
        xlamd = 5.0;

        iter = 0;
        flamd = 1.0;
    }
    ~Fit2d() {
        delete &rast;
    }

    int wide;                           // patch of image being fitted is wide*wide
    std::vector<Raster> &rast;          // The pixels being fitted, thought of as a raster scan
    int nin;
    int nref;
    int iter;
    int tooclose;
    Vector param;
    Matrix alpha;
    Vector beta;
    Vector elold;
    Vector elnew;
    Vector sgold;
    Vector sgnew;
    Matrix alold;
    Matrix alnew;
    Vector beold;
    Vector benew;
    double chisq;
    double chold;
    double stold;
    double chnew;
    double stnew;
    double flamd;
    double flamdok;
    double xlamd;
    double ratiomin;
    int nitmax;
    double dc2zmin;
    double dc2zmax;
    double lost;
    double toosmall;                    // minimum acceptable width
    int ncols;                          // number of columns in image
    int nrows;                          // number of rows in image
    int status;                         // status of processing
};

/************************************************************************************************************/
/**
 * Set alpha and beta
 * From DCALC2.F
 */
static void dcalc2(Fit2d *fit,                  // the struct carrying all the fit information
                   Fit2d::Vector const &el,     // current parameters of the fit
                   Fit2d::Matrix *alpha,        // desired alpha matrix
                   Fit2d::Vector *beta          // desired beta vector
                  ) {
/*
 * Pass arguments and initialize
 */
    double const a = el[FittedModel::PEAK];
    double const b = el[FittedModel::SKY];
    double const x0 = el[FittedModel::X0];
    double const y0 = el[FittedModel::Y0];
    double const s = el[FittedModel::SIGMA];
    
    if (a <= 0.0) {
        fit->status = FittedModel::BAD_A;
        return;
    } else if (s <= 0.0) {
        fit->status = FittedModel::BAD_WIDTH;
        return;
    } else {
        fit->status = 0;
    }
    
    beta->setZero();
    alpha->setZero();

    fit->nin = 0;
    fit->chnew = 0.0;
/*
 * Examine all pixels in raster
 */
    for (int i = 0, nrast = fit->rast.size(); i != nrast; i++) {
        double const dx = (fit->rast[i].x - x0)/s;
        double const dy = (fit->rast[i].y - y0)/s;
        double const d = hypot(dx, dy);
        
        if (d >= fit->dc2zmin && d <= fit->dc2zmax) {
            double const arg = exp(-d*d/2.0);
            double const funct = a*arg + b;

            Fit2d::Vector df;
            df << arg, 1.0, a*arg*dx/s, a*arg*dy/s, a*arg*(dx*dx/s + dy*dy/s);

            double const r = fit->rast[i].zo - funct; // residual
            fit->chnew += r*r;
            *beta += r*df;
            *alpha += df*df.transpose();    // outer product

            fit->nin++;
        }
    }
    int nu = fit->nin - (FittedModel::NPARAM + 1); // number of degrees of freedom
    if (nu <= 0) {
        fit->status = FittedModel::TOO_FEW;
        return;
    }
    fit->stnew = sqrt(fit->chnew/nu);
}

/************************************************************************************************************/
/*
 * From Bevington's CURFIT
 */
static void curf2(Fit2d *fit) {
/*
 * Initialization
 */
    fit->status = 0;
    if (fit->iter == 0) {
        dcalc2(fit, fit->elnew, &fit->alnew, &fit->benew);
        if (fit->status != 0) {
            return;
        }
        fit->nref = fit->nin;
        fit->stnew = sqrt(fit->chnew/(fit->nref - (FittedModel::NPARAM + 1)));
    }
/*
 * Save Current Parameters
 */
    fit->chold = fit->chnew;
    fit->stold = fit->stnew;
    fit->elold = fit->elnew;
    fit->sgold = fit->sgnew;
    fit->beold = fit->benew;
    fit->alold = fit->alnew;
/*
 * Previous Call To DCALC Used To Fill Current
 */
    int const NPLIM = 4;
    for (int poor = 0; poor != NPLIM; ++poor) {
        for (int j = 0; j != FittedModel::NPARAM; ++j) {
            if (fit->alnew(j, j) == 0.0) {
                fit->status = FittedModel::DIAGONAL;
                return;
            }
            for (int k = 0; k != FittedModel::NPARAM; ++k) {
                fit->alpha(j, k) = fit->alnew(j, k)/sqrt(fit->alnew(j, j)*fit->alnew(k, k));
            }
            fit->alpha(j, j) = 1.0 + fit->flamd;
        }
/*
 * Inversion.
 */
#if 0
        fit->alpha.computeInverse();    // has no way to check if inverse succeeded
#else
        Eigen::FullPivLU<Fit2d::Matrix> alphaLU(fit->alpha);
        if (!alphaLU.isInvertible()) {
            fit->status = -1;
            return;
        }
        fit->alpha = alphaLU.inverse();
#endif

/*
 * New Elements
 */
        for (int j = 0; j != FittedModel::NPARAM; ++j) {
            for (int k = 0; k != FittedModel::NPARAM; ++k) {
                fit->elnew[j] += fit->benew[k]*fit->alpha(j,k)/sqrt(fit->alnew(j, j)*fit->alnew(k, k));
            }
        }
/*
 * Compute Current Chi-Squared And Next ALPHA+BETA
 */
        dcalc2(fit, fit->elnew, &fit->alnew, &fit->benew);
        if (fit->status == FittedModel::TOO_FEW) {
            return;
        }
        fit->stnew = sqrt(fit->chnew/(fit->nref - (FittedModel::NPARAM + 1)));

        for (int j = 0; j != FittedModel::NPARAM; ++j) {
            fit->sgnew[j] = fit->stnew*sqrt(fabs(fit->alpha(j, j)/fit->alnew(j, j)));
        }
/*
 * Quick Return If Solution Got Better
 */
        if (fit->status == 0.0 && fit->stnew <= fit->stold) {
            fit->flamd = fit->flamd/fit->xlamd;
            return;
        }
/*
 * Undo Solution And Try Again
 */
        fit->flamd = 3.0*fit->flamd;
        fit->chnew = fit->chold;
        fit->stnew = fit->stold;
        fit->elnew = fit->elold;
        fit->sgnew = fit->sgold;
        fit->benew = fit->beold;
        fit->alnew = fit->alold;
    }
}

/************************************************************************************************************/
/*
 * Test of convergence or fatal errors
 */
static void cond2(Fit2d *fit) {
/*
 * Ignore CURF errors for these conditions
 */
    if (fit->iter <= 3) {
        fit->status = 0;
        return;
    }
    if (fit->flamd < fit->flamdok) {
        fit->status = FittedModel::ALMOST;
        return;
    }
/*
 * Catch Fatal Errors
 */
    if (fit->status < 0) {
        return;
    }
    if (fit->nin < FittedModel::NPARAM + 1) {
        fit->status = FittedModel::TOO_FEW;
        return;
    }
    if (fit->chnew <= 0.0) {
        fit->status = FittedModel::CHI_SQUARED;
        return;
    }
    if (fit->elnew[FittedModel::X0] < 0.0 || fit->elnew[FittedModel::X0] > fit->ncols ||
        fit->elnew[FittedModel::Y0] < 0.0 || fit->elnew[FittedModel::Y0] > fit->nrows) {
        fit->status = FittedModel::RANGE;
        return;
    }
    if (fit->elnew[FittedModel::SIGMA] < 0.0) {
        fit->status = FittedModel::BAD_WIDTH;
        return;
    }
    
    double const ex = fabs(fit->param[FittedModel::X0] - fit->elnew[FittedModel::X0]);
    double const ey = fabs(fit->param[FittedModel::Y0] - fit->elnew[FittedModel::Y0]);
    if (ex > fit->lost || ey > fit->lost) {
        fit->status = FittedModel::LOST;
        return;
    }
/*
 * Test for convergence
 */
    double const ratio = (fit->stnew - fit->stold)/fit->stnew;
    if (ratio > fit->ratiomin) {
        fit->status = (fit->flamd < 0.001) ? FittedModel::CONVERGE : FittedModel::POOR;
    } else if (fit->iter > fit->nitmax) {
        fit->status = FittedModel::ITERATE;
    } else {
        fit->status = 0;
    }
}

/************************************************************************************************************/
/*
 *  First guess for 2-D Gaussian
 */
template<typename PixelT>
static void fg2(afwImage::Image<PixelT> const& im,     ///< The image
                double x0, double y0,   ///< Initial guess for the position
                Fit2d *fit              ///< Information needed for the fit
               ) {
/*
 * Sanity Check
 */
    int ix0 = static_cast<int>(x0 + 0.5);
    int iy0 = static_cast<int>(y0 + 0.5);
    if(ix0 < fit->tooclose || im.getWidth() - ix0 < fit->tooclose ||
       iy0 < fit->tooclose || im.getHeight() - iy0 < fit->tooclose) {
        fit->status = FittedModel::BAD_GUESS;
        return;
    }
    int jx0 = ix0;
    int jy0 = iy0;
    double peak = im(jx0, jy0);
/*
 * Extract interesting portion
 */
    double const w = fit->wide/2;
    int xl = static_cast<int>(ix0 - w + 0.5);
    if (xl < 0) {
        xl = 0;
    }
    int xh = static_cast<int>(xl + 2*w + 0.5);
    if (xh > im.getWidth()) {
        xh = im.getWidth();
    }
    int yl = static_cast<int>(iy0 - w + 0.5);
    if (yl < 0) {
        yl = 0;
    }
    int yh = static_cast<int>(yl + 2*w + 0.5);
    if (yh > im.getHeight()) {
        yh = im.getHeight();
    }

    double sky = im(xl, yl);
    int put = 0;
    for (int y = yl; y != yh; ++y) {
        for (int x = xl; x != xh; ++x) {
            fit->rast[put].x = x;
            fit->rast[put].y = y;
            fit->rast[put].zo = im(x, y);
#if USE_WEIGHT
            fit->rast[put].zweight = 1.0;
#endif
            if (im(x, y) < sky) {
                sky = im(x, y);
            }
            double const ex = x - ix0;
            double const ey = y - iy0;
            double const er = hypot(ex, ey);
            if (er <= fit->lost) {
                if (im(x, y) > peak) {
                    jx0 = x;
                    jy0 = y;
                    peak = im(x, y);
                }
            }
            put++;
        }
    }
    fit->rast.resize(put);
    ix0 = jx0;
    iy0 = jy0;
/*
 * Look For FWHM
 */
    double xmin = xl;
    double xmax = xh - 1;
    double ymin = yl;
    double ymax = yh - 1;
    double const test = 0.5*(im(ix0, iy0) + sky); // pixel value of half-maximum above sky
    for (int x = ix0; x < xh - 1; ++x) {
        if (im(x + 1, iy0) <= test) {
            double a = test - im(x, iy0);
            double b = im(x + 1, iy0) - im(x, iy0);

            xmax = x + ((b == 0.0) ? 0.5 : a/b);
            break;
        }
    }
    for (int x = ix0; x > 0; --x) {
        if (im(x - 1, iy0) <= test) {
            double a = test - im(x, iy0);
            double b = im(x - 1, iy0) - im(x, iy0);

            xmin = x - ((b == 0.0) ? 0.5 : a/b);
            break;
        }
    }
    for (int y = iy0; y < yh - 1; ++y) {
        if (im(ix0, y + 1) <= test) {
            double a = test - im(ix0, y);
            double b = im(ix0, y + 1) - im(ix0, y);

            ymax = y + ((b == 0.0) ? 0.5 : a/b);
            break;
        }
    }
    for (int y = iy0; y > 0; --y) {
        if (im(ix0, y - 1) <= test) {
            double a = test - im(ix0, y);
            double b = im(ix0, y - 1) - im(ix0, y);

            ymin = y - ((b == 0.0) ? 0.5 : a/b);
            break;
        }
    }
/*
 * Assemble the fitting vector
 */
    fit->param[FittedModel::PEAK] = im(ix0, iy0) - sky;
    fit->param[FittedModel::SKY] = sky;
    fit->param[FittedModel::X0] = x0;
    fit->param[FittedModel::Y0] = y0;
    fit->param[FittedModel::SIGMA] = 0.5*((xmax - xmin)+(ymax - ymin))/2.354;
    if (fit->param[FittedModel::SIGMA] < fit->toosmall) {
        fit->param[FittedModel::SIGMA] = fit->toosmall;
    }
    fit->elnew = fit->param;

    fit->status = 0;
}

/************************************************************************************************************/
/**
 *  Compute centroids with 2-D Gaussian fitter
 */
template<typename PixelT>
FittedModel twodg(afwImage::Image<PixelT> const& im,           ///< The image to fit
                  double x0,                    ///< Initial guess for position, column
                  double y0                     ///<                             and row
                 ) {
/*
 * Initialize the fitting structure
 */
    Fit2d fit(im);
/*
 * Initialization
 */
    fg2(im, x0, y0, &fit);
    if (fit.status != 0) {
        std::vector<double> params(FittedModel::NPARAM);
        std::copy(&fit.elnew[0], &fit.elnew[0] + fit.elnew.size(), params.begin());

        return FittedModel(fit.status, params);
    }
/*
 * Find best-fit model
 */
    for (fit.iter = 0; fit.status == 0; ++fit.iter) {
        curf2(&fit);
        cond2(&fit);
    }

#if 0
    twodgMinuit(&fit);
#endif

    std::vector<double> params(FittedModel::NPARAM);
    std::copy(&fit.elnew[0], &fit.elnew[0] + fit.elnew.size(), params.begin());

    return FittedModel(fit.status, params, fit.iter, fit.flamd, fit.chnew);
}
//
// Explicit instantiations
//
#define MAKE_TWODG(IMAGE_T) \
    template FittedModel twodg(IMAGE_T const& im, double x0, double y0)

MAKE_TWODG(afwImage::Image<float>);
MAKE_TWODG(afwImage::Image<double>);
MAKE_TWODG(afwImage::Image<int>);
