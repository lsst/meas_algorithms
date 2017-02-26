
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

/*!
 * @brief Implementation of PSF code
 *
 * @file
 *
 * @ingroup algorithms
 */
#include <cmath>
#include "lsst/base.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/geom/Angle.h"
#include "lsst/afw/geom/SpanSet.h"
#include "lsst/afw/image/ImagePca.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/math/SpatialCell.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/table/Source.h"
#include "lsst/meas/base/GaussianCentroid.h"
#include "lsst/meas/algorithms/PSF.h"

/************************************************************************************************************/

namespace afwDetection = lsst::afw::detection;
namespace afwGeom = lsst::afw::geom;
namespace afwImage = lsst::afw::image;
namespace afwMath = lsst::afw::math;

namespace lsst {
namespace meas {
namespace algorithms {
/**
 * @brief Constructor for PsfAttributes
 */
PsfAttributes::PsfAttributes(
        CONST_PTR(lsst::afw::detection::Psf) psf, ///< The psf whose attributes we want
        int const iX,                       ///< the x position in the frame we want the attributes at
        int const iY                        ///< the y position in the frame we want the attributes at
                            )
{
    // N.b. (iX, iY) are ints so that we know this image is centered in the central pixel of _psfImage
    _psfImage = psf->computeImage(afwGeom::PointD(iX, iY));
}

/**
 * @brief Constructor for PsfAttributes
 */
PsfAttributes::PsfAttributes(
        CONST_PTR(lsst::afw::detection::Psf) psf, ///< The psf whose attributes we want
        lsst::afw::geom::Point2I const& cen       ///< the position in the frame we want the attributes at
                            ) :
    // N.b. cen is a PointI so that we know this image is centered in the central pixel of _psfImage
    _psfImage(psf->computeImage(afwGeom::PointD(cen)))
{
}

namespace {

/*
 * Return an estimate of <r> == <sqrt(x^2 + y^2)> for an image (i.e. sum(I*r)/sum(I))
 *
 * For a Gaussian N(0, alpha^2),  <r> = sqrt(pi/2) alpha
 */
template<typename ImageT>
double
computeFirstMoment(ImageT const& image,        // the data to process
                   float const xCen, float const yCen // centre of object
                  )
{
    double sum = 0.0;
    double norm = 0.0;
    for (int iY = 0; iY != image->getHeight(); ++iY) {
        int iX = 0;
        for (afwImage::Image<double>::x_iterator ptr = image->row_begin(iY),
                                                 end = image->row_end(iY); ptr != end; ++ptr, ++iX) {
            double const x = iX - xCen;
            double const y = iY - yCen;
            double const r = std::sqrt( x*x + y*y );
            double const m = (*ptr)*r;
            norm += *ptr;
            sum += m;
        }
    }

    std::string errmsg("");
    if (sum < 0.0) {
        errmsg = "sum(I*r) is negative.  ";
    }
    if (norm <= 0.0) {
        errmsg += "sum(I) is <= 0.";
    }
    if (errmsg != "") {
        throw LSST_EXCEPT(lsst::pex::exceptions::DomainError, errmsg);
    }

    return sum/norm;
}

/*
 * Return an estimate of <r^2> == <x^2 + y^2> for an image (i.e. sum(I*r^2)/sum(I))
 *
 * For a Gaussian N(0, alpha^2),  <r^2> = 2 alpha^2
 */
template<typename ImageT>
double
computeSecondMoment(ImageT const& image,        // the data to process
                    float const xCen, float const yCen // centre of object
                   )
{
    double sum = 0.0;
    double norm = 0.0;
    for (int iY = 0; iY != image->getHeight(); ++iY) {
        int iX = 0;
        for (afwImage::Image<double>::x_iterator ptr = image->row_begin(iY),
                                                 end = image->row_end(iY); ptr != end; ++ptr, ++iX) {
            double const x = iX - xCen;
            double const y = iY - yCen;
            double const r2 = x*x + y*y;
            double const m = (*ptr)*r2;
             norm += *ptr;
            sum += m;
        }
    }

    std::string errmsg("");
    if (sum < 0.0) {
        errmsg = "sum(I*r*r) is negative.  ";
    }
    if (norm <= 0.0) {
        errmsg += "sum(I) is <= 0.";
    }
    if (errmsg != "") {
        throw LSST_EXCEPT(lsst::pex::exceptions::DomainError, errmsg);
    }

    return sum/norm;
}

/*****************************************************************************/
/*
 * Calculate weighted moments of an object up to 2nd order
 */
template<typename ImageT>
std::pair<bool, double>
calcmom(ImageT const& image,                // the image data
        float const xCen, float const yCen, // centre of object
        double const w11                    // weights
       )
{
    assert(w11 >= 0);                   /* i.e. it was set */
    if (fabs(w11) > 1e6) {
        return std::make_pair(false, std::numeric_limits<double>::quiet_NaN());
    }

    double sum = 0, sumrr = 0.0;

    for (int i = 0; i < image.getHeight(); ++i) {
        float const y = i - yCen;
        float const y2 = y*y;

        typename ImageT::x_iterator ptr = image.row_begin(i);
        for (int j = 0; j < image.getWidth(); ++j, ++ptr) {
            float const x = j - xCen;
            float const x2 = x*x;
            float const expon = (x2 + y2)*w11;

            if (expon <= 14.0) {
                float const weight = exp(-0.5*expon);
                float const tmod = *ptr;
                float const ymod = tmod*weight;
                sum += ymod;
                sumrr += (x2 + y2)*ymod;
            }
        }
    }

    if (sum <= 0 || sumrr < 0) {
        return std::make_pair(false, std::numeric_limits<double>::quiet_NaN());
    }

    return std::make_pair(true, 0.5*sumrr/sum); // 0.5:  1-D moment
}

/*
 * Return an estimate of <r^2> == <x^2 + y^2> for an image using adaptive moments
 *
 * For a Gaussian N(0, alpha^2),  <r^2> = 2 alpha^2
 *
 * This is basically the SdssShape code simplified for a circularly symmetrical case.  I don't want to call
 * the shape code here as this class may well be moving to afw along with Psf
 */
template<typename ImageT>
double
computeSecondMomentAdaptive(ImageT const& image,        // the data to process
                            float const xCen, float const yCen // centre of object
                      )
{
    int const MAXIT = 100;              // \todo from Policy XXX
    float const TOL = 0.0001;
    double w11 = 0.5;                   // current weight for moments
    float sigma11_ow_old = 1e6;         // previous version of sigma11_ow

    bool unweighted = false;            // do we need to use an unweighted moment?
    int iter = 0;                       // iteration number
    for (; iter < MAXIT; ++iter) {
        std::pair<bool, double> moments = calcmom(*image, xCen, yCen, w11);

        if (not moments.first) {
            unweighted = true;
            break;
        }
/*
 * Did we converge?
 */
        float const sigma11_ow = moments.second; // quadratic moments of weight*object

        if (iter > 0 && fabs(sigma11_ow/sigma11_ow_old - 1.0) < TOL) {
            break;                              /* yes; we converged */
        }

        sigma11_ow_old = sigma11_ow;
/*
 * Didn't converge, calculate new values for weighting function
 *
 * The product of two Gaussians is a Gaussian, the inverse-variances add
 *
 * We know sigma11_ow and sigma11W == 1/w11, the variances of the weighted object
 * and of the weights themselves.  We can estimate the object's variance as
 *   1/sigma11_ow - 1/sigma11W
 * and, as we want to find a set of weights with the _same_ covariance as the
 * object we take this to be the an estimate of our correct weights.
 *
 * N.b. This assumes that the object is roughly Gaussian.
 * Consider the object:
 *   O == delta(x + p) + delta(x - p)
 * the covariance of the weighted object is equal to that of the unweighted
 * object, and this prescription fails badly.  If we detect this, we set
 * unweighted, and calculate the UNweighted moments
 * instead.
 */
        w11 = 1/sigma11_ow - w11;       // inverse of new sigma11_ow

        if (w11 <= 0) {                 // product-of-Gaussians assumption failed
            unweighted = true;
            break;
        }
    }
/*
 * Problems; try calculating the un-weighted moments
 */
    double sigma11W;                    // quadratic moment of the weighting function

    if (iter == MAXIT || unweighted) {
        w11 = 0;                        // => unweighted
        std::pair<bool, double> moments = calcmom(*image, xCen, yCen, w11);

        if (moments.first) {
            sigma11W = moments.second;  // estimate of object moment
        } else {
            sigma11W = 1/12.0;          // a single pixel
        }
    } else {
        sigma11W = 1/w11;
    }

    return 2*sigma11W;                  // 2:  <x^2> + <y^2>
}

}

/**
 * @brief Compute the 'sigma' value for an equivalent gaussian psf.
 *
 */
double PsfAttributes::computeGaussianWidth(PsfAttributes::Method how) const {
    /*
     * Estimate the PSF's center.  This *really* needs to be rewritten to avoid using MeasureSources;
     * we shouldn't need to instantiate source objects just to measure an adaptive centroid!
     */
    afwImage::MaskedImage<double> mi = afwImage::MaskedImage<double>(_psfImage);
    typedef afwImage::Exposure<double> Exposure;
    Exposure::Ptr exposure = makeExposure(mi);
    auto foot = std::make_shared<afwDetection::Footprint>(std::make_shared<afwGeom::SpanSet>(exposure->getBBox(
        afwImage::LOCAL)));

    afwGeom::Point2D center(_psfImage->getX0() + _psfImage->getWidth()/2,
                            _psfImage->getY0() + _psfImage->getHeight()/2);
    double x(center.getX());
    double y(center.getY());
    afwGeom::Point2D fittedCenter = base::GaussianCentroidAlgorithm::fitCentroid(*_psfImage, x, y);
    double const xCen = fittedCenter.getX();
    double const yCen = fittedCenter.getY();
    switch (how) {
      case ADAPTIVE_MOMENT:
        return ::sqrt(0.5*computeSecondMomentAdaptive(_psfImage, xCen, yCen));
      case FIRST_MOMENT:
          return ::sqrt(2.0/afwGeom::PI)*computeFirstMoment(_psfImage, xCen, yCen);
      case SECOND_MOMENT:
        return ::sqrt(0.5*computeSecondMoment(_psfImage, xCen, yCen));
      case NOISE_EQUIVALENT:
        return ::sqrt(computeEffectiveArea()/(4*afwGeom::PI));
      case BICKERTON:
        {
            double sum = 0.0;
            double norm = 0.0;
            for (int iY = 0; iY != _psfImage->getHeight(); ++iY) {
                int iX = 0;
                for (afwImage::Image<double>::x_iterator ptr = _psfImage->row_begin(iY),
                                                         end = _psfImage->row_end(iY); ptr != end;
                     ++ptr, ++iX) {
                    double const x = iX - xCen;
                    double const y = iY - yCen;
                    double const r = std::sqrt( x*x + y*y );
                    double const m = (*ptr)*r;
                    norm += (*ptr)*(*ptr);
                    sum += m*m;
                }
            }
            return sqrt(sum/norm);
        }
      default:
        abort();
    }
}

/**
 * @brief Compute the effective area of the psf ( sum(I)^2/sum(I^2) )
 *
 */
double PsfAttributes::computeEffectiveArea() const {

    double sum = 0.0;
    double sumsqr = 0.0;
    for (int iY = 0; iY != _psfImage->getHeight(); ++iY) {
        afwImage::Image<double>::x_iterator end = _psfImage->row_end(iY);
        for (afwImage::Image<double>::x_iterator ptr = _psfImage->row_begin(iY); ptr != end; ++ptr) {
            sum += *ptr;
            sumsqr += (*ptr)*(*ptr);
        }
    }
    return sum*sum/sumsqr;
}

}}}
