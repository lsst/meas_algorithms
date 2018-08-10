// -*- lsst-c++ -*-

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

#include "lsst/geom/AffineTransform.h"
#include "lsst/geom/Box.h"
#include "lsst/pex/exceptions.h"
#include "lsst/meas/algorithms/WarpedPsf.h"
#include "lsst/afw/math/warpExposure.h"
#include "lsst/afw/image/Image.h"

namespace lsst {
namespace meas {
namespace algorithms {

namespace {

inline double min4(double a, double b, double c, double d) {
    return std::min(std::min(a, b), std::min(c, d));
}

inline double max4(double a, double b, double c, double d) {
    return std::max(std::max(a, b), std::max(c, d));
}

// TODO: make this routine externally callable and more generic using templates
//  (also useful in e.g. math/offsetImage.cc)
PTR(afw::detection::Psf::Image) zeroPadImage(afw::detection::Psf::Image const &im, int xPad, int yPad) {
    int nx = im.getWidth();
    int ny = im.getHeight();

    PTR(afw::detection::Psf::Image)
    out = std::make_shared<afw::detection::Psf::Image>(nx + 2 * xPad, ny + 2 * yPad);
    out->setXY0(im.getX0() - xPad, im.getY0() - yPad);

    geom::Box2I box(geom::Point2I(xPad, yPad), geom::Extent2I(nx, ny));
    out->assign(im, box, afw::image::LOCAL);

    return out;
}

geom::Box2I computeBBoxFromTransform(geom::Box2I const bbox, geom::AffineTransform const &t) {
    static const int dst_padding = 0;

    // This is the maximum scaling I can imagine we'd want - it's approximately what you'd
    // get from trying to coadd 2"-pixel data (e.g. 2MASS) onto a 0.01"-pixel grid (e.g.
    // from JWST).  Anything beyond that is probably a bug elsewhere in the code, and
    // catching that can save us from segfaults and catastrophic memory consumption.
    static const double maxTransformCoeff = 200.0;

    if (t.getLinear().getMatrix().lpNorm<Eigen::Infinity>() > maxTransformCoeff) {
        throw LSST_EXCEPT(pex::exceptions::RangeError, "Unexpectedly large transform passed to WarpedPsf");
    }

    // min/max coordinate values in input image
    int in_xlo = bbox.getMinX();
    int in_xhi = bbox.getMinX() + bbox.getWidth() - 1;
    int in_ylo = bbox.getMinY();
    int in_yhi = bbox.getMinY() + bbox.getHeight() - 1;

    // corners of output image
    geom::Point2D c00 = t(geom::Point2D(in_xlo, in_ylo));
    geom::Point2D c01 = t(geom::Point2D(in_xlo, in_yhi));
    geom::Point2D c10 = t(geom::Point2D(in_xhi, in_ylo));
    geom::Point2D c11 = t(geom::Point2D(in_xhi, in_yhi));

    //
    // bounding box for output image
    //
    int out_xlo = floor(min4(c00.getX(), c01.getX(), c10.getX(), c11.getX())) - dst_padding;
    int out_ylo = floor(min4(c00.getY(), c01.getY(), c10.getY(), c11.getY())) - dst_padding;
    int out_xhi = ceil(max4(c00.getX(), c01.getX(), c10.getX(), c11.getX())) + dst_padding;
    int out_yhi = ceil(max4(c00.getY(), c01.getY(), c10.getY(), c11.getY())) + dst_padding;

    geom::Box2I ret = geom::Box2I(geom::Point2I(out_xlo, out_ylo),
                                  geom::Extent2I(out_xhi - out_xlo + 1, out_yhi - out_ylo + 1));
    return ret;
}

/**
 * @brief Alternate interface to afw::math::warpImage()
 * in which the caller does not need to precompute the output bounding box.
 *
 * This version takes an affine transform instead of an arbitrary xy transform.
 *
 * @param[in] im  Image to warp
 * @param[in] srcToDest  Affine transformation from source pixels to destination pixels in the forward
 *                  direction; the warping code only uses the inverse direction
 * @param[in] wc  Warping parameters
 *
 * The input image is assumed zero-padded.
 */
PTR(afw::detection::Psf::Image)
warpAffine(afw::detection::Psf::Image const &im, geom::AffineTransform const &srcToDest,
           afw::math::WarpingControl const &wc) {
    std::shared_ptr<afw::geom::TransformPoint2ToPoint2> srcToDestTransform =
            afw::geom::makeTransform(srcToDest);

    afw::math::SeparableKernel const &kernel = *wc.getWarpingKernel();
    geom::Point2I const &center = kernel.getCtr();
    int const xPad = std::max(center.getX(), kernel.getWidth() - center.getX());
    int const yPad = std::max(center.getY(), kernel.getHeight() - center.getY());

    // allocate output image
    geom::Box2I bbox = computeBBoxFromTransform(im.getBBox(), srcToDest);
    PTR(afw::detection::Psf::Image) ret = std::make_shared<afw::detection::Psf::Image>(bbox);

    // zero-pad input image
    PTR(afw::detection::Psf::Image) im_padded = zeroPadImage(im, xPad, yPad);

    // warp it!
    afw::math::warpImage(*ret, *im_padded, *srcToDestTransform, wc, 0.0);
    return ret;
}

}  // namespace

WarpedPsf::WarpedPsf(PTR(afw::detection::Psf const) undistortedPsf,
                     PTR(afw::geom::TransformPoint2ToPoint2 const) distortion,
                     CONST_PTR(afw::math::WarpingControl) control)
        : ImagePsf(false),
          _undistortedPsf(undistortedPsf),
          _distortion(distortion),
          _warpingControl(control) {
    _init();
}

WarpedPsf::WarpedPsf(PTR(afw::detection::Psf const) undistortedPsf,
                     PTR(afw::geom::TransformPoint2ToPoint2 const) distortion, std::string const &kernelName,
                     unsigned int cache)
        : ImagePsf(false),
          _undistortedPsf(undistortedPsf),
          _distortion(distortion),
          _warpingControl(new afw::math::WarpingControl(kernelName, "", cache)) {
    _init();
}

void WarpedPsf::_init() {
    if (!_undistortedPsf) {
        throw LSST_EXCEPT(pex::exceptions::LogicError,
                          "Undistorted Psf passed to WarpedPsf must not be None/NULL");
    }
    if (!_distortion) {
        throw LSST_EXCEPT(pex::exceptions::LogicError, "Transform passed to WarpedPsf must not be None/NULL");
    }
    if (!_warpingControl) {
        throw LSST_EXCEPT(pex::exceptions::LogicError,
                          "WarpingControl passed to WarpedPsf must not be None/NULL");
    }
}

geom::Point2D WarpedPsf::getAveragePosition() const {
    return _distortion->applyForward(_undistortedPsf->getAveragePosition());
}

PTR(afw::detection::Psf) WarpedPsf::clone() const {
    return std::make_shared<WarpedPsf>(_undistortedPsf->clone(), _distortion, _warpingControl);
}

PTR(afw::detection::Psf) WarpedPsf::resized(int width, int height) const {
    // For a given set of requested dimensions and distortion, it is not guaranteed that a
    // _undistortedPsf would exist to manifest those dimensions after distortion
    // Not possible to implement with member data currently in WarpedPsf
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Not Implemented");
}

PTR(afw::detection::Psf::Image)
WarpedPsf::doComputeKernelImage(geom::Point2D const &position, afw::image::Color const &color) const {
    geom::AffineTransform t = afw::geom::linearizeTransform(*_distortion->inverted(), position);
    geom::Point2D tp = t(position);

    PTR(Image) im = _undistortedPsf->computeKernelImage(tp, color);

    // Go to the warped coordinate system with 'p' at the origin
    auto srcToDest = geom::AffineTransform(t.inverted().getLinear());
    PTR(afw::detection::Psf::Psf::Image) ret = warpAffine(*im, srcToDest, *_warpingControl);

    double normFactor = 1.0;
    //
    // Normalize the output image to sum 1
    // FIXME defining a member function Image::getSum() would be convenient here and in other places
    //
    normFactor = 0.0;
    for (int y = 0; y != ret->getHeight(); ++y) {
        afw::detection::Psf::Image::x_iterator imEnd = ret->row_end(y);
        for (afw::detection::Psf::Image::x_iterator imPtr = ret->row_begin(y); imPtr != imEnd; imPtr++) {
            normFactor += *imPtr;
        }
    }
    if (normFactor == 0.0) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "psf image has sum 0");
    }
    *ret /= normFactor;
    return ret;
}

geom::Box2I WarpedPsf::doComputeBBox(geom::Point2D const &position, afw::image::Color const &color) const {
    geom::AffineTransform t = afw::geom::linearizeTransform(*_distortion->inverted(), position);
    geom::Point2D tp = t(position);
    geom::Box2I bboxUndistorted = _undistortedPsf->computeBBox(tp, color);
    geom::Box2I ret =
            computeBBoxFromTransform(bboxUndistorted, geom::AffineTransform(t.inverted().getLinear()));
    return ret;
}

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst
