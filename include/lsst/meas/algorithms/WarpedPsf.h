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

#include "lsst/geom/Box.h"
#include "lsst/afw/geom/Transform.h"
#include "lsst/afw/math/warpExposure.h"
#include "lsst/meas/algorithms/ImagePsf.h"

#ifndef LSST_AFW_DETECTION_WARPEDPSF_H
#define LSST_AFW_DETECTION_WARPEDPSF_H

namespace lsst {
namespace meas {
namespace algorithms {

/**
 * @brief A Psf class that maps an arbitrary Psf through a coordinate transformation
 *
 * If K_0(x,x') is the unwarped PSF, and f is the coordinate transform, then the
 * warped PSF is defined by
 *
 *   K(f(x),f(x')) = K_0(x,x')      (*)
 *
 * We linearize the coordinate transform in the vicinity of the point where the
 * PSF is computed.  The definition (*) does not include the Jacobian of the
 * transformation, since the afw convention is that PSF's are normalized to
 * have integral 1 anyway.
 */
class WarpedPsf : public ImagePsf {
public:
    /**
     * @brief Construct WarpedPsf from unwarped psf and distortion.
     *
     * If p is the nominal pixel position, and p' is the true position on the sky, then our
     * convention for the transform is that p' = distortion.applyForward(p)
     */
    WarpedPsf(CONST_PTR(afw::detection::Psf) undistortedPsf,
              CONST_PTR(afw::geom::TransformPoint2ToPoint2) distortion,
              CONST_PTR(afw::math::WarpingControl) control);
    WarpedPsf(CONST_PTR(afw::detection::Psf) undistortedPsf,
              CONST_PTR(afw::geom::TransformPoint2ToPoint2) distortion,
              std::string const& kernelName = "lanczos3", unsigned int cache = 10000);

    /**
     *  @brief Return the average of the positions of the stars that went into this Psf.
     *
     *  For WarpedPsf, this is just the transform of the undistorted Psf's average position.
     */
    virtual geom::Point2D getAveragePosition() const;

    /// Polymorphic deep copy.  Usually unnecessary, as Psfs are immutable.
    virtual PTR(afw::detection::Psf) clone() const;

    /// Return a clone with specified kernel dimensions
    virtual PTR(afw::detection::Psf) resized(int width, int height) const;

protected:
    virtual PTR(afw::detection::Psf::Image)
            doComputeKernelImage(geom::Point2D const& position, afw::image::Color const& color) const;

protected:
    PTR(afw::detection::Psf const) _undistortedPsf;
    PTR(afw::geom::TransformPoint2ToPoint2 const) _distortion;

private:
    void _init();
    CONST_PTR(afw::math::WarpingControl) _warpingControl;

    virtual geom::Box2I doComputeBBox(geom::Point2D const& position, afw::image::Color const& color) const;
};

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst

#endif  // LSST_MEAS_ALGORITHMS_WARPEDPSF_H
