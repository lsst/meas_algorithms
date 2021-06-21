// -*- lsst-c++ -*-

/*
 * This file is part of meas_algorithms.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
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
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef LSST_MEAS_ALGORITHMS_WARPEDPSF_H
#define LSST_MEAS_ALGORITHMS_WARPEDPSF_H

#include "lsst/geom/Box.h"
#include "lsst/afw/geom/Transform.h"
#include "lsst/afw/math/warpExposure.h"
#include "lsst/meas/algorithms/ImagePsf.h"

namespace lsst {
namespace meas {
namespace algorithms {

/**
 * A Psf class that maps an arbitrary Psf through a coordinate transformation
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
     * Construct WarpedPsf from unwarped psf and distortion.
     *
     * If p is the nominal pixel position, and p' is the true position on the sky, then our
     * convention for the transform is that p' = distortion.applyForward(p)
     */
    WarpedPsf(std::shared_ptr<afw::detection::Psf const> undistortedPsf,
              std::shared_ptr<afw::geom::TransformPoint2ToPoint2 const> distortion,
              std::shared_ptr<afw::math::WarpingControl const> control);
    WarpedPsf(std::shared_ptr<afw::detection::Psf const> undistortedPsf,
              std::shared_ptr<afw::geom::TransformPoint2ToPoint2 const> distortion,
              std::string const& kernelName = "lanczos3", unsigned int cache = 10000);

    /**
     *  Return the average of the positions of the stars that went into this Psf.
     *
     *  For WarpedPsf, this is just the transform of the undistorted Psf's average position.
     */
    geom::Point2D getAveragePosition() const override;

    /// Polymorphic deep copy.  Usually unnecessary, as Psfs are immutable.
    std::shared_ptr<afw::detection::Psf> clone() const override;

    /// Return a clone with specified kernel dimensions
    std::shared_ptr<afw::detection::Psf> resized(int width, int height) const override;

protected:
    std::shared_ptr<afw::detection::Psf::Image> doComputeKernelImage(
            geom::Point2D const& position, afw::image::Color const& color) const override;

    std::shared_ptr<afw::detection::Psf const> _undistortedPsf;
    std::shared_ptr<afw::geom::TransformPoint2ToPoint2 const> _distortion;

private:
    void _init();
    std::shared_ptr<afw::math::WarpingControl const> _warpingControl;

    geom::Box2I doComputeBBox(geom::Point2D const& position, afw::image::Color const& color) const override;
};

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst

#endif  // LSST_MEAS_ALGORITHMS_WARPEDPSF_H
