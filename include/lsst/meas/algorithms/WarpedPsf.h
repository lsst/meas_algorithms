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

#include "lsst/afw/geom/Point.h"
#include "lsst/afw/geom/Extent.h"
#include "lsst/afw/geom/XYTransform.h"
#include "lsst/meas/algorithms/ImagePsf.h"

#ifndef LSST_AFW_DETECTION_WARPEDPSF_H
#define LSST_AFW_DETECTION_WARPEDPSF_H

namespace lsst { namespace meas { namespace algorithms {

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
     * convention for the transform is that p' = distortion.forwardTransform(p)
     */
    WarpedPsf(CONST_PTR(afw::detection::Psf) undistortedPsf, CONST_PTR(afw::geom::XYTransform) distortion);

    virtual afw::geom::Point2D getAveragePosition() const;

    virtual PTR(afw::detection::Psf) clone() const;

protected:

    virtual PTR(afw::detection::Psf::Image) doComputeKernelImage(
        afw::geom::Point2D const & position, afw::image::Color const & color
    ) const;

protected:
    CONST_PTR(afw::detection::Psf) _undistortedPsf;
    CONST_PTR(afw::geom::XYTransform) _distortion;
};

}}} // namespace lsst::meas::algorithms

#endif  // LSST_MEAS_ALGORITHMS_WARPEDPSF_H
