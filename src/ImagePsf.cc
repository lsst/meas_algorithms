// -*- LSST-C++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
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

#include "lsst/afw/image/MaskedImage.h"
#include "lsst/meas/algorithms/ImagePsf.h"
#include "lsst/meas/base/SdssShape.h"
#include "lsst/meas/base/ApertureFlux.h"

namespace lsst { namespace meas { namespace algorithms {

double ImagePsf::doComputeApertureFlux(
    double radius, afw::geom::Point2D const & position, afw::image::Color const & color
) const {
    afw::image::Image<double> const & image(*computeKernelImage(position, color, INTERNAL));

    afw::geom::Point2D const center(0.0, 0.0);
    afw::geom::ellipses::Axes const axes(radius, radius);
    base::ApertureFluxResult result = base::ApertureFluxAlgorithm::computeSincFlux(
        image,
        afw::geom::ellipses::Ellipse(axes, center),
        base::ApertureFluxControl()
    );
    return result.flux;
}

afw::geom::ellipses::Quadrupole ImagePsf::doComputeShape(
    afw::geom::Point2D const & position, afw::image::Color const & color
) const {
    PTR(Image) image = computeKernelImage(position, color, INTERNAL);
    return meas::base::SdssShapeAlgorithm::computeAdaptiveMoments(
        *image,
        afw::geom::Point2D(0.0, 0.0)  // image has origin at the center
    ).getShape();
}

}}} // namespace lsst::meas::algorithms
