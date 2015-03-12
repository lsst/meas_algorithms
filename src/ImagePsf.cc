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
#include "lsst/meas/algorithms/Photometry.h"
#include "lsst/meas/base/detail/SdssShapeImpl.h"

namespace lsst { namespace meas { namespace algorithms {

double ImagePsf::doComputeApertureFlux(
    double radius, afw::geom::Point2D const & position, afw::image::Color const & color
) const {
    afw::image::MaskedImage<double> mi(computeKernelImage(position, color, INTERNAL));

    afw::geom::Point2D const center(0.0, 0.0);
    afw::geom::ellipses::Axes const axes(radius, radius);

    std::pair<double,double> result =
        photometry::calculateSincApertureFlux(mi, afw::geom::ellipses::Ellipse(axes, center));
    
    return result.first;
}

afw::geom::ellipses::Quadrupole ImagePsf::doComputeShape(
    afw::geom::Point2D const & position, afw::image::Color const & color
) const {
    base::detail::SdssShapeImpl shape;
    PTR(Image) image = computeKernelImage(position, color, INTERNAL);
    // n.b. getAdaptiveMoments doesn't account for xy0, so we have to do it manually
    base::detail::getAdaptiveMoments(
        *image,
        0.0, -image->getX0(), -image->getY0(), 1.0,   // background, x, y, shiftmax
        &shape
    );
    return afw::geom::ellipses::Quadrupole(shape.getIxx(), shape.getIyy(), shape.getIxy());
}

}}} // namespace lsst::meas::algorithms
