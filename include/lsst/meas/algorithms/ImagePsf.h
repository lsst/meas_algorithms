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
#ifndef LSST_MEAS_ALGORITHMS_ImagePsf_h_INCLUDED
#define LSST_MEAS_ALGORITHMS_ImagePsf_h_INCLUDED

#include "lsst/geom/Point.h"
#include "lsst/afw/detection/Psf.h"

namespace lsst {
namespace meas {
namespace algorithms {

/**
 *  @brief An intermediate base class for Psfs that use an image representation.
 *
 *  ImagePsf exists only to provide implementations of doComputeApertureFlux and doComputeShape
 *  for its derived classes.  These implementations use the SincFlux and SdssShape algorithms
 *  defined in meas_algorithms, and hence could not be included with the Psf base class in afw.
 */
class ImagePsf : public afw::table::io::PersistableFacade<ImagePsf>, public afw::detection::Psf {
protected:
    explicit ImagePsf(bool isFixed = false) : afw::detection::Psf(isFixed) {}

    virtual double doComputeApertureFlux(double radius, geom::Point2D const& position,
                                         afw::image::Color const& color) const;

    virtual afw::geom::ellipses::Quadrupole doComputeShape(geom::Point2D const& position,
                                                           afw::image::Color const& color) const;
};

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst

#endif  // !LSST_MEAS_ALGORITHMS_ImagePsf_h_INCLUDED
