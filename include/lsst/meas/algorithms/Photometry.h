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
 
#if !defined(LSST_MEAS_ALGORITHMS_PHOTOMETRY_H)
#define LSST_MEAS_ALGORITHMS_PHOTOMETRY_H 1
//!
// Utility routines for photometry
//
#include <string>
#include "lsst/base.h"

namespace lsst {
namespace afw {
    namespace image {
        template<typename T> class Image;
    }
}
namespace meas {
namespace algorithms {
namespace photometry {

template<typename MaskedImageT>
std::pair<double, double>
calculateSincApertureFlux(MaskedImageT const& mimage, double const xcen, double const ycen,
                          double const r1, double const r2,
                          double const angle=0.0, double const ellipticity=0.0
                         );
/**
 * Calculate the flux in a filled elliptical aperture
 */
template<typename MaskedImageT>
std::pair<double, double>
calculateSincApertureFlux(MaskedImageT const& mimage, ///< Image to measure
                          double const xcen,          ///< object's column position
                          double const ycen,  ///< object's row position
                          double const radius,    ///< major axis of aperture; pixels
                          double const angle, ///< angle of major axis, measured +ve from x-axis; radians
                          double const ellipticity ///< Desired ellipticity
                         )
{
    return calculateSincApertureFlux(mimage, xcen, ycen, 0.0, radius, angle, ellipticity);
}

}}}} // namespace lsst::meas::algorithms::photometry
#endif
