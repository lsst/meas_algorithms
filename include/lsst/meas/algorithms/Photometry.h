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

unsigned int const SINC_COEFFS_MAX_CACHE = 64;   ///< Maximum number of coefficients to cache

/**
 * Calculate the flux in an elliptical annulus
 *
 * Caching of the coefficients is only performed for circular apertures, and
 * when explicitly allowed.  Turning caching on when the exact annulus radii
 * do not frequently recur (e.g., when the radii are determined from moments of
 * the objects) will needlessly consume memory.
 */
template<typename MaskedImageT>
std::pair<double, double>
calculateSincApertureFlux(MaskedImageT const& mimage, ///< Image to measure
                          double const xcen,          ///< Center in x
                          double const ycen,          ///< Center in y
                          double const r1,            ///< Inner radius of annulus
                          double const r2,            ///< Outer radius of annulus
                          double const angle=0.0,     ///< angle of major axis, measured +ve from x; radians
                          double const ellipticity=0.0, ///< Desired ellipticity
                          bool const allowCache=false ///< Allow caching of coefficients?
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
                          double const ellipticity, ///< Desired ellipticity
                          bool const allowCache=false ///< Allow caching of coefficients?
                         )
{
    return calculateSincApertureFlux(mimage, xcen, ycen, 0.0, radius, angle, ellipticity);
}

}}}} // namespace lsst::meas::algorithms::photometry
#endif
