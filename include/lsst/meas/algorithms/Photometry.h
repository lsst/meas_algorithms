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
#include "lsst/afw/image/MaskedImage.h"

namespace lsst {
namespace afw {
namespace geom {
namespace ellipses {
    class Ellipse;
}}}
namespace meas {
namespace algorithms {
namespace photometry {

/*
 * A comparison function that doesn't require equality closer than machine epsilon
 */
template <typename T>
struct fuzzyCompare {
    bool operator()(T x, T y) const {
        if (isEqual(x, y)) {
            return false;
        }
        return (x - y < 0) ? true : false;
    }
    bool isEqual(T x, T y) const {
        return ::fabs(x - y) < std::numeric_limits<T>::epsilon();
    }
};

/// A singleton to calculate and cache the coefficients for sinc photometry
///
/// Caching is only performed for circular apertures (because elliptical
/// apertures are assumed to be generated dynamically, and hence not expected
/// to recur).  Caching must be explicitly requested for a particular circular
/// aperture (using the 'cache' method).
///
/// Unless otherwise noted, apertures are elliptical annuli.  The outer boundary
/// is defined by an ellipse, while the inner boundary is a scaled version of
/// the outer ellipse.
template<typename PixelT>
class SincCoeffs {
public:
    typedef afw::image::Image<PixelT> CoeffT;

    /// Cache the coefficients for a particular aperture
    ///
    /// The aperture is a circular annulus.
    static void cache(float r1,         ///< Inner radius
                      float r2          ///< Outer radius
        );

    /// Get the coefficients for an aperture
    ///
    /// Coefficients are retrieved from the cache, if available; otherwise they
    /// will be generated.
    static CONST_PTR(CoeffT)
    get(afw::geom::ellipses::Axes const& axes, ///< Ellipse defining outer boundary
        float const innerFactor=0.0            ///< Scale to apply to ellipse for inner boundary
        );

    /// Calculate the coefficients for an aperture
    static PTR(CoeffT)
    calculate(afw::geom::ellipses::Axes const& axes, ///< Ellipse defining outer boundary
              double const innerFactor=0.0           ///< Scale to apply to ellipse for inner boundary
        );

private:
    typedef std::map<float, PTR(CoeffT), fuzzyCompare<float> > CoeffMap;
    typedef std::map<float, CoeffMap, fuzzyCompare<float> > CoeffMapMap;
    SincCoeffs() : _cache() {};
    SincCoeffs(SincCoeffs const&); // unimplemented: singleton
    void operator=(SincCoeffs const&); // unimplemented: singleton
    static SincCoeffs& getInstance();

    /// Search the cache for coefficients for an aperture
    ///
    /// If the coefficients are not cached, a null shared_ptr will be returned.
    CONST_PTR(CoeffT)
    _lookup(afw::geom::ellipses::Axes const& axes, double const innerFactor=0.0) const;

    CoeffMapMap _cache;                 ///< Cache of coefficients
};


/**
 * Calculate the flux in an elliptical annulus
 *
 * The outer boundary is defined by an ellipse, while the inner boundary is a
 * scaled version of the outer ellipse.
 */
template<typename MaskedImageT>
std::pair<double, double>
calculateSincApertureFlux(MaskedImageT const& mimage, ///< Image to measure
                          afw::geom::ellipses::Ellipse const& ellipse, ///< Outer aperture
                          double const innerFactor=0.0 ///< Scale to apply to ellipse for inner boundary
                         );

}}}} // namespace lsst::meas::algorithms::photometry
#endif
