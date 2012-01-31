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
 
#if !defined(LSST_MEAS_ALGORITHMS_EXPOSURE_PATCH_H)
#define LSST_MEAS_ALGORITHMS_EXPOSURE_PATCH_H

#include "lsst/base.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/geom/AffineTransform.h"

namespace lsst { namespace meas { namespace algorithms {

/// A convenience container for the exposure, peak and footprint that will be measured.
///
/// This is more useful than a std::pair or similar.
template<typename ExposureT>
class ExposurePatch {
public:
    typedef unsigned char FlagT;        ///< Type for flags
    typedef PTR(ExposurePatch) Ptr;
    typedef CONST_PTR(ExposurePatch) ConstPtr;

    /// Constructor
    ExposurePatch(CONST_PTR(ExposureT) exp, ///< Exposure of interest
                  CONST_PTR(afw::detection::Footprint) foot, ///< Footprint on exposure
                  afw::geom::Point2D const& center           ///< Center of object on exposure
        ): _exp(exp), _foot(foot), _center(center), _fromStandard(), _toStandard() {}
    ExposurePatch(CONST_PTR(ExposureT) exp, ///< Exposure of interest
                  afw::detection::Footprint const& standardFoot, ///< Footprint on some other exposure
                  afw::geom::Point2D const& standardCenter,  ///< Center on that other exposure
                  afw::image::Wcs const& standardWcs         ///< WCS for that other exposure
        ) : _exp(exp) {
        afw::image::Wcs const& expWcs = *exp->getWcs();
        afw::coord::Coord::ConstPtr sky = standardWcs.pixelToSky(standardCenter);
        const_cast<CONST_PTR(afw::detection::Footprint)&>(_foot) = standardFoot.transform(standardWcs, expWcs,
                                                                                          exp->getBBox());
        const_cast<afw::geom::Point2D&>(_center) = expWcs.skyToPixel(sky);
        const_cast<afw::geom::AffineTransform&>(_fromStandard) = standardWcs.linearizePixelToSky(sky) *
            expWcs.linearizeSkyToPixel(sky);
        const_cast<afw::geom::AffineTransform&>(_toStandard) = expWcs.linearizePixelToSky(sky) *
            standardWcs.linearizeSkyToPixel(sky);
    }

    /// Accessors
    CONST_PTR(ExposureT) const getExposure() const { return _exp; }
    CONST_PTR(afw::detection::Footprint) const getFootprint() const { return _foot; }
    afw::geom::Point2D const& getCenter() const { return _center; }
    afw::geom::AffineTransform const& fromStandard() const { return _fromStandard; }
    afw::geom::AffineTransform const& toStandard() const { return _toStandard; }

    /// Modifiers
    void setCenter(afw::geom::Point2D const& center) { _center = center; }

private:
    CONST_PTR(ExposureT) const _exp;    ///< Exposure to be measured
    CONST_PTR(afw::detection::Footprint) const _foot; ///< Footprint to be measured
    afw::geom::Point2D _center;   ///< Center of source on exposure
    afw::geom::AffineTransform const _fromStandard; ///< Transform from standard WCS
    afw::geom::AffineTransform const _toStandard; ///< Transform to standard WCS
};

/// Factory function for ExposurePatch
template<typename ExposureT>
PTR(ExposurePatch<ExposureT>) makeExposurePatch(
    CONST_PTR(ExposureT) exp,
    CONST_PTR(afw::detection::Footprint) foot,
    afw::geom::Point2D const& center
    ) {
    return boost::make_shared<ExposurePatch<ExposureT> >(exp, foot, center);
}
template<typename ExposureT>
PTR(ExposurePatch<ExposureT>) makeExposurePatch(
    CONST_PTR(ExposureT) exp,
    afw::detection::Footprint const& standardFoot,
    afw::geom::Point2D const& standardCenter,
    afw::image::Wcs const& standardWcs
    ) {
    return boost::make_shared<ExposurePatch<ExposureT> >(exp, standardFoot, standardCenter, standardWcs);
}

}}} // namespace

#endif
