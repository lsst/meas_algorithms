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
#include "lsst/geom/AffineTransform.h"
#include "lsst/geom/Point.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/geom/SkyWcs.h"

namespace lsst {
namespace meas {
namespace algorithms {

/// A convenience container for the exposure, peak and footprint that will be measured.
///
/// This is more useful than a std::pair or similar.
template <typename ExposureT>
class ExposurePatch {
public:
    typedef unsigned char FlagT;  ///< Type for flags
    typedef std::shared_ptr<ExposurePatch> Ptr;
    typedef std::shared_ptr<ExposurePatch const> ConstPtr;

    /// Constructor
    ExposurePatch(std::shared_ptr<ExposureT const> exp,                   ///< Exposure of interest
                  std::shared_ptr<afw::detection::Footprint const> foot,  ///< Footprint on exposure
                  geom::Point2D const& center                 ///< Center of object on exposure
                  )
            : _exp(exp), _foot(foot), _center(center), _fromStandard(), _toStandard() {}
    ExposurePatch(std::shared_ptr<ExposureT const> exp,                       ///< Exposure of interest
                  afw::detection::Footprint const& standardFoot,  ///< Footprint on some other exposure
                  geom::Point2D const& standardCenter,            ///< Center on that other exposure
                  afw::geom::SkyWcs const& standardWcs            ///< WCS for that other exposure
                  )
            : _exp(exp) {
        afw::geom::SkyWcs const& expWcs = *exp->getWcs();
        auto const sky = standardWcs.pixelToSky(standardCenter);
        const_cast<std::shared_ptr<afw::detection::Footprint const>&>(_foot) =
                standardFoot.transform(standardWcs, expWcs, exp->getBBox());
        const_cast<geom::Point2D&>(_center) = expWcs.skyToPixel(sky);
        const_cast<geom::AffineTransform&>(_fromStandard) =
                standardWcs.linearizePixelToSky(sky) * expWcs.linearizeSkyToPixel(sky);
        const_cast<geom::AffineTransform&>(_toStandard) =
                expWcs.linearizePixelToSky(sky) * standardWcs.linearizeSkyToPixel(sky);
    }

    /// Accessors
    std::shared_ptr<ExposureT const> const getExposure() const { return _exp; }
    std::shared_ptr<afw::detection::Footprint const> const getFootprint() const { return _foot; }
    geom::Point2D const& getCenter() const { return _center; }
    geom::AffineTransform const& fromStandard() const { return _fromStandard; }
    geom::AffineTransform const& toStandard() const { return _toStandard; }

    /// Modifiers
    void setCenter(geom::Point2D const& center) { _center = center; }

private:
    std::shared_ptr<ExposureT const> const _exp;                   ///< Exposure to be measured
    std::shared_ptr<afw::detection::Footprint const> const _foot;  ///< Footprint to be measured
    geom::Point2D _center;                             ///< Center of source on exposure
    geom::AffineTransform const _fromStandard;         ///< Transform from standard WCS
    geom::AffineTransform const _toStandard;           ///< Transform to standard WCS
};

/// Factory function for ExposurePatch
template <typename ExposureT>
std::shared_ptr<ExposurePatch<ExposureT>>
makeExposurePatch(std::shared_ptr<ExposureT const> exp, std::shared_ptr<afw::detection::Footprint const> foot,
                  geom::Point2D const& center) {
    return std::make_shared<ExposurePatch<ExposureT> >(exp, foot, center);
}
template <typename ExposureT>
std::shared_ptr<ExposurePatch<ExposureT>>
makeExposurePatch(std::shared_ptr<ExposureT const> exp, afw::detection::Footprint const& standardFoot,
                  geom::Point2D const& standardCenter, afw::geom::SkyWcs const& standardWcs) {
    return std::make_shared<ExposurePatch<ExposureT> >(exp, standardFoot, standardCenter, standardWcs);
}

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst

#endif
