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
 
#if !defined(LSST_MEAS_ALGORITHMS_ALGORITHM_H)
#define LSST_MEAS_ALGORITHMS_ALGORITHM_H

#include "boost/noncopyable.hpp"

#include "lsst/base.h"
#include "lsst/pex/logging/Log.h"
#include "lsst/pex/policy.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/detection/Source.h"
#include "lsst/afw/image/Filter.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/geom/AffineTransform.h"

#include "lsst/meas/algorithms/Flags.h"

namespace pexLog = lsst::pex::logging;
namespace pexPolicy = lsst::pex::policy;
namespace afwDet = lsst::afw::detection;
namespace afwImage = lsst::afw::image;
namespace afwGeom = lsst::afw::geom;
namespace afwCoord = lsst::afw::coord;

namespace lsst { namespace meas { namespace algorithms {

/// A convenience container for the exposure, peak and footprint that will be measured.
///
/// This is more useful than a std::pair or similar.
template<typename ExposureT>
class ExposurePatch : private boost::noncopyable {
public:
    typedef unsigned char FlagT;
    typedef PTR(ExposurePatch) Ptr;
    typedef CONST_PTR(ExposurePatch) ConstPtr;

    /// Flag values, indicating which measurement is bad
    enum { NONE           = 0x00,     /// None bad
           EDGE           = 0x01,     /// Footprint overlaps an edge
           INTERP         = 0x02,     /// Footprint includes interpolated pixels
           INTERP_CENTER  = 0x04,     /// Peak pixel is interpolated
           SAT            = 0x08,     /// Footprint includes saturated pixels
           SAT_CENTER     = 0x10,     /// Peak pixel is saturated
           ASTROMETRY     = 0x20,     /// Bad astrometry
           SHAPE          = 0x40,     /// Bad shapes
           PHOTOMETRY     = 0x80,     /// Bad photometry
           ALL            = 0xFF      /// All are bad
    };

    /// Constructor
    explicit ExposurePatch(CONST_PTR(ExposureT) exp,
                           CONST_PTR(afwDet::Footprint) foot,
                           afwGeom::Point2D const& center
        ): _exp(exp), _foot(foot), _center(center), _fromStandard(), _toStandard(), _flags(NONE) {}
    explicit ExposurePatch(CONST_PTR(ExposureT) exp,
                           CONST_PTR(afwDet::Footprint) foot,
                           afwGeom::Point2D const& standardCenter,
                           afwImage::Wcs const& standardWcs
        ) : _exp(exp), _foot(foot), _center(), _fromStandard(), _toStandard(), _flags(NONE) {
        afwImage::Wcs const& expWcs = *exp->getWcs();
        afwCoord::Coord::ConstPtr sky = standardWcs.pixelToSky(standardCenter);
        const_cast<afwGeom::Point2D&>(_center) = expWcs.skyToPixel(sky);
        const_cast<afwGeom::AffineTransform&>(_fromStandard) = standardWcs.linearizePixelToSky(sky) *
            expWcs.linearizeSkyToPixel(sky);
        const_cast<afwGeom::AffineTransform&>(_toStandard) = expWcs.linearizePixelToSky(sky) *
            standardWcs.linearizeSkyToPixel(sky);
    }

    /// Accessors
    CONST_PTR(ExposureT) const getExposure() const { return _exp; }
    CONST_PTR(afwDet::Footprint) const getFootprint() const { return _foot; }
    afwGeom::Point2D const& getCenter() const { return _center; }
    afwGeom::AffineTransform const& fromStandard() const { return _fromStandard; }
    afwGeom::AffineTransform const& toStandard() const { return _toStandard; }
    bool getFlags() const { return _flags; }

    /// Modifiers
    void setFlags(FlagT flags) { _flags = flags; }
    void orFlag(FlagT flags) { _flags |= flags; }
    void setCenter(afwGeom::Point2D const& center) { _center = center; }

    /// Flag translator
    ///
    /// Converts ExposurePatch flags to Source flags
    static boost::int64_t sourceFlags(int epFlags) {
        boost::int64_t sFlags = 0;
        if (epFlags & EDGE) { sFlags |= Flags::EDGE; }
        if (epFlags & INTERP) { sFlags |= Flags::INTERP; }
        if (epFlags & INTERP_CENTER) { sFlags |= Flags::INTERP_CENTER; }
        if (epFlags & SAT) { sFlags |= Flags::SATUR; }
        return sFlags;
    }

private:
    CONST_PTR(ExposureT) const _exp;          // Exposure to be measured
    CONST_PTR(afwDet::Footprint) const _foot; // Footprint to be measured
    afwGeom::Point2D _center;           // Center of source on exposure
    afwGeom::AffineTransform const _fromStandard; // Transform from standard WCS
    afwGeom::AffineTransform const _toStandard; // Transform to standard WCS
    FlagT _flags;                       // Flags indicating which measurement is bad
};

/// Factory function for ExposurePatch
template<typename ExposureT>
PTR(ExposurePatch<ExposureT>) makeExposurePatch(CONST_PTR(ExposureT) exp, CONST_PTR(afwDet::Footprint) foot,
                                                afwGeom::Point2D const& center) {
    return boost::make_shared<ExposurePatch<ExposureT> >(exp, foot, center);
}
template<typename ExposureT>
PTR(ExposurePatch<ExposureT>) makeExposurePatch(CONST_PTR(ExposureT) exp, CONST_PTR(afwDet::Footprint) foot,
                                                afwGeom::Point2D const& standardCenter,
                                                afwImage::Wcs const& standardWcs) {
    return boost::make_shared<ExposurePatch<ExposureT> >(exp, foot, standardCenter, standardWcs);
}

/// Base class for algorithms for measuring MeasurementT (e.g., Photometry)
template<typename MeasurementT, typename ExposureT>
class Algorithm : private boost::noncopyable {
public:
    /// Constructor
    Algorithm() {}

    /// Destructor
    virtual ~Algorithm() {}

    /// Measure a single value from a single image.
    ///
    /// Returns leaf of MeasurementT (single measurement).
    ///
    /// Pure-virtual, so subclasses MUST define: it is the essence of the
    /// measurement, as the other measure functions can (but need not) be
    /// defined from it.
    virtual PTR(MeasurementT) measureSingle(afwDet::Source const& target,
                                            afwDet::Source const& source,
                                            ExposurePatch<ExposureT> const& patch) const = 0;
    
    /// Measure a single value from multiple images.
    ///
    /// Returns leaf of MeasurementT (single measurement).
    ///
    /// Because it is a 'group' of images (images in the same filter), we can
    /// assume they share some characteristics (e.g., center, shape).
    ///
    /// Defaults to iteratively calling measureOne. However, if the measurement
    /// cannot be obtained by merely averaging the outputs of a single
    /// measurement, e.g., some measured parameters are made across all
    /// exposures as part of the measurement (e.g., a common shape), then the
    /// Algorithm needs to define this method.
    virtual PTR(MeasurementT) measureMultiple(afwDet::Source const& target,
                                              afwDet::Source const& source,
                                              std::vector<CONST_PTR(ExposurePatch<ExposureT>)> const& patches
        ) const {
        typedef std::vector<CONST_PTR(ExposurePatch<ExposureT>)> PatchVector;
        PTR(MeasurementT) meas(new MeasurementT());
        for (typename PatchVector::const_iterator iter = patches.begin(); iter != patches.end(); ++iter) {
            PTR(MeasurementT) val;
            try {
                CONST_PTR(ExposurePatch<ExposureT>) patch = *iter;
                val = measureSingle(target, source, *patch);
            } catch (lsst::pex::exceptions::Exception const& e) {
#if 0
                std::cerr << (boost::format("Measuring single %s at (%d,%d): %s") %
                              getName() % source.getXAstrom() % source.getYAstrom() %
                              e.what()).str() << std::endl;
#endif
                val = measureNull();
            }
            val->setAlgorithm(getName());
            meas->add(val);
        }
        return meas->average();
    }

    /// Return a null measurement
    ///
    /// This is called when we hit an exception.
    virtual PTR(MeasurementT) measureNull(void) const {
        return MeasurementT::null();
    }
    
    /// Configure the algorithm
    virtual void configure(pexPolicy::Policy const&) {};

    /// Name of the algorithm
    virtual std::string getName() const = 0;

    /// Clone algorithm
    virtual PTR(Algorithm<MeasurementT, ExposureT>) clone() const = 0;
};


}}} // namespace

#endif
