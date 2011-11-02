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
 
#if !defined(LSST_MEAS_ALGORITHMS_DETAIL_MEASURE_H)
#define LSST_MEAS_ALGORITHMS_DETAIL_MEASURE_H

#include <list>
#include <cmath>
#include "lsst/base.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/detection.h"
#include "lsst/afw/detection/Measurement.h"
#include "lsst/afw/detection/Astrometry.h"
#include "lsst/afw/detection/Photometry.h"
#include "lsst/afw/detection/Shape.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/geom/AffineTransform.h"
#include "lsst/meas/algorithms/MeasureQuantity.h"
#include "lsst/meas/algorithms/Flags.h"
#include "lsst/meas/algorithms/Algorithm.h"

namespace lsst {

namespace pex {
    namespace policy {
        class Policy;
    }
}

namespace meas {
namespace algorithms {

// Forward declarations
template<typename MeasurementT, typename ExposureT> class MeasureQuantity;
template<typename ExposureT> class FootprintCentroid;
template<typename MaskedImageT> class FootprintBits;

namespace detail {

namespace pexExceptions = lsst::pex::exceptions;
namespace pexPolicy = lsst::pex::policy;
namespace pexLog = lsst::pex::logging;
namespace afwDet = lsst::afw::detection;
namespace afwGeom = lsst::afw::geom;

template<typename ExposureT>
void checkFootprint(ExposurePatch<ExposureT>& patch, 
                    typename ExposureT::MaskedImageT::Mask::Pixel const bits // Bits in footprint
    );

/// Measuring sources on the same image on which they were detected.
template<typename ExposureT>
struct SingleMeasurer {
    typedef ExposurePatch<ExposureT> ExposureContainerT;
    typedef ExposurePatch<ExposureT> const ConstExposureContainerT;
    /// Check pixels in the footprint, setting appropriate flags, and get a rough starting x,y position
    static void check(afwDet::Source& source, ExposurePatch<ExposureT>& patch) {
        FootprintCentroid<typename ExposureT::MaskedImageT> centroider(patch.getExposure()->getMaskedImage());
        centroider.apply(*patch.getFootprint());
        double const x = centroider.getX();
        double const y = centroider.getY();
        source.setXAstrom(x);
        source.setYAstrom(y);
        patch.setCenter(afwGeom::Point2D(x, y));
        checkFootprint(patch, centroider.getBits());
    }

    /// Make the exposure container carry const members
    static ExposurePatch<ExposureT> const& constify(ExposurePatch<ExposureT> const& patch) {
        /// No change needed in this case
        return patch;
    }

    /// Make the measurement
    template<typename MeasurementT>
    static PTR(MeasurementT) measure(CONST_PTR(MeasureQuantity<MeasurementT, ExposureT>) measurer,
                                     afwDet::Source& target, afwDet::Source const& source,
                                     ExposurePatch<ExposureT> const& patch) {
        return measurer->measureSingle(target, source, patch);
    }

    /// Execute the algorithm
    template<typename MeasurementT>
    static PTR(MeasurementT) algorithm(CONST_PTR(Algorithm<MeasurementT, ExposureT>) alg,
                                       afwDet::Source const& target, afwDet::Source const& source,
                                       ExposurePatch<ExposureT> const& patch) {
        return alg->measureSingle(target, source, patch);
    }                                       

    /// Update the astrometry
    static void updateAstrom(afwDet::Source const& target, ExposurePatch<ExposureT>& patch) {
        afwGeom::Point2D const center(target.getXAstrom(), target.getYAstrom());
        patch.setCenter(patch.fromStandard()(center));
    }

    /// Get combined flags
    static boost::int64_t flags(ExposurePatch<ExposureT> const& patch) {
        return patch.getFlags();
    }
};

/// Measuring a single source on multiple images
template<typename ExposureT>
struct MultipleMeasurer {
    typedef std::vector<PTR(ExposurePatch<ExposureT>)> ExposureContainerT;
    typedef std::vector<CONST_PTR(ExposurePatch<ExposureT>)> const ConstExposureContainerT;
    static void check(afwDet::Source& source, 
                      std::vector<PTR(ExposurePatch<ExposureT>)>& patches) {
        for (size_t i = 0; i < patches.size(); ++i) {
            PTR(ExposurePatch<ExposureT>) p = patches[i];
            FootprintBits<typename ExposureT::MaskedImageT> bitsFunctor(p->getExposure()->getMaskedImage());
            bitsFunctor.apply(*p->getFootprint());
            checkFootprint(*p, bitsFunctor.getBits());
        }
    }
    static std::vector<CONST_PTR(ExposurePatch<ExposureT>)> const&
    constify(std::vector<PTR(ExposurePatch<ExposureT>)> const& patches) {
        /// The compiler doesn't know how to automatically convert
        /// std::vector<PTR(T)> to std::vector<CONST_PTR(T)> because the way the
        /// template system works means that in theory the two may be
        /// specialised differently.  This is an explicit conversion.
        ///
        /// see e.g., http://stackoverflow.com/questions/2102244/vector-and-const
        return reinterpret_cast<std::vector<CONST_PTR(ExposurePatch<ExposureT>)> const&>(patches);
    }
    template<typename MeasurementT>
    static PTR(MeasurementT) measure(CONST_PTR(MeasureQuantity<MeasurementT, ExposureT>) measurer,
                                     afwDet::Source& target, afwDet::Source const& source,
                                     std::vector<CONST_PTR(ExposurePatch<ExposureT>)> const& patches) {
        return measurer->measureMultiple(target, source, patches);
    }
    
    template<typename MeasurementT>
    static PTR(MeasurementT) algorithm(CONST_PTR(Algorithm<MeasurementT, ExposureT>) alg,
                                       afwDet::Source const& target, afwDet::Source const& source,
                                       std::vector<CONST_PTR(ExposurePatch<ExposureT>)> const& patches) {
        return alg->measureMultiple(target, source, patches);
    }                                       

    static void updateAstrom(afwDet::Source const& target,
                             std::vector<PTR(ExposurePatch<ExposureT>)>& patches) {
        afwGeom::Point2D const center(target.getXAstrom(), target.getYAstrom());
        for (size_t i = 0; i < patches.size(); ++i) {
            patches[i]->setCenter(patches[i]->fromStandard()(center));
        }
    }

    static boost::int64_t flags(std::vector<CONST_PTR(ExposurePatch<ExposureT>)> const& patches) {
        boost::int64_t eFlags = ExposurePatch<ExposureT>::ALL;
        for (typename ConstExposureContainerT::const_iterator iter = patches.begin(); 
             iter != patches.end(); ++iter) {
            eFlags &= (*iter)->getFlags();
        }
        return eFlags;
    }
};



}}}} // lsst::meas::algorithms::detail
#endif
