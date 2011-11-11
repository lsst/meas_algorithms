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
 
#if !defined(LSST_MEAS_ALGORITHMS_MEASURE_H)
#define LSST_MEAS_ALGORITHMS_MEASURE_H

//!
// Measure properties of an image selected by a Footprint
//
#include <list>
#include <cmath>
#include "lsst/base.h"
#include "boost/cstdint.hpp"
#include "boost/type_traits.hpp"
#include "boost/noncopyable.hpp"
#include "lsst/pex/logging/Log.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/image/ImageUtils.h"
#include "lsst/afw/detection.h"
#include "lsst/afw/detection/Measurement.h"
#include "lsst/afw/detection/Astrometry.h"
#include "lsst/afw/detection/Photometry.h"
#include "lsst/afw/detection/Shape.h"
#include "lsst/meas/algorithms/MeasureQuantity.h"
#include "lsst/meas/algorithms/Flags.h"
#include "lsst/meas/algorithms/detail/Measure.h"

namespace lsst {
namespace pex {
    namespace policy {
        class Policy;
    }
}
namespace afw {
    namespace detection {
        class Psf;
    }
}
namespace meas {
namespace algorithms {

/// High-level class to perform source measurement
///
/// Iterates over the various measurement types (Astrometry, Shape, Photometry).
template<typename ExposureT>
class MeasureSources : private boost::noncopyable {
public:
    typedef PTR(MeasureSources) Ptr;
    typedef CONST_PTR(MeasureSources) ConstPtr;
    typedef MeasureAstrometry<ExposureT> MeasureAstrometryT;
    typedef MeasureShape<ExposureT> MeasureShapeT;
    typedef MeasurePhotometry<ExposureT> MeasurePhotometryT;

    MeasureSources(lsst::pex::policy::Policy const& policy ///< Policy to describe processing
        );
    virtual ~MeasureSources() {}
    
    
    /// Return the Policy used to describe processing
    lsst::pex::policy::Policy const& getPolicy() const { return _policy; }
    /// Return the log
    lsst::pex::logging::Log &getLog() const { return *_moLog; }
    /// return the astrometric measurer
    typename MeasureAstrometryT::Ptr getMeasureAstrom() const { return _measureAstrom; }
    /// return the photometric measurer
    typename MeasurePhotometryT::Ptr getMeasurePhotom() const { return _measurePhotom; }
    /// return the shape measurer
    typename MeasureShapeT::Ptr getMeasureShape() const { return _measureShape; }

    /// Measure a single exposure
    virtual void measure(
        lsst::afw::detection::Source& target, ///< Input/output source: has footprint, receives measurements
        CONST_PTR(ExposureT) exp              ///< Exposure to measure
        ) const;
    virtual void measure(
        lsst::afw::detection::Source& target, ///< Input/output source: has footprint, receives measurements
        CONST_PTR(ExposureT) exp,             ///< Exposure to measure
        afwGeom::Point2D const& center        ///< Position to measure (in image frame)
        ) const;
    virtual void measure(
        lsst::afw::detection::Source& target, ///< Output source: receives measurements
        lsst::afw::detection::Source const& source, ///< Input source: has footprint and previous measurements
        lsst::afw::image::Wcs const& wcs,           ///< WCS for input source
        CONST_PTR(ExposureT) exp                    ///< Exposure to measure
        ) const;
    /// Measure multiple images
    virtual void measure(
        lsst::afw::detection::Source& target, ///< Output source: receives measurements
        lsst::afw::detection::Source const& source, ///< Input source: has footprint and previous measurements
        lsst::afw::image::Wcs const& wcs,           ///< WCS for input source
        std::vector<CONST_PTR(ExposureT)> const& exposures ///< Exposures to measure
        ) const;

private:
    lsst::pex::policy::Policy _policy;   // Policy to describe processing
    PTR(lsst::pex::logging::Log) _moLog; // log for measureObjects

    /*
     * Objects that know how to measure the object's properties
     */
    typename MeasureAstrometry<ExposureT>::Ptr _measureAstrom;
    typename MeasurePhotometry<ExposureT>::Ptr _measurePhotom;
    typename MeasureShape<ExposureT>::Ptr      _measureShape;
};

/**
 * Factory function to return a MeasureSources of the correct type (cf. std::make_pair)
 */
template<typename ExposureT>
PTR(MeasureSources<ExposureT>) makeMeasureSources(ExposureT const& exp,
                                                  lsst::pex::policy::Policy const& policy) {
    return boost::make_shared<MeasureSources<ExposureT> >(policy);
}
       
}}}
#endif // LSST_MEAS_ALGORITHMS_MEASURE_H
