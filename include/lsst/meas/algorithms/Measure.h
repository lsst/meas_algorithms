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
#include "lsst/pex/logging/Log.h"
#include "lsst/pex/policy.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/image/ImageUtils.h"
#include "lsst/afw/detection.h"
#include "lsst/meas/algorithms/Algorithm.h"

namespace lsst {
namespace afw {
    namespace detection {
        class Psf;
    }
}
namespace meas {
namespace algorithms {

/// High-level class to perform source measurement
///
template<typename ExposureT>
class MeasureSources {
public:

    explicit MeasureSources(pex::policy::Policy const& policy = pex::policy::Policy());
    
    /**
     *  Return the Policy used to describe processing
     *
     *  This no longer includes the description of which algorithms to run;
     *  those should be added using getMeasureAstrom().addAlgorithm(), etc.
     *  This is done automatically by MeasureSourcesConfig.makeMeasureSources()
     *  in Python.
     */
    pex::policy::Policy const& getPolicy() const { return _policy; }

    /// Return the log
    pex::logging::Log & getLog() const { return *_log; }

    /// Return the schema defined by the registered algorithms.
    afw::table::Schema getSchema() const { return _schema; }

    /// Set the schema.  The given schema must be a superset of the current schema.
    void setSchema(afw::table::Schema const & schema) {
        if (!schema.contains(_schema)) {
            throw LSST_EXCEPT(
                lsst::pex::exceptions::LogicErrorException,
                "New schema for MeasureSources must contain all fields in the original schema."
            );
        }
        _schema = schema;
    }

    /// Add a new algorithm and register its fields with the schema.
    void addAlgorithm(AlgorithmControl const & ctrl) {
        _algorithms.push_back(ctrl.makeAlgorithm<ExposureT>(_schema));
    }

    /**
     *  @brief Construct a source table.
     *
     *  The table's measurement slots will be set using the policy passed to MeasureSources
     *  upon construction.
     */
    PTR(afw::table::SourceTable) makeTable() const;

    /**
     *  @brief Measure a single source on a single exposure.
     *
     *  The table associated with the source must have been created by makeTable().
     */
    void apply(
        afw::table::SourceRecord & target, ///< Input/output source: has footprint, receives measurements
        CONST_PTR(ExposureT) exp              ///< Exposure to measure
    ) const;

    /**
     *  @brief Measure a single source on a single exposure.
     *
     *  The table associated with the source must have been created by makeTable().
     */
    void apply(
        afw::table::SourceRecord & target, ///< Input/output source: has footprint, receives measurements
        CONST_PTR(ExposureT) exp,          ///< Exposure to measure,
        afw::geom::Point2D const & center  ///< Initial centroid to use.
    ) const;

private:

    void _apply(
        afw::table::SourceRecord & target,
        ExposurePatch<ExposureT> & patch
    ) const;

    pex::policy::Policy _policy;   // Policy to describe processing
    PTR(pex::logging::Log) _log; // log for measureObjects
    afw::table::Schema _schema;
    afw::table::Key<double> _sgKey;
    std::list<PTR(Algorithm<ExposureT>)> _algorithms;
};

/**
 * Factory function to return a MeasureSources of the correct type (cf. std::make_pair)
 */
template<typename ExposureT>
PTR(MeasureSources<ExposureT>) makeMeasureSources(ExposureT const& exp,
                                                  pex::policy::Policy const& policy) {
    return boost::make_shared<MeasureSources<ExposureT> >(policy);
}
       
}}}
#endif // LSST_MEAS_ALGORITHMS_MEASURE_H
