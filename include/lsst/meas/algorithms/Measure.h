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
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/image/ImageUtils.h"
#include "lsst/afw/detection.h"
#include "lsst/afw/detection/Measurement.h"
#include "lsst/afw/detection/Astrometry.h"
#include "lsst/afw/detection/Photometry.h"
#include "lsst/afw/detection/Shape.h"
#include "lsst/meas/algorithms/Flags.h"

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
class MeasureSources {
public:

    MeasureSources(pex::policy::Policy const& policy);

    virtual ~MeasureSources() {}
    
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
    pex::logging::Log & getLog() const { return *_moLog; }

    /// Measure a single exposure
    virtual void measure(
        afw::table::SourceRecord & target, ///< Input/output source: has footprint, receives measurements
        CONST_PTR(ExposureT) exp              ///< Exposure to measure
    ) const;

private:
    pex::policy::Policy _policy;   // Policy to describe processing
    PTR(pex::logging::Log) _moLog; // log for measureObjects
    
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
