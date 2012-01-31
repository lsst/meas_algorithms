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
#include <set>
#include <cmath>

#include "boost/cstdint.hpp"
#include "boost/type_traits.hpp"
#include "boost/array.hpp"

#include "lsst/base.h"
#include "lsst/pex/logging/Log.h"
#include "lsst/pex/config.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/image/ImageUtils.h"
#include "lsst/afw/detection.h"
#include "lsst/meas/algorithms/Algorithm.h"

namespace lsst { namespace meas { namespace algorithms {

struct SourceSlotControl {
    LSST_CONTROL_FIELD(
        centroid, std::string, "the name of the centroiding algorithm used to set source (x,y)"
    );
    LSST_CONTROL_FIELD(
        shape, std::string, "the name of the algorithm used to set source ellipticity and size parameters"
    );
    LSST_CONTROL_FIELD(
        apFlux, std::string, "the name of the algorithm used to set the source aperture flux slot"
    );
    LSST_CONTROL_FIELD(
        modelFlux, std::string, "the name of the algorithm used to set the source model flux"
    );
    LSST_CONTROL_FIELD(
        psfFlux, std::string, "the name of the algorithm used to set the source PSF flux"
    );
    LSST_CONTROL_FIELD(
        instFlux, std::string, "the name of the algorithm used to set the source instrumental flux"
    );

    void apply(PTR(afw::table::SourceTable) const & table) const;

    SourceSlotControl() :
        centroid("centroid.sdss"), shape("shape.sdss"),
        apFlux("flux.sinc"), modelFlux("flux.gaussian"), psfFlux("flux.psf"), instFlux("flux.gaussian")
    {}

};

struct ClassificationControl {
    LSST_CONTROL_FIELD(sg_fac1, double, "First S/G parameter; critical ratio of inst to psf flux");
    LSST_CONTROL_FIELD(sg_fac2, double, "Second S/G parameter; correction for instFlux error");
    LSST_CONTROL_FIELD(sg_fac3, double, "Third S/G parameter; correction for psfFlux error");

    ClassificationControl() :
        sg_fac1(0.925), sg_fac2(0.0), sg_fac3(0.0)
    {}
};

// Flags are now indices to bits, not bitmasks.
class MeasureSources {
public:

    enum FlagEnum {
        EDGE=0,
        INTERPOLATED,
        INTERPOLATED_CENTER,
        SATURATED,
        SATURATED_CENTER,
        PEAKCENTER,
        N_MEASURE_SOURCES_FLAGS
    };

    typedef boost::array< afw::table::Key<afw::table::Flag>, N_MEASURE_SOURCES_FLAGS > FlagKeyArray;

    MeasureSources();

    explicit MeasureSources(Schema const & schema);

    /// Return the schema keys for the flags set by MeasureSources itself (not its algorithms).
    FlagKeyArray const & getFlagKeys() const { return _flagKeys; }

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

    /// Add an algorithm defined by a control object.
    void addAlgorithm(AlgorithmControl const & algorithmControl) {
        _controlSet.insert(algorithmControl.clone());
    }

    template <typename ExposureT>
    void apply(afw::table::SourceVector const & sources, ExposureT const & exposure) const;

private:

    class CompareAlgorithmControl {
    public:
        bool operator()(PTR(AlgorithmControl) const & a, PTR(AlgorithmControl) const & b) const {
            return a->order < b->order;
        }
    };

    typedef std::multiset<PTR(AlgorithmControl),CompareAlgorithmControl> AlgorithmControlSet;

    void _initialize();

    afw::table::Schema _schema;
    FlagKeyArray _flagKeys;
    AlgorithmControlSet _controlSet;
};

}}}
#endif // LSST_MEAS_ALGORITHMS_MEASURE_H
