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

class MeasureSources {
public:

    typedef std::list<CONST_PTR(Algorithm)> AlgorithmList;

    /// @brief Create a new MeasureSources object starting with the minimal source schema.
    MeasureSources();

    /// @brief Create a new MeasureSources object starting with the given schema.
    explicit MeasureSources(afw::table::Schema const & schema);

    /// @brief Return the schema defined by the registered algorithms.
    afw::table::Schema getSchema() const { return _schema; }

    /// @brief Set the schema.  The given schema must be a superset of the current schema.
    void setSchema(afw::table::Schema const & schema) {
        if (!schema.contains(_schema)) {
            throw LSST_EXCEPT(
                lsst::pex::exceptions::LogicErrorException,
                "New schema for MeasureSources must contain all fields in the original schema."
            );
        }
        _schema = schema;
    }

    /**
     *  @brief Add an algorithm defined by a control object.
     *
     *  Algorithms are not sorted by their 'order' control field when added.  This must be done in advance,
     *  because we would also like to guarantee that an algorithm's dependencies are registered
     *  to the schema before it is, and this registration happens when the algorithm is added here.
     */
    void addAlgorithm(AlgorithmControl const & algorithmControl) {
        _algorithms.push_back(algorithmControl.makeAlgorithm(_schema));
    }

    /**
     *  @brief Return the list of algorithms.
     *
     *  The order is the same as the order the algorithms were added in, which is also the order they
     *  will be executed.
     */
    AlgorithmList const & getAlgorithms() const { return _algorithms; }

    /**
     *  @brief Apply the registered algorithms to the given source.
     *
     *  This overload simply passes an explicit center to the algorithms (which may or may not use it).
     */
    template <typename PixelT>
    void apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    /**
     *  @brief Apply the registered algorithms to the given source.
     *
     *  This overload uses the zeroth peak in the source's footprint to determine the center passed
     *  to the algorithms.
     */
    template <typename PixelT>
    void apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure
    ) const;

private:
    afw::table::Schema _schema;
    PTR(pex::logging::Log) _log;
    AlgorithmList _algorithms;
};

}}}
#endif // LSST_MEAS_ALGORITHMS_MEASURE_H
