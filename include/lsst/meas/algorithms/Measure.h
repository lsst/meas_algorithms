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
#include "lsst/meas/algorithms/CentroidControl.h"

namespace lsst { namespace meas { namespace algorithms {

class MeasureSources {
public:

    typedef std::list<CONST_PTR(Algorithm)> AlgorithmList;

    /// @brief Create a new MeasureSources object starting with the minimal source schema.
    MeasureSources();

    /// @brief Create a new MeasureSources object starting with the given schema.
    explicit MeasureSources(afw::table::Schema const & schema);

    /// @brief Return the schema defined by the registered algorithms.
    afw::table::Schema getSchema() const { return _schema; }

    /// @brief Set the schema.  The new schema must be a superset of the current schema.
    void setSchema(afw::table::Schema const & schema);

    /**
     *  @brief Add an algorithm defined by a control object.
     *
     *  Algorithms are not sorted by their 'order' control field when added.  This must be done in advance,
     *  because we would also like to guarantee that an algorithm's dependencies are registered
     *  to the schema before it is, and this registration happens when the algorithm is added here.
     *
     *  Algorithms may be registered multiple times, but never with the same name; this allows the same
     *  algorithm to be run multiple times with different parameters.
     */
    void addAlgorithm(AlgorithmControl const & algorithmControl);

    /**
     *  @brief Set the centroid algorithm run before all other algorithms to refine the center point.
     *
     *  The given centroid algorithm will inserted at the front of the list, and if it runs successfully
     *  its output will be used as 'center' argument passed to all subsequent algorithms.
     *
     *  If another centroid algorithm has already been set, this will be moved in front of it, 
     *  and the old algorithm will just be run like any other algorithm.
     *
     *  This also registers a flag in the schema, 'flags.badcentroid', that will be set if the
     *  centroid algorithm did not succeed and hence the center passed to subsequent algorithms was
     *  the user's center argument or peak value (depending on which overload of apply is called).
     */
    void setCentroider(CentroidControl const & centroidControl);

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
     *  This overload passes a user-defined center to the algorithms (unless setCentroider has
     *  been called, in which case only that algorithm will be passed the given center).
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
     *  to the algorithms (unless setCentroider has been called, in which case only that
     *  algorithm will use the zeroth peak).
     */
    template <typename PixelT>
    void apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure
    ) const;

private:
    afw::table::Schema _schema;
    afw::table::Key<afw::table::Flag> _badCentroidKey;
    PTR(pex::logging::Log) _log;
    AlgorithmList _algorithms;
    PTR(CentroidAlgorithm) _centroider;
};

}}}
#endif // LSST_MEAS_ALGORITHMS_MEASURE_H
