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

class MeasureSourcesBuilder;

class MeasureSources {
public:
    typedef MeasureSourcesBuilder Builder;
    typedef std::list<CONST_PTR(Algorithm)> AlgorithmList;

    /**
     *  @brief Return the list of algorithms.
     *
     *  The order is the same as the order the algorithms will be executed.  The special centroider,
     *  if present, will always be the first item in this list.
     */
    AlgorithmList const & getAlgorithms() const { return _algorithms; }

    /**
     *  @brief Apply the registered algorithms to the given source.
     *
     *  This overload passes a user-defined center to the algorithms.  If refineCenter is true
     *  and setCentroider has been called, that algorithm will be used to refine the centroid.
     *
     *  If the centroider fails (by setting its flag or returning NaNs), the original center
     *  will be used for subsequent measurements and 'flags.badcentroid' will be set.
     */
    template <typename PixelT>
    void apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center,
        bool refineCenter=true
    ) const;

    /**
     *  @brief Apply the registered algorithms to the given source.
     *
     *  This overload uses the peak position as the center for the algorithms; if setCentroider
     *  has been called, that algorithm will refine the centroid.
     *
     *  If the centroider fails (by setting its flag or returning NaNs), the peak
     *  will be used for subsequent measurements and 'flags.badcentroid' will be set.
     */
    template <typename PixelT>
    void apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure
    ) const;

    /**
     *  @brief Apply the registered algorithms to the given source.
     *
     *  This version of apply() uses the sky coordinates as the center for the algorithms;
     *  no refinement of the centroid is made.
     */
    template <typename PixelT>
    void applyWithCoord(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure
    ) const {
        apply(source, exposure, exposure.getWcs()->skyToPixel(source.getCoord()), false);
    }

    /**
     *  @brief Apply the registered algorithms to the given source.
     *
     *  This version of apply() uses a previous centroid as the center for the algorithms;
     *  no refinement of the centroid is made.
     */
    template <typename PixelT>
    void applyWithPixel(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure
    ) const {
        apply(source, exposure, afw::geom::Point2D(source.getX(), source.getY()), false);
    }

    /**
     *  @brief Apply the registered algorithms to the given source.
     *
     *  This overload uses a reference Source to provide the center and footprint.  This is
     *  intended for forced photometry (measuring a Source on one image based on its detection
     *  in another).
     */
    template <typename PixelT>
    void apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::table::SourceRecord const& reference,
        CONST_PTR(afw::image::Wcs) referenceWcs = CONST_PTR(afw::image::Wcs)()
    ) const;

private:

    MeasureSources() {}

    friend class MeasureSourcesBuilder;

    afw::table::Key<afw::table::Flag> _badCentroidKey;
    PTR(pex::logging::Log) _log;
    AlgorithmList _algorithms;
    PTR(CentroidAlgorithm) _centroider;
};

class MeasureSourcesBuilder {
public:

    /**
     *  @brief Add an algorithm defined by a control object.
     *
     *  Algorithms may be registered multiple times, but never with the same name; this allows the same
     *  algorithm to be run multiple times with different parameters.
     */
    MeasureSourcesBuilder & addAlgorithm(AlgorithmControl const & algorithmControl);

    /**
     *  @brief Set the centroid algorithm run before all other algorithms to refine the center point.
     *
     *  The given centroid algorithm will inserted at the front of the list, and if it runs successfully
     *  its output will be used as 'center' argument passed to all subsequent algorithms.
     *
     *  This also registers a flag in the schema, 'flags.badcentroid', that will be set if the
     *  centroid algorithm did not succeed and hence the center passed to subsequent algorithms was
     *  the user's center argument or peak value (depending on which overload of apply is called).
     */
    MeasureSourcesBuilder & setCentroider(CentroidControl const & centroidControl);
    
    /**
     *  @brief Build a MeasureSources object.
     *
     *  @param[in,out] schema              Schema where algorithm fields will be registered.
     *  @param[in,out] metadata            Algorithms may put table-level descriptive information that
     *                                     describes their configuration here.
     */
    MeasureSources build(
        afw::table::Schema & schema,
        PTR(daf::base::PropertyList) const & metadata = PTR(daf::base::PropertyList)()
    ) const;

    explicit MeasureSourcesBuilder(std::string const & prefix = "") : _prefix(prefix) {}

private:

    struct ComparePriority {
        bool operator()(CONST_PTR(AlgorithmControl) const & a, CONST_PTR(AlgorithmControl) const & b) const {
            return (a->priority < b->priority) ||
                (a->priority == b->priority && a->name < b->name);
        }
    };

    typedef std::set<CONST_PTR(AlgorithmControl),ComparePriority> ControlSet;
    
    std::string _prefix;
    CONST_PTR(CentroidControl) _centroider;
    ControlSet _ctrls;
};

}}}
#endif // LSST_MEAS_ALGORITHMS_MEASURE_H
