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

#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Log.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/algorithms/Measure.h"

namespace lsst {
namespace meas {
namespace algorithms {

MeasureSourcesBuilder & MeasureSourcesBuilder::addAlgorithm(AlgorithmControl const & algorithmControl) {
    PTR(AlgorithmControl) p = algorithmControl.clone();
    p->name = _prefix + p->name;
    _ctrls.insert(p);
    return *this;
}

MeasureSourcesBuilder & MeasureSourcesBuilder::setCentroider(CentroidControl const & centroidControl) {
    PTR(CentroidControl) p = centroidControl.clone();
    p->name = _prefix + p->name;
    _centroider = p;
    return *this;
}

MeasureSources MeasureSourcesBuilder::build(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const & metadata
) const {
    MeasureSources r;
    r._log.reset(pex::logging::Log::getDefaultLog().createChildLog("meas.algorithms.MeasureSource", pex::logging::Log::INFO));
    if (_centroider) {
        r._badCentroidKey = schema.addField<afw::table::Flag>(
            _prefix + "flags.badcentroid",
            "the centroid algorithm used to feed centers to other algorithms failed"
        );
        r._centroider = _centroider->makeAlgorithm(schema, metadata);
        r._algorithms.push_back(r._centroider);
    }
    for (ControlSet::const_iterator i = _ctrls.begin(); i != _ctrls.end(); ++i) {
        r._algorithms.push_back((**i).makeAlgorithm(schema, metadata));
    }
    return r;
}

namespace {

/// Apply algorithm to measure source on image
///
/// Common code for executing an algorithm and catching exceptions.
template<typename PixelT>
void applyAlgorithm(
    Algorithm const& algorithm,         ///< Algorithm to apply
    afw::table::SourceRecord & source,  ///< Source to measure
    afw::image::Exposure<PixelT> const & exposure, ///< Exposure on which to measure
    afw::geom::Point2D const& center,              ///< Center to use
    PTR(pex::logging::Log) log                     ///< Log for errors
    )
{
    try {
        algorithm.apply(source, exposure, center);
    } catch (pex::exceptions::Exception const& e) {
        // Swallow all exceptions, because one bad measurement shouldn't affect all others
        log->log(pex::logging::Log::DEBUG, boost::format("Measuring %s at (%d,%d): %s") %
                 algorithm.getControl().name % center.getX() % center.getY() % e.what());
    } catch (...) {
        log->log(pex::logging::Log::WARN, 
                 boost::format("Measuring %s at (%d,%d): Unknown non-LSST exception.") %
                 algorithm.getControl().name % center.getX() % center.getY());
    }
}

} // namespace

template <typename PixelT>
void MeasureSources::apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center,
    bool refineCenter
) const {
    afw::geom::Point2D c(center);
    CONST_PTR(afw::table::SourceTable) table = source.getTable();
    for (
        AlgorithmList::const_iterator i = _algorithms.begin();
        i != _algorithms.end();
        ++i
    ) {
        applyAlgorithm(**i, source, exposure, c, _log);
        if (refineCenter && *i == _centroider) { // should only match the first alg, but test is cheap
            if (source.get(_centroider->getKeys().flag)) {
                source.set(_badCentroidKey, true);
            } else {
                c = source.get(_centroider->getKeys().meas);
            }
        }
    }
}

template <typename PixelT>
void MeasureSources::apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure
) const {
    CONST_PTR(afw::detection::Footprint) foot = source.getFootprint();
    if (!foot) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, 
                          (boost::format("No footprint for source %d") % source.getId()).str());
    }
    afw::detection::Footprint::PeakList const& peakList = foot->getPeaks();
    if (peakList.size() == 0) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, 
                          (boost::format("No peak for source %d") % source.getId()).str());
    }
    PTR(afw::detection::Peak) peak = peakList[0];
    // set the initial centroid in the patch using the peak, then refine it.
    afw::geom::Point2D center(peak->getFx(), peak->getFy());
    apply(source, exposure, center, true);
}


template <typename PixelT>
void MeasureSources::apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::table::SourceRecord const& reference,
    CONST_PTR(afw::image::Wcs) referenceWcs
) const {
    CONST_PTR(afw::image::Wcs) wcs = exposure.getWcs();

    if (referenceWcs) {
        CONST_PTR(afw::detection::Footprint) refFoot = reference.getFootprint();
        if (!refFoot) {
            throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, 
                              (boost::format("No footprint for reference %d") % reference.getId()).str());
        }
        source.setFootprint(refFoot->transform(*referenceWcs, *wcs, exposure.getBBox()));
    }

    source.set(afw::table::SourceTable::getCoordKey(), reference.getCoord());
    applyWithCoord(source, exposure);
}

#define INSTANTIATE(TYPE) \
template void MeasureSources::apply(afw::table::SourceRecord &, afw::image::Exposure<TYPE> const &) const; \
template void MeasureSources::apply(afw::table::SourceRecord &, afw::image::Exposure<TYPE> const &, \
                                    afw::geom::Point2D const &, bool) const; \
template void MeasureSources::applyWithCoord(afw::table::SourceRecord &, \
                                             afw::image::Exposure<TYPE> const &) const; \
template void MeasureSources::applyWithPixel(afw::table::SourceRecord &, \
                                             afw::image::Exposure<TYPE> const &) const; \
template void MeasureSources::apply(afw::table::SourceRecord &, afw::image::Exposure<TYPE> const &, \
                                    afw::table::SourceRecord const&, CONST_PTR(afw::image::Wcs)) const;

INSTANTIATE(float);
INSTANTIATE(double);

#undef INSTANTIATE

}}} // namespace lsst::meas::algorithms
