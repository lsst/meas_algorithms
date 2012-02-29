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
    ControlSet::iterator i = _ctrls.insert(p);
    return *this;
}

MeasureSourcesBuilder & MeasureSourcesBuilder::setCentroider(CentroidControl const & centroidControl) {
    PTR(CentroidControl) p = centroidControl.clone();
    p->name = _prefix + p->name;
    _centroider = centroidControl.clone();
    return *this;
}

MeasureSources MeasureSourcesBuilder::build(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const & metadata
) const {
    MeasureSources r;
    r._log.reset(pex::logging::Log::getDefaultLog().createChildLog("meas.algorithms.MeasureSource",
                                                                   pex::logging::Log::INFO));
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

template <typename PixelT>
void MeasureSources::apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center
) const {
    afw::geom::Point2D c(center);
    CONST_PTR(afw::table::SourceTable) table = source.getTable();
    for (
        AlgorithmList::const_iterator i = _algorithms.begin();
        i != _algorithms.end();
        ++i
    ) {
        try {
            (**i).apply(source, exposure, c);
            if (*i == _centroider) { // should only match the first alg, but test is cheap
                if (source.get(_centroider->getKeys().flag)) {
                    source.set(_badCentroidKey, true);
                } else {
                    c = source.get(_centroider->getKeys().meas);
                }
            }
        } catch (pex::exceptions::Exception const & e) {
            // Swallow all exceptions, because one bad measurement shouldn't affect all others
            _log->log(pex::logging::Log::DEBUG, boost::format("Measuring %s at (%d,%d): %s") %
                      (**i).getControl().name % center.getX() % center.getY() % e.what());
        } catch (...) {
            _log->log(pex::logging::Log::WARN, 
                      boost::format("Measuring %s at (%d,%d): Unknown non-LSST exception.") %
                      (**i).getControl().name % center.getX() % center.getY());
        }
    }
}

template <typename PixelT>
void MeasureSources::apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure
) const {
    CONST_PTR(afw::detection::Footprint) foot = source.getFootprint();
    afw::detection::Footprint::PeakList const& peakList = foot->getPeaks();
    if (peakList.size() == 0) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, 
                          (boost::format("No peak for source %d") % source.getId()).str());
    }
    PTR(afw::detection::Peak) peak = peakList[0];
    // set the initial centroid in the patch using the peak, then refine it.
    afw::geom::Point2D center(peak->getFx(), peak->getFy());
    apply(source, exposure, center);
}

template void MeasureSources::apply(afw::table::SourceRecord &, afw::image::Exposure<float> const &) const;
template void MeasureSources::apply(afw::table::SourceRecord &, afw::image::Exposure<double> const &) const;

template void MeasureSources::apply(
    afw::table::SourceRecord &, afw::image::Exposure<float> const &, afw::geom::Point2D const &
) const;

template void MeasureSources::apply(
    afw::table::SourceRecord &, afw::image::Exposure<double> const &, afw::geom::Point2D const &
) const;

}}} // namespace lsst::meas::algorithms
