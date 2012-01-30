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
 

/**
 * @file
 */

#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/meas/algorithms/Measure.h"
#include "lsst/meas/algorithms/CentroidControl.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace afwDet = lsst::afw::detection;
namespace afwImage = lsst::afw::image;

namespace lsst {
namespace meas {
namespace algorithms {

class SillyCentroidControl : public AlgorithmControl {
public:
    SillyCentroidControl() : AlgorithmControl("centroid.silly") {}

private:
    LSST_ALGORITHM_CONTROL_PRIVATE_DECL()
};

namespace {

/**
 * @brief A class that knows how to calculate centroids by guessing the wrong answer
 */
template<typename ExposureT>
class SillyCentroid : public Algorithm<ExposureT>
{
public:
    typedef Algorithm<ExposureT> AlgorithmT;

    SillyCentroid(SillyCentroidControl const & ctrl, afw::table::Schema & schema) :
        AlgorithmT(ctrl),
        _keys(addCentroidFields(schema, ctrl.name, "silly centroid docs"))
    {}

    virtual void apply(afw::table::SourceRecord &, ExposurePatch<ExposureT> const&) const;

private:
    afw::table::KeyTuple<afw::table::Centroid> _keys;
};

/**
 * @brief Given an image and a pixel position, return a Centroid offset by (1, 1) from initial position
 */
template<typename ExposureT>
void SillyCentroid<ExposureT>::apply(
    afw::table::SourceRecord & source,
    ExposurePatch<ExposureT> const& patch
) const {
    source.set(_keys.meas, patch.getCenter() + afw::geom::Extent2D(1, 1));
    source.set(_keys.flag, true);
}

} // anonymous

LSST_ALGORITHM_CONTROL_PRIVATE_IMPL(SillyCentroidControl, SillyCentroid)

}}}
