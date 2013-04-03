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

namespace test {
namespace foo {
namespace bar {

class SillyCentroidControl : public lsst::meas::algorithms::CentroidControl {
public:
    LSST_CONTROL_FIELD(param, int, "Difference to apply to y");

    SillyCentroidControl() : lsst::meas::algorithms::CentroidControl("centroid.silly"), param(0) {}
private:
    virtual PTR(lsst::meas::algorithms::AlgorithmControl) _clone() const { 
        return boost::make_shared<SillyCentroidControl>(*this);
    }

    virtual PTR(lsst::meas::algorithms::Algorithm) _makeAlgorithm(
        lsst::afw::table::Schema & schema,         
        PTR(lsst::daf::base::PropertyList) const & metadata,
        lsst::meas::algorithms::AlgorithmMap const & other
    ) const;

};

#ifndef SWIG
namespace {

/**
 * @brief A class that knows how to calculate centroids by guessing the wrong answer
 */
class SillyCentroid : public lsst::meas::algorithms::CentroidAlgorithm {
public:

    SillyCentroid(
        SillyCentroidControl const & ctrl,
        lsst::afw::table::Schema & schema,
        lsst::meas::algorithms::AlgorithmMap const & others
    ) : lsst::meas::algorithms::CentroidAlgorithm(ctrl, schema, "silly centroid docs"),
        _nOthers(others.size()), _param(ctrl.param)
    {}

private:
    
    template <typename PixelT>
    void _apply(
        lsst::afw::table::SourceRecord & source,
        lsst::afw::image::Exposure<PixelT> const & exposure,
        lsst::afw::geom::Point2D const & center
    ) const;
    
    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(SillyCentroid);

    int _nOthers;
    int _param;
};

/**
 * Given an image and a pixel position, return a Centroid offset by (nOthers, param) from initial position,
 * where nOthers is the number of algorithms registered before this one.  This is just a trick to
 * return something about the previously registered algorithms to the Python test code when it
 * doesn't have access to the algorithm class.
 */
template <typename PixelT>
void SillyCentroid::_apply(
    lsst::afw::table::SourceRecord & source,
    lsst::afw::image::Exposure<PixelT> const & exposure,
    lsst::afw::geom::Point2D const & center
) const {
    source.set(getKeys().meas, center + lsst::afw::geom::Extent2D(_nOthers, _param));
    source.set(getKeys().flag, false);
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(SillyCentroid);

} // anonymous

PTR(lsst::meas::algorithms::Algorithm) SillyCentroidControl::_makeAlgorithm(
    lsst::afw::table::Schema & schema, PTR(lsst::daf::base::PropertyList) const & metadata,
    lsst::meas::algorithms::AlgorithmMap const & others) const
{
    return boost::make_shared<SillyCentroid>(*this, boost::ref(schema), others);
}
#endif

}}}
