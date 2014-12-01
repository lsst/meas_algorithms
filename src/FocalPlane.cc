/*
 * Lsst Data Management System
 * Copyright 2008-2013 LSST Corporation.
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

// -*- LSST-C++ -*-
#include <cmath>
#include "lsst/afw/cameraGeom/Detector.h"
#include "lsst/meas/algorithms/FocalPlane.h"


namespace lsst {
namespace meas {
namespace algorithms {

namespace {

/// A class that knows how to calculate focal plane position
class FocalPlane : public algorithms::Algorithm {
public:
    typedef afw::table::Point<double> FieldT;

    FocalPlane(FocalPlaneControl const& ctrl, afw::table::Schema& schema) :
        algorithms::Algorithm(ctrl),
        _key(schema.addField<FieldT>(ctrl.name, "Focal plane position"))
        {}

private:

    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(FocalPlane);

    afw::table::Key<FieldT> _key;
};

template <typename PixelT>
void FocalPlane::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const& exposure,
    afw::geom::Point2D const & center
) const {
    CONST_PTR(afw::cameraGeom::Detector) det = exposure.getDetector();
    afw::geom::Point2D fp;
    double const NaN = std::numeric_limits<double>::quiet_NaN();
    if (!det) {
        fp = afw::geom::Point2D(NaN, NaN);
    } else {
        fp = det->getPositionFromPixel(center).getMm();
    }
    source.set(_key, fp);
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(FocalPlane);

} // anonymous namespace


PTR(algorithms::Algorithm) FocalPlaneControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const &
) const {
    return boost::make_shared<FocalPlane>(*this, boost::ref(schema));
}


}}} // namespace lsst::meas::algorithms
