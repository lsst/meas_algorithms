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

#include "lsst/afw/image.h"
#include "lsst/afw/table.h"
#include "lsst/meas/algorithms/CentroidControl.h"
#include "lsst/meas/algorithms/RecordCentroid.h"

namespace lsst {
namespace meas {
namespace algorithms {

namespace {

/**
 * @brief A class that records the centroid used in measurement.
 * 
 * Note that this is not a real CentroidAlgorithm, and hence cannot be used in a "slot" (i.e., getX(), etc).
 * In order to use in a "slot", we would need to also record a covariance matrix, which it's not clear how to
 * populate, and hence it would be a waste of space.
 *
 */
class RecordCentroid : public Algorithm {
public:

    RecordCentroid(RecordCentroidControl const & ctrl, afw::table::Schema & schema) :
        Algorithm(ctrl),
        _key(schema.addField<afw::table::Centroid::MeasTag>("centroid.record",
                                                            "simply record centroid used for measurement",
                                                            "pixels"))
        {}

private:

    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    afw::table::Key<afw::table::Point<double> > _key;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(RecordCentroid);
};

template<typename PixelT>
void RecordCentroid::_apply(
    afw::table::SourceRecord & source, 
    afw::image::Exposure<PixelT> const&,
    afw::geom::Point2D const& center
) const {
    source.set(_key, center);
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(RecordCentroid);

} // anonymous

PTR(AlgorithmControl) RecordCentroidControl::_clone() const {
    return boost::make_shared<RecordCentroidControl>(*this);
}

PTR(Algorithm) RecordCentroidControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const &
) const {
    return boost::make_shared<RecordCentroid>(*this, boost::ref(schema));
}

}}}  // namespace lsst::meas::algorithms
