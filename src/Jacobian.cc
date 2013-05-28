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
#include "lsst/meas/algorithms/Jacobian.h"


namespace lsst {
namespace meas {
namespace algorithms {

namespace {

/// A class that knows how to calculate Jacobians
class Jacobian : public algorithms::Algorithm {
public:

    Jacobian(JacobianControl const& ctrl, afw::table::Schema& schema) :
        algorithms::Algorithm(ctrl),
        _key(schema.addField<double>(ctrl.name, "Jacobian correction")),
        _scale(::pow(ctrl.pixelScale, -2))
        {}

private:

    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(Jacobian);

    afw::table::Key<double> _key;
    float _scale;                       // Scale to apply to pixel area
};

template <typename PixelT>
void Jacobian::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const& exposure,
    afw::geom::Point2D const & center
) const {
    source.set(_key,
               ::fabs(_scale*exposure.getWcs()->linearizePixelToSky(center, afw::geom::arcseconds).getLinear().
                      computeDeterminant()));
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(Jacobian);

} // anonymous namespace


PTR(algorithms::Algorithm) JacobianControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const &
) const {
    return boost::make_shared<Jacobian>(*this, boost::ref(schema));
}


}}} // namespace lsst::meas::algorithms
