// -*- LSST-C++ -*-
/*
 * LSST Data Management System
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

#include "lsst/utils/ieee.h"
#include "lsst/meas/algorithms/CorrectFluxes.h"
#include "lsst/afw/detection/Psf.h"

namespace lsst { namespace meas { namespace algorithms {

namespace {

class CorrectFluxes : public Algorithm {
public:

    CorrectFluxes(
        CorrectFluxesControl const & ctrl,
        afw::table::Schema & schema,
        AlgorithmMap const & others
    );

private:

    void _apply(
        afw::table::SourceRecord & source,
        PTR(afw::detection::Psf const) psf,
        afw::geom::Point2D const & center
    ) const;

    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const {
        return this->_apply(source, exposure.getPsf(), center);
    }

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(CorrectFluxes);

};

CorrectFluxes::CorrectFluxes(
    CorrectFluxesControl const & ctrl,
    afw::table::Schema & schema,
    AlgorithmMap const & others
) : Algorithm(ctrl) {}

void CorrectFluxes::_apply(
    afw::table::SourceRecord & source,
    PTR(afw::detection::Psf const) psf,
    afw::geom::Point2D const & center
) const {
    CorrectFluxesControl const & ctrl = static_cast<CorrectFluxesControl const &>(this->getControl());
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(CorrectFluxes);

} // anonymous

PTR(AlgorithmControl) CorrectFluxesControl::_clone() const {
    return boost::make_shared<CorrectFluxesControl>(*this);
}

PTR(Algorithm) CorrectFluxesControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const & metadata,
    AlgorithmMap const & others,
    bool isForced
) const {
    return boost::make_shared<CorrectFluxes>(*this, boost::ref(schema), others);
}

}}} // namespace lsst::meas::algorithms
