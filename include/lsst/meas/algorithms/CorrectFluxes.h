// -*- LSST-C++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2014 LSST Corporation.
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

#ifndef LSST_MEAS_ALGORITHMS_CorrectFluxes_h_INCLUDED
#define LSST_MEAS_ALGORITHMS_CorrectFluxes_h_INCLUDED

#include <vector>
#include <set>

#include "lsst/meas/algorithms/Algorithm.h"

namespace lsst { namespace meas { namespace algorithms {

typedef std::set<std::string> ApCorrRegistry;

// Return a singleton set that contains the names of all flux fields that should be aperture-corrected
ApCorrRegistry & getApCorrRegistry();

/**
 *  @brief A control object for a pluggable algorithm that scales fluxes to a common system.
 */
class CorrectFluxesControl : public AlgorithmControl {
public:

    LSST_CONTROL_FIELD(
        doFlagApCorrFailures, bool,
        "Whether to set the general failure flag for a flux when it cannot be aperture-corrected"
    );

    LSST_CONTROL_FIELD(
        doRecordApCorr, bool,
        "Whether to save the per-source per-flux aperture corrections and their errors"
    );

    LSST_CONTROL_FIELD(
        ignored, std::vector<std::string>,
        ("List of flux fields that should not be corrected (otherwise all fields in "
         "getApCorrRegistry() will be)")
    );

    CorrectFluxesControl() :
        AlgorithmControl("correctfluxes", 3.0),
        doFlagApCorrFailures(true),
        doRecordApCorr(true),
        ignored()
    {}

    PTR(CorrectFluxesControl) clone() const {
        return boost::static_pointer_cast<CorrectFluxesControl>(_clone());
    }

private:

    virtual PTR(AlgorithmControl) _clone() const;

    virtual PTR(Algorithm) _makeAlgorithm(
        afw::table::Schema & schema,
        PTR(daf::base::PropertyList) const & metadata,
        AlgorithmMap const & others,
        bool isForced
    ) const;

};

}}} // namespace lsst::meas::algorithms

#endif // !LSST_MEAS_ALGORITHMS_CorrectFluxes_h_INCLUDED
