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

#ifndef LSST_MEAS_ALGORITHMS_CorrectFluxes_h_INCLUDED
#define LSST_MEAS_ALGORITHMS_CorrectFluxes_h_INCLUDED

#include "lsst/meas/algorithms/Algorithm.h"
#include "lsst/meas/algorithms/ScaledFlux.h"

namespace lsst { namespace meas { namespace algorithms {

/**
 *  @brief A control object for a pluggable algorithm that scales fluxes to a common system.
 *
 *  Unlike most pluggable algorithms, CorrectFluxes mostly operates on existing fields, instead
 *  of creating its own.  Upon construction, it extracts Keys from all ScaledFlux algorithms
 *  that have been registered, one of which is considered the "canonical flux", which all the
 *  ScaledFluxes will be tied to for point sources.
 *
 *  When applied to a source, in its default configuration CorrectFluxes will:
 *   - ties all the ScaledFluxes to the "canonical" flux by multiplying them by the canonical flux's
 *     "PSF factor" (see ScaledFlux) divided by their own PSF factor;
 *   - computes an aperture correction as the result of measuring that aperture on the PSF model
 *     at the location of the source, divided by the canonical PSF factor;
 *   - applies the aperture correction to all ScaledFluxes, by multiplying each flux by the
 *     aperture correction.
 *
 *  @note Generally, the canonical flux will be the PSF flux, which is just the dot product of the PSF
 *  model with the data.  If the PSF model is properly normalized, the PSF factor for the canonical
 *  flux should thus be exactly one, but we do not assume this here.
 */
class CorrectFluxesControl : public AlgorithmControl {
public:

    LSST_CONTROL_FIELD(
        doApCorr, bool,
        "Whether to compute and apply the aperture correction from the PSF model"
    );
    LSST_CONTROL_FIELD(
        doFlagApCorrFailures, bool,
        "Whether to set the general failure flag for a flux when it cannot be aperture-corrected"
    );
    LSST_CONTROL_FIELD(
        doTieScaledFluxes, bool,
        "Whether to tie all ScaledFluxes to the canonical flux"
    );
    LSST_CONTROL_FIELD(
        doFlagTieFailures, bool,
        "Whether to set the general failure flag for a flux when it cannot be tied to the canonical flux"
    );
    LSST_CONTROL_FIELD(
        apCorrRadius, double,
        "Aperture to correct fluxes to in pixels"
    );
    LSST_CONTROL_FIELD(
        canonicalFluxName, std::string,
        "Name of algorithm to directly aperture-correct and tie other ScaledFluxes to"
    );
    LSST_CONTROL_FIELD(
        canonicalFluxIndex, int,
        "Index of canonical flux within canonical flux algorithm"
    );

    CorrectFluxesControl() :
        AlgorithmControl("correctfluxes", 3.0),
        doApCorr(true),
        doFlagApCorrFailures(true),
        doTieScaledFluxes(true),
        doFlagTieFailures(true),
        apCorrRadius(7.0),
        canonicalFluxName("flux.psf"),
        canonicalFluxIndex(0)
    {}

    PTR(CorrectFluxesControl) clone() const {
        return boost::static_pointer_cast<CorrectFluxesControl>(_clone());
    }

private:

    virtual PTR(AlgorithmControl) _clone() const;

    virtual PTR(Algorithm) _makeAlgorithm(
        afw::table::Schema & schema,
        PTR(daf::base::PropertyList) const & metadata,
        AlgorithmMap const & others
    ) const;

};

}}} // namespace lsst::meas::algorithms

#endif // !LSST_MEAS_ALGORITHMS_CorrectFluxes_h_INCLUDED
