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
 
#ifndef LSST_MEAS_ALGORITHMS_NormalizableFlux_H
#define LSST_MEAS_ALGORITHMS_NormalizableFlux_H
//!
// Control/algorithm hierarchy for flux measurement.
//

#include "lsst/base.h"
#include "lsst/pex/config.h"
#include "lsst/meas/algorithms/Algorithm.h"

namespace lsst {
namespace meas {
namespace algorithms {

/**
 *  @brief An abstract base class for flux algorithms that need to be scaled to a common system.
 *
 *  This is intended as an extra (i.e multiple-inheritance) base class for flux measurement algorithms,
 *  but it doesn't inherit from Algorithm itself so there should be no diamond hierarchy problems.
 *
 *  Algorithms that inherit from NormalizableFlux must also apply themselves to the PSF model, storing
 *  the "flux" and an error flag using keys available through the NormalizableFlux interface.  The
 *  actual fluxes and flux errors will then be corrected using this "PSF factor" measurement and
 *  the measured aperture corrections at the aperture corrections stage.
 *
 *  Inheriting from NormalizableFlux is now what marks a flux measurement as one that requires aperture
 *  corrections, unlike the previous system, in which fluxes to correct were part of the aperture
 *  correction configuration.
 */
class NormalizableFlux {
public:

    typedef std::pair< afw::table::Key<float>, afw::table::Key<afw::table::Flag> > PsfFactorKeyPair;

    /// @brief Return the number of normalizable flux measurements provided by this algorithm.
    virtual int getFluxCount() const { return 0; }

    /// @brief Return the keys that correspond to the nth flux measurement.
    virtual afw::table::KeyTuple<afw::table::Flux> getFluxKeys(int n=0) const = 0;

    /// @brief Return the keys that correspond to the PSF factor for the nth flux measurement.
    virtual PsfFactorKeyPair getPsfFactorKeys(int n=0) const = 0;

    virtual ~NormalizableFlux() {}
};

}}}// namespace lsst::meas::algorithms

#endif // !LSST_MEAS_ALGORITHMS_NormalizableFlux_H
