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
 
#ifndef LSST_MEAS_ALGORITHMS_ASTROMETRYCONTROL_H
#define LSST_MEAS_ALGORITHMS_ASTROMETRYCONTROL_H
//!
// Control (and secretly, factory) object hierarchy for astrometry algorithms.
//

#include "lsst/base.h"
#include "lsst/pex/config.h"
#include "lsst/meas/algorithms/Algorithm.h"

namespace lsst {
namespace meas {
namespace algorithms {

typedef AlgorithmControl<afw::detection::Astrometry> AstrometryControl;

/**
 *  @brief C++ control object for Gaussian astrometry.
 *
 *  @sa GaussianAstrometryConfig.
 */
class GaussianAstrometryControl : public AstrometryControl {
public:

    GaussianAstrometryControl() {}

private:
    LSST_ALGORITHM_CONTROL_PRIVATE_DECL()
};

/**
 *  @brief C++ control object for naive astrometry.
 *
 *  @sa NaiveAstrometryConfig.
 */
class NaiveAstrometryControl : public AstrometryControl {
public:

    LSST_CONTROL_FIELD(background, double, "FIXME! NEVER DOCUMENTED!");

    NaiveAstrometryControl() : background(0.0) {}

private:
    LSST_ALGORITHM_CONTROL_PRIVATE_DECL()
};

/**
 *  @brief C++ control object for SDSS astrometry.
 *
 *  @sa SdssAstrometryConfig.
 */
class SdssAstrometryControl : public AstrometryControl {
public:

    LSST_CONTROL_FIELD(binmax, int, "FIXME! NEVER DOCUMENTED!");
    LSST_CONTROL_FIELD(peakMin, double, "FIXME! NEVER DOCUMENTED!");
    LSST_CONTROL_FIELD(wfac, double, "FIXME! NEVER DOCUMENTED!");

    SdssAstrometryControl() : binmax(16), peakMin(-1.0), wfac(1.5) {}

private:
    LSST_ALGORITHM_CONTROL_PRIVATE_DECL()
};

}}}// namespace lsst::meas::algorithms

#endif // !LSST_MEAS_ALGORITHMS_ASTROMETRYCONTROL_H
