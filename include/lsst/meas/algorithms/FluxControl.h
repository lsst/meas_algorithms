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
 
#ifndef LSST_MEAS_ALGORITHMS_FLUXCONTROL_H
#define LSST_MEAS_ALGORITHMS_FLUXCONTROL_H
//!
// Control (and secretly, factory) object hierarchy for flux algorithms.
//

#include "lsst/base.h"
#include "lsst/pex/config.h"
#include "lsst/meas/algorithms/Algorithm.h"

namespace lsst {
namespace meas {
namespace algorithms {

/**
 *  @brief C++ control object for aperture flux.
 *
 *  @sa ApertureFluxConfig.
 */
class ApertureFluxControl : public AlgorithmControl {
public:

    LSST_CONTROL_FIELD(radii, std::vector<double>, "vector of radii for apertures (in pixels)");

    ApertureFluxControl() : AlgorithmControl("flux.aperture"), radii() {}

private:
    LSST_ALGORITHM_CONTROL_PRIVATE_DECL()
};

/**
 *  @brief C++ control object for Gaussian flux.
 *
 *  @sa GaussianFluxConfig.
 */
class GaussianFluxControl : public AlgorithmControl {
public:

    LSST_CONTROL_FIELD(fixed, bool,
                       "if true, use existing shape and centroid measurements instead of fitting");
    LSST_CONTROL_FIELD(background, double, "FIXME! NEVER DOCUMENTED!");
    LSST_CONTROL_FIELD(shiftmax, double, "FIXME! NEVER DOCUMENTED!");
    LSST_CONTROL_FIELD(centroid, std::string, "name of centroid field to use if fixed is true");
    LSST_CONTROL_FIELD(shape, std::string, "name of shape field to use if fixed is true");

    GaussianFluxControl() : 
        AlgorithmControl("flux.gaussian"), fixed(false), background(0.0), shiftmax(10.0),
        centroid("shape.sdss.centroid"), shape("shape.sdss")
    {}

private:
    LSST_ALGORITHM_CONTROL_PRIVATE_DECL()
};

/**
 *  @brief C++ control object for naive flux.
 *
 *  @sa NaiveFluxConfig.
 */
class NaiveFluxControl : public AlgorithmControl {
public:

    LSST_CONTROL_FIELD(radius, double, "FIXME! NEVER DOCUMENTED!");

    NaiveFluxControl() : AlgorithmControl("flux.naive"), radius(9.0) {}

private:
    LSST_ALGORITHM_CONTROL_PRIVATE_DECL()
};

/**
 *  @brief C++ control object for PSF flux.
 *
 *  @sa PsfFluxConfig.
 */
class PsfFluxControl : public AlgorithmControl {
public:

    PsfFluxControl() : AlgorithmControl("flux.psf") {}

private:
    LSST_ALGORITHM_CONTROL_PRIVATE_DECL()
};

/**
 *  @brief C++ control object for sinc aperture flux.
 *
 *  @sa SincFluxConfig.
 */
class SincFluxControl : public AlgorithmControl {
public:

    LSST_CONTROL_FIELD(radius1, double, "FIXME! NEVER DOCUMENTED!");
    LSST_CONTROL_FIELD(radius2, double, "FIXME! NEVER DOCUMENTED!");
    LSST_CONTROL_FIELD(angle, double, "FIXME! NEVER DOCUMENTED!");
    LSST_CONTROL_FIELD(ellipticity, double, "FIXME! NEVER DOCUMENTED!");

    SincFluxControl() : 
        AlgorithmControl("flux.sinc"), radius1(0.0), radius2(0.0), angle(0.0), ellipticity(0.0) {}

private:
    LSST_ALGORITHM_CONTROL_PRIVATE_DECL()
};

}}}// namespace lsst::meas::algorithms

#endif // !LSST_MEAS_ALGORITHMS_FLUXCONTROL_H
