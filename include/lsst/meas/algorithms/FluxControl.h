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
// Control/algorithm hierarchy for flux measurement.
//

#include "lsst/base.h"
#include "lsst/pex/config.h"
#include "lsst/meas/algorithms/Algorithm.h"

namespace lsst {
namespace meas {
namespace algorithms {

class FluxControl;

/**
 *  @brief Intermediate base class for algorithms that compute a flux.
 */
class FluxAlgorithm : public Algorithm {
public:

    /**
     *  @brief Tuple type that holds the keys that define a standard flux algorithm.
     *
     *  Algorithms are encouraged to add additional flags as appropriate, but these are required.
     */
    typedef afw::table::KeyTuple<afw::table::Flux> KeyTuple;

    /// @copydoc Algorithm::getControl
    FluxControl const & getControl() const;

    /// @brief Return the standard flux keys registered by this algorithm.
    KeyTuple const & getKeys() const { return _keys; }

protected:

    /// @brief Initialize with a manually-constructed key tuple.
    FluxAlgorithm(FluxControl const & ctrl, KeyTuple const & keys);

    /// @brief Initialize using afw::table::addFlux field to fill out repetitive descriptions.
    FluxAlgorithm(FluxControl const & ctrl, afw::table::Schema & schema, char const * doc);

private:
    KeyTuple _keys;
};

class FluxControl : public AlgorithmControl {
public:

    PTR(FluxControl) clone() const { return boost::static_pointer_cast<FluxControl>(_clone()); }

    PTR(FluxAlgorithm) makeAlgorithm(
        afw::table::Schema & schema,
        PTR(daf::base::PropertyList) const & metadata = PTR(daf::base::PropertyList)(),
        AlgorithmControlMap const & others = AlgorithmControlMap(),
        bool isForced = false
    ) const {
        return boost::static_pointer_cast<FluxAlgorithm>(_makeAlgorithm(schema, metadata, others, isForced));
    }

protected:
    explicit FluxControl(std::string const & name_, double priority=2.0) :
        AlgorithmControl(name_, priority) {}
};

inline FluxAlgorithm::FluxAlgorithm(FluxControl const & ctrl, KeyTuple const & keys) :
    Algorithm(ctrl), _keys(keys)
{}

inline FluxAlgorithm::FluxAlgorithm(
    FluxControl const & ctrl, afw::table::Schema & schema, char const * doc
) :
    Algorithm(ctrl), _keys(afw::table::addFluxFields(schema, ctrl.name, doc))
{}

inline FluxControl const & FluxAlgorithm::getControl() const {
    return static_cast<FluxControl const &>(Algorithm::getControl());
}

/**
 *  @brief C++ control object for aperture flux.
 *
 *  Does not inherit from flux control because it measures an array of fluxes, rather
 *  than a single one; we could have another intermediate base class for that.
 *
 *  @sa ApertureFluxConfig.
 */
class ApertureFluxControl : public AlgorithmControl {
public:

    LSST_CONTROL_FIELD(radii, std::vector<double>, "vector of radii for apertures (in pixels)");

    ApertureFluxControl() : AlgorithmControl("flux.aperture", 2.0), radii() {}

private:
    virtual PTR(AlgorithmControl) _clone() const;
    virtual PTR(Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};

/**
 *  @brief C++ control object for Gaussian flux.
 *
 *  @sa GaussianFluxConfig.
 */
class GaussianFluxControl : public FluxControl {
public:

    LSST_CONTROL_FIELD(fixed, bool,
                       "if true, use existing shape and centroid measurements instead of fitting");
    LSST_CONTROL_FIELD(background, double, "FIXME! NEVER DOCUMENTED!");
    LSST_CONTROL_FIELD(shiftmax, double, "FIXME! NEVER DOCUMENTED!");
    LSST_CONTROL_FIELD(centroid, std::string, "name of centroid field to use if fixed is true");
    LSST_CONTROL_FIELD(shape, std::string, "name of shape field to use if fixed is true");

    GaussianFluxControl() : 
        FluxControl("flux.gaussian"), fixed(false), background(0.0), shiftmax(10.0),
        centroid("shape.sdss.centroid"), shape("shape.sdss")
    {}

private:
    virtual PTR(AlgorithmControl) _clone() const;
    virtual PTR(Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};

/**
 *  @brief C++ control object for naive flux.
 *
 *  @sa NaiveFluxConfig.
 */
class NaiveFluxControl : public FluxControl {
public:

    LSST_CONTROL_FIELD(radius, double, "FIXME! NEVER DOCUMENTED!");

    NaiveFluxControl() : FluxControl("flux.naive"), radius(7.0) {}

private:
    virtual PTR(AlgorithmControl) _clone() const;
    virtual PTR(Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};

/**
 *  @brief C++ control object for PSF flux.
 *
 *  @sa PsfFluxConfig.
 */
class PsfFluxControl : public FluxControl {
public:

    PsfFluxControl() : FluxControl("flux.psf") {}

private:
    virtual PTR(AlgorithmControl) _clone() const;
    virtual PTR(Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};

/**
 *  @brief C++ control object for sinc aperture flux.
 *
 *  @sa SincFluxConfig.
 */
class SincFluxControl : public FluxControl {
public:

    LSST_CONTROL_FIELD(radius1, double, "major axis of inner boundary (pixels)");
    LSST_CONTROL_FIELD(radius2, double, "major axis of outer boundary (pixels)");
    LSST_CONTROL_FIELD(angle, double, "measured from x anti-clockwise; radians");
    LSST_CONTROL_FIELD(ellipticity, double, "1 - b/a");

    SincFluxControl() : 
        FluxControl("flux.sinc"), radius1(0.0), radius2(7.0), angle(0.0), ellipticity(0.0) {}

private:
    virtual PTR(AlgorithmControl) _clone() const;
    virtual PTR(Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};

}}}// namespace lsst::meas::algorithms

#endif // !LSST_MEAS_ALGORITHMS_FLUXCONTROL_H
