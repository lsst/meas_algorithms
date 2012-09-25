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
 
#ifndef LSST_MEAS_ALGORITHMS_CENTROIDCONTROL_H
#define LSST_MEAS_ALGORITHMS_CENTROIDCONTROL_H
//!
// Control/algorithm hierarchy for centroid measurement.
//

#include "lsst/base.h"
#include "lsst/pex/config.h"
#include "lsst/meas/algorithms/Algorithm.h"

namespace lsst {
namespace meas {
namespace algorithms {

class CentroidControl;

/**
 *  @brief Intermediate base class for algorithms that compute a centroid.
 */
class CentroidAlgorithm : public Algorithm {
public:

    /**
     *  @brief Tuple type that holds the keys that define a standard centroid algorithm.
     *
     *  Algorithms are encouraged to add additional flags as appropriate, but these are required.
     */
    typedef afw::table::KeyTuple<afw::table::Centroid> KeyTuple;

    /// @copydoc Algorithm::getControl
    CentroidControl const & getControl() const;

    /// @brief Return the standard centroid keys registered by this algorithm.
    KeyTuple const & getKeys() const { return _keys; }

protected:

    /// @brief Initialize with a manually-constructed key tuple.
    CentroidAlgorithm(CentroidControl const & ctrl, KeyTuple const & keys);

    /// @brief Initialize using afw::table::addCentroid field to fill out repetitive descriptions.
    CentroidAlgorithm(CentroidControl const & ctrl, afw::table::Schema & schema, char const * doc);

private:
    KeyTuple _keys;
};

class CentroidControl : public AlgorithmControl {
public:

    PTR(CentroidControl) clone() const { return boost::static_pointer_cast<CentroidControl>(_clone()); }

    PTR(CentroidAlgorithm) makeAlgorithm(
        afw::table::Schema & schema,
        PTR(daf::base::PropertyList) const & metadata = PTR(daf::base::PropertyList)(),
        AlgorithmControlMap const & others = AlgorithmControlMap(),
        bool isForced = false
    ) const {
        return boost::static_pointer_cast<CentroidAlgorithm>(
            _makeAlgorithm(schema, metadata, others, isForced)
        );
    }

protected:
    explicit CentroidControl(std::string const & name_, double priority=0.0) :
        AlgorithmControl(name_, priority) {}
};

inline CentroidAlgorithm::CentroidAlgorithm(CentroidControl const & ctrl, KeyTuple const & keys) :
    Algorithm(ctrl), _keys(keys)
{}

inline CentroidAlgorithm::CentroidAlgorithm(
    CentroidControl const & ctrl, afw::table::Schema & schema, char const * doc
) :
    Algorithm(ctrl), _keys(afw::table::addCentroidFields(schema, ctrl.name, doc))
{}

inline CentroidControl const & CentroidAlgorithm::getControl() const {
    return static_cast<CentroidControl const &>(Algorithm::getControl());
}

/**
 *  @brief C++ control object for Gaussian centroid.
 */
class GaussianCentroidControl : public CentroidControl {
public:

    GaussianCentroidControl() : CentroidControl("centroid.gaussian") {}

private:
    virtual PTR(AlgorithmControl) _clone() const;
    virtual PTR(Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};

/**
 *  @brief C++ control object for naive centroid.
 */
class NaiveCentroidControl : public CentroidControl {
public:

    LSST_CONTROL_FIELD(background, double, "FIXME! NEVER DOCUMENTED!");

    NaiveCentroidControl() : CentroidControl("centroid.naive"), background(0.0) {}

private:
    virtual PTR(AlgorithmControl) _clone() const;
    virtual PTR(Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};

/**
 *  @brief C++ control object for SDSS centroid.
 */
class SdssCentroidControl : public CentroidControl {
public:

    LSST_CONTROL_FIELD(binmax, int, "maximum allowed binning");
    LSST_CONTROL_FIELD(peakMin, double, "if the peak's less thatn this insist on binning at least once");
    LSST_CONTROL_FIELD(wfac, double, "fiddle factor for adjusting the binning");

    SdssCentroidControl() : CentroidControl("centroid.sdss"), binmax(16), peakMin(-1.0), wfac(1.5) {}

private:
    virtual PTR(AlgorithmControl) _clone() const;
    virtual PTR(Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};

}}}// namespace lsst::meas::algorithms

#endif // !LSST_MEAS_ALGORITHMS_CENTROIDCONTROL_H
