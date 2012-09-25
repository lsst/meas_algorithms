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
 
#ifndef LSST_MEAS_ALGORITHMS_SHAPECONTROL_H
#define LSST_MEAS_ALGORITHMS_SHAPECONTROL_H
//!
// Control/algorithm hierarchy for shape measurements.
//

#include "lsst/base.h"
#include "lsst/pex/config.h"
#include "lsst/meas/algorithms/Algorithm.h"

namespace lsst {
namespace meas {
namespace algorithms {

class ShapeControl;

/**
 *  @brief Intermediate base class for algorithms that compute a shape.
 */
class ShapeAlgorithm : public Algorithm {
public:

    /**
     *  @brief Tuple type that holds the keys that define a standard shape algorithm.
     *
     *  Algorithms are encouraged to add additional flags as appropriate, but these are required.
     */
    typedef afw::table::KeyTuple<afw::table::Shape> KeyTuple;

    /// @copydoc Algorithm::getControl
    ShapeControl const & getControl() const;

    /// @brief Return the standard shape keys registered by this algorithm.
    KeyTuple const & getKeys() const { return _keys; }

protected:

    /// @brief Initialize with a manually-constructed key tuple.
    ShapeAlgorithm(ShapeControl const & ctrl, KeyTuple const & keys);

    /// @brief Initialize using afw::table::addShape field to fill out repetitive descriptions.
    ShapeAlgorithm(ShapeControl const & ctrl, afw::table::Schema & schema, char const * doc);

private:
    KeyTuple _keys;
};

class ShapeControl : public AlgorithmControl {
public:

    PTR(ShapeControl) clone() const { return boost::static_pointer_cast<ShapeControl>(_clone()); }

    PTR(ShapeAlgorithm) makeAlgorithm(
        afw::table::Schema & schema,
        PTR(daf::base::PropertyList) const & metadata = PTR(daf::base::PropertyList)(),
        AlgorithmControlMap const & others = AlgorithmControlMap(),
        bool isForced = false
    ) const {
        return boost::static_pointer_cast<ShapeAlgorithm>(_makeAlgorithm(schema, metadata, others, isForced));
    }

protected:
    explicit ShapeControl(std::string const & name_) : AlgorithmControl(name_, 1.0) {}
};

inline ShapeAlgorithm::ShapeAlgorithm(ShapeControl const & ctrl, KeyTuple const & keys) :
    Algorithm(ctrl), _keys(keys)
{}

inline ShapeAlgorithm::ShapeAlgorithm(
    ShapeControl const & ctrl, afw::table::Schema & schema, char const * doc
) :
    Algorithm(ctrl), _keys(afw::table::addShapeFields(schema, ctrl.name, doc))
{}

inline ShapeControl const & ShapeAlgorithm::getControl() const {
    return static_cast<ShapeControl const &>(Algorithm::getControl());
}

/**
 *  @brief C++ control object for SDSS shape.
 *
 *  @sa SdssShapeConfig.
 */
class SdssShapeControl : public ShapeControl {
public:

    LSST_CONTROL_FIELD(background, double, "FIXME! NEVER DOCUMENTED!");

    SdssShapeControl() : ShapeControl("shape.sdss"), background(0.0) {}

private:
    virtual PTR(AlgorithmControl) _clone() const;
    virtual PTR(Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};

}}}// namespace lsst::meas::algorithms

#endif // !LSST_MEAS_ALGORITHMS_SHAPECONTROL_H
