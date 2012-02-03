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
// Control (and secretly, factory) object hierarchy for shape algorithms.
//

#include "lsst/base.h"
#include "lsst/pex/config.h"
#include "lsst/meas/algorithms/Algorithm.h"

namespace lsst {
namespace meas {
namespace algorithms {

typedef AlgorithmControl<afw::detection::Shape> ShapeControl;

/**
 *  @brief C++ control object for SDSS shape.
 *
 *  @sa SdssShapeConfig.
 */
class SdssShapeControl : public ShapeControl {
public:

    LSST_CONTROL_FIELD(background, double, "FIXME! NEVER DOCUMENTED!");

    SdssShapeControl() : background(0.0) {}

private:
    LSST_ALGORITHM_CONTROL_PRIVATE_DECL()
};

}}}// namespace lsst::meas::algorithms

#endif // !LSST_MEAS_ALGORITHMS_SHAPECONTROL_H
