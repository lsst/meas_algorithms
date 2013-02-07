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
 
#ifndef LSST_MEAS_ALGORITHMS_SDSSSHAPECONTROL_H
#define LSST_MEAS_ALGORITHMS_SDSSSHAPECONTROL_H
//!
// SDSS shape measurement interface
//
#include "lsst/meas/algorithms/ShapeControl.h"
#include "lsst/meas/algorithms/detail/SdssShape.h"

namespace lsst {
namespace meas {
namespace algorithms {

/**
 *  @brief C++ control object for SDSS shape.
 *
 *  @sa SdssShapeConfig.
 */
class SdssShapeControl : public ShapeControl {
public:

    LSST_CONTROL_FIELD(background, double, "Additional value to add to background");
    LSST_CONTROL_FIELD(maxIter, int, "Maximum number of iterations");
    LSST_CONTROL_FIELD(tol1, float, "Convergence tolerance for e1,e2");
    LSST_CONTROL_FIELD(tol2, float, "Convergence tolerance for FWHM");

    SdssShapeControl() : ShapeControl("shape.sdss"), background(0.0), maxIter(detail::sdssShapeMaxIter),
                         tol1(detail::sdssShapeTol1), tol2(detail::sdssShapeTol2) {}

private:
    virtual PTR(AlgorithmControl) _clone() const;
    virtual PTR(Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};

}}}// namespace lsst::meas::algorithms

#endif // !LSST_MEAS_ALGORITHMS_SDSSSHAPECONTROL_H
