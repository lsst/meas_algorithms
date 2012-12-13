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
 

/**
 * @file
 */

#include "lsst/meas/algorithms/Measure.h"
#include "lsst/meas/algorithms/CentroidControl.h"

namespace lsst {
namespace meas {
namespace algorithms {

class SillyCentroidControl : public CentroidControl {
public:
    SillyCentroidControl(int dY_=1) : CentroidControl("centroid.silly"), dY(dY_) {}
    LSST_CONTROL_FIELD(dY, int, "Number of pixels to offset the centroid in y");
private:
    virtual PTR(AlgorithmControl) _clone() const { return boost::make_shared<SillyCentroidControl>(*this); }

    virtual PTR(Algorithm) _makeAlgorithm(
        afw::table::Schema & schema,         
        PTR(daf::base::PropertyList) const & metadata,
        AlgorithmControlMap const & other) const;
};

}}}
