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
 
#ifndef LSST_MEAS_ALGORITHMS_SKYCOORD_H
#define LSST_MEAS_ALGORITHMS_SKYCOORD_H

#include "lsst/meas/algorithms/Algorithm.h"

namespace lsst { namespace meas { namespace algorithms {

/**
 *  @brief Control/factory for an algorithm that sets ra/dec from x/y.
 *
 *  The algorithm class itself adds nothing to the public interface of its base class, so
 *  it is declared only in the source file.
 *
 *  This is almost certainly not the right way to set ra/dec for sources in production,
 *  but it's very useful during development.
 *
 *  Unlike most algorithms, SkyCoord does "own" any fields in the Schema; it sets the ra/dec
 *  fields that are included in the minimal source schema, and uses the centroid slot as an
 *  input (note that this means the slots must be defined on the table before running the
 *  algorithm).
 */
class SkyCoordControl : public AlgorithmControl {
public:

    SkyCoordControl() : AlgorithmControl("skycoord", 5.0) {}

    PTR(SkyCoordControl) clone() const {
        return boost::static_pointer_cast<SkyCoordControl>(_clone());
    }
    
private:

    virtual PTR(AlgorithmControl) _clone() const;

    virtual PTR(Algorithm) _makeAlgorithm(
        afw::table::Schema & schema,
        PTR(daf::base::PropertyList) const & metadata
    ) const;
    
};

}}} // namespace lsst::meas::algorithms

#endif
