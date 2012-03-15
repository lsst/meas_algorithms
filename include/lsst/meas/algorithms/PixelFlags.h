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
 
#if !defined(LSST_MEAS_ALGORITHMS_PIXELFLAGS_H)
#define LSST_MEAS_ALGORITHMS_PIXELFLAGS_H

#include "lsst/meas/algorithms/Algorithm.h"

namespace lsst { namespace meas { namespace algorithms {

/**
 *  @brief Control/factory for the algorithm that sets source bits from pixel mask bits.
 *
 *  The algorithm class itself adds nothing to the public interface of its base class, so
 *  it is declared only in the source file.
 */
class PixelFlagControl : public AlgorithmControl {
public:

    PixelFlagControl() : AlgorithmControl("flags.pixel", 0.0) {}

    PTR(PixelFlagControl) clone() const { return boost::static_pointer_cast<PixelFlagControl>(_clone()); }
    
private:

    virtual PTR(AlgorithmControl) _clone() const;

    virtual PTR(Algorithm) _makeAlgorithm(
        afw::table::Schema & schema,
        PTR(daf::base::PropertyList) const & metadata
    ) const;
    
};


}}} // namespace

#endif
