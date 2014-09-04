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
 
#ifndef LSST_MEAS_ALGORITHMS_CLASSIFICATION_H
#define LSST_MEAS_ALGORITHMS_CLASSIFICATION_H

#include "lsst/meas/algorithms/Algorithm.h"

namespace lsst { namespace meas { namespace algorithms {

/**
 *  @brief Control/factory for the algorithm that does star/galaxy classification
 *
 *  The algorithm class itself adds nothing to the public interface of its base class, so
 *  it is declared only in the source file.
 *
 *  This algorithm is based entirely on the ratio of Model to PSF fluxes; it never fails unless 
 *  one of these fails, so it does not have its own failure flag.
 *
 *  @todo this class needs a more specific name, especially now that classifiers are pluggable
 */
class ClassificationControl : public AlgorithmControl {
public:

    LSST_CONTROL_FIELD(fluxRatio, double, "critical ratio of model to psf flux");
    LSST_CONTROL_FIELD(modelErrFactor, double, "correction factor for modelFlux error");
    LSST_CONTROL_FIELD(psfErrFactor, double, "correction factor for psfFlux error");

    ClassificationControl() :
        AlgorithmControl("classification.extendedness", 5.0),
        fluxRatio(0.95), modelErrFactor(0.0), psfErrFactor(0.0)
    {}

    PTR(ClassificationControl) clone() const {
        return boost::static_pointer_cast<ClassificationControl>(_clone());
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
