// -*- lsst-c++ -*-
#ifndef LSST_MEAS_ALGORITHMS_FOCALPLANE_H
#define LSST_MEAS_ALGORITHMS_FOCALPLANE_H

#include "lsst/meas/algorithms/Algorithm.h"

namespace lsst { namespace meas { namespace algorithms {

/**
 *  @brief C++ control object for focal plane position
 */
class FocalPlaneControl : public algorithms::AlgorithmControl {
public:
    FocalPlaneControl() : algorithms::AlgorithmControl("focalplane", 3.0) {}

private:
    virtual PTR(algorithms::AlgorithmControl) _clone() const {
        return boost::make_shared<FocalPlaneControl>(*this);
    }
    virtual PTR(algorithms::Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;

};

}}} // namespace lsst::meas::algorithms

#endif // !LSST_MEAS_ALGORITHMS_FOCALPLANE_H
