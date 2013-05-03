// -*- lsst-c++ -*-
#ifndef LSST_MEAS_ALGORITHMS_JACOBIAN_H
#define LSST_MEAS_ALGORITHMS_JACOBIAN_H

#include "lsst/meas/algorithms/Algorithm.h"

namespace lsst { namespace meas { namespace algorithms {

/**
 *  @brief C++ control object for Jacobian correction
 *
 *  @sa JacobianConfig.
 */
class JacobianControl : public algorithms::AlgorithmControl {
public:
    JacobianControl() : algorithms::AlgorithmControl("jacobian", 3.0), pixelScale(0.5) {}
    LSST_CONTROL_FIELD(pixelScale, float, "Nominal pixel size (arcsec)");

private:
    virtual PTR(algorithms::AlgorithmControl) _clone() const {
        return boost::make_shared<JacobianControl>(*this);
    }
    virtual PTR(algorithms::Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;

};

}}} // namespace lsst::meas::algorithms

#endif // !LSST_MEAS_ALGORITHMS_JACOBIAN_H
