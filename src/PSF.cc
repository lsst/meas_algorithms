/*!
 * \brief Implementation of PSF code
 *
 * \file
 *
 * \ingroup algorithms
 */
#include <typeinfo>
#include <cmath>
#include "lsst/meas/algorithms/PSF.h"

namespace algorithms = lsst::meas::algorithms;

algorithms::PSF::PSF(lsst::afw::math::Kernel::PtrT kernel ///< The Kernel corresponding to this PSF
                    ) : lsst::daf::data::LsstBase(typeid(this)),
                        _kernel(kernel),
                        _A(1)
{
    ;
}

/// PSF's destructor; declared pure virtual, but we still need an implementation
algorithms::PSF::~PSF() {}

///
/// Set the PSF's kernel
///
void algorithms::PSF::setKernel(lsst::afw::math::Kernel::PtrT kernel) {
    _kernel = kernel;
}

///
/// Return the PSF's kernel
///
lsst::afw::math::Kernel::PtrT algorithms::PSF::getKernel() {
    return _kernel;
}

///
/// Return the PSF's kernel
///
boost::shared_ptr<const lsst::afw::math::Kernel> algorithms::PSF::getKernel() const {
    return boost::shared_ptr<const lsst::afw::math::Kernel>(_kernel);
}
