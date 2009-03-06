/*!
 * Represent a PSF as a circularly symmetrical double Gaussian
 *
 * \file
 *
 * \ingroup algorithms
 */
#include <cmath>
#include "lsst/pex/exceptions.h"
#include "lsst/meas/algorithms/detail/dgPSF.h"
#include "lsst/meas/algorithms/detail/PSFImpl.h"

namespace lsst { namespace meas { namespace algorithms {

/************************************************************************************************************/
/**
 * Constructor for a dgPSF
 */
dgPSF::dgPSF(int size,                  ///< Kernel should have dimensions (size*size)
             double sigma1,             ///< Width of inner Gaussian
             double sigma2,             ///< Width of outer Gaussian
             double b                   ///< Central amplitude of outer Gaussian (inner amplitude == 1)
            ) :
    PSF(),
    _sigma1(sigma1), _sigma2(sigma2), _b(b)
{
    static bool first = true;
    if (first) {
        dgPSF::registerType("DGPSF", DGPSF);
        first = false;
    }

    if (b == 0.0 && sigma2 == 0.0) {
        _sigma2 = 1.0;                  // avoid 0/0 at centre of PSF
    }

    if (_sigma1 == 0 || _sigma2 == 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::DomainErrorException,
                          (boost::format("sigma may not be 0: %g, %g") % _sigma1 % _sigma2).str());
    }
    
    if (size > 0) {
        lsst::afw::math::DoubleGaussianFunction2<double> dg(_sigma1, _sigma2, _b);
        setKernel(lsst::afw::math::Kernel::PtrT(new lsst::afw::math::AnalyticKernel(size, size, dg)));
    }
}

/// \brief Evaluate the PSF at (dx, dy) (relative to the centre), taking the central amplitude to be 1.0
double dgPSF::doGetValue(double const dx,            ///< Desired column (relative to centre)
                         double const dy             ///< Desired row
                        ) const {
    double const r2 = dx*dx + dy*dy;
    double const psf1 = exp(-r2/(2*_sigma1*_sigma1));
    if (_b == 0) {
        return psf1;
    }
    
    double const psf2 = exp(-r2/(2*_sigma2*_sigma2));

    return (psf1 + _b*psf2)/(1 + _b);
}

//
// We need to make an instance here so as to register it with Centroider
//
// \cond
namespace {                                                 \
    PSF* foo = new dgPSF(0, 1.0);
}

// \endcond
}}}
