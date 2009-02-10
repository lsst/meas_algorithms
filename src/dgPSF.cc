/*!
 * Represent a PSF as a circularly symmetrical double Gaussian
 *
 * \file
 *
 * \ingroup algorithms
 */
#include <cmath>
#include "lsst/meas/algorithms/PSF.h"

namespace lsst { namespace meas { namespace algorithms {

/************************************************************************************************************/
/**
 * Constructor for a dgPSF
 */
dgPSF::dgPSF(double sigma1,             ///< Width of inner Gaussian
             double sigma2,             ///< Width of outer Gaussian
             double b,                  ///< Central amplitude of outer Gaussian (inner amplitude == 1)
             int size                   ///< Kernel should have dimensions (size*size)
            ) :
    PSF(),
    _sigma1(sigma1),
    _sigma2(sigma2),
    _b(b)
{
    if (size > 0) {
        lsst::afw::math::DoubleGaussianFunction2<double> dg(sigma1, sigma2, b);
        setKernel(lsst::afw::math::Kernel::PtrT(new lsst::afw::math::AnalyticKernel(size, size, dg)));
    }
}

/// \brief Evaluate the PSF at (dx, dy) (relative to the centre), taking the central amplitude to be 1.0
double dgPSF::getValue(double const dx,            ///< Desired column (relative to centre)
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

}}}
