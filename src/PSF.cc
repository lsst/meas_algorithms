/*!
 * \brief Implementation of PSF code
 *
 * \file
 *
 * \ingroup detection
 */
#include <typeinfo>
#include <cmath>
#include "lsst/detection/PSF.h"

using namespace lsst::detection;

lsst::detection::PSF::PSF() :
    lsst::daf::data::LsstBase(typeid(this)),
    _A(1)
{
    ;
}

/*
 * Represent a PSF as a circularly symmetrical double Gaussian
 */
lsst::detection::dgPSF::dgPSF(
                              double sigma1, ///< Width of inner Gaussian
                              double sigma2, ///< Width of outer Gaussian
                              double b  ///< Central amplitude of outer Gaussian (inner amplitude == 1)
                      ) :
    PSF(),
    _sigma1(sigma1),
    _sigma2(sigma2),
    _b(b)
{
    ;
}

std::string lsst::detection::dgPSF::toString() const {
    return "Hello world";
}

// \brief Evaluate the PSF at (col, row), taking the central amplitude to be 1.0
double lsst::detection::dgPSF::getValue(double const dcol, ///< Desired column (relative to centre)
                                        double const drow ///< Desired row
                                       ) const {
    double const r2 = dcol*dcol + drow*drow;
    double const psf1 = exp(-r2/(2*_sigma1*_sigma1));
    if (_b == 0) {
        return psf1;
    }
    
    double const psf2 = exp(-r2/(2*_sigma2*_sigma2));

    return (psf1 + _b*psf2)/(1 + _b);
}
