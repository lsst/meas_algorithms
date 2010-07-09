// -*- LSST-C++ -*-
/*!
 * Represent a PSF as a circularly symmetrical single Gaussian
 *
 * \file
 *
 * \ingroup algorithms
 */
#include <cmath>
#include "lsst/pex/exceptions.h"
#include "lsst/meas/algorithms/detail/sgPsf.h"
#include "lsst/afw/image/ImageUtils.h"

namespace afwMath = lsst::afw::math;

namespace lsst {
namespace meas {
namespace algorithms {

/************************************************************************************************************/
/**
 * Constructor for a sgPsf
 */
sgPsf::sgPsf(int width,                         ///< Number of columns in realisations of PSF
             int height,                        ///< Number of rows in realisations of PSF
             double sigma,                       ///< Width of Gaussian
             double,        ///< needed to match PSF.h definition but unused
             double         ///< needed to match PSF.h definition but unused
            ) :
    KernelPsf(), _sigma(sigma)
{
    if (sigma <= 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::DomainErrorException,
                          (boost::format("sigma may not be 0: %g") % sigma).str());
    }
    
    if (width > 0) {
        afwMath::GaussianFunction1<double> sg(sigma);
        setKernel(afwMath::SeparableKernel::Ptr(new afwMath::SeparableKernel(width, height, sg, sg)));
    }
}

//
// We need to make an instance here so as to register it
//
// \cond
namespace {
    volatile bool isInstance =
        lsst::afw::detection::Psf::registerMe<sgPsf,
                                              boost::tuple<int, int, double,double,double> >("SingleGaussian");
}

// \endcond
}}}
