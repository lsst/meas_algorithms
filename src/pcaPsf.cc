// -*- LSST-C++ -*-
/*!
 * \brief Represent a PSF as a linear combination of PCA (== Karhunen-Loeve) basis functions
 *
 * \file
 *
 * \ingroup algorithms
 */
#include <cmath>
#include <numeric>
#include "lsst/base.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/image/ImageUtils.h"
#include "lsst/afw/math/Statistics.h"
#include "lsst/meas/algorithms/detail/pcaPsf.h"

namespace afwDetection = lsst::afw::detection;
namespace afwImage = lsst::afw::image;
namespace lsst {
namespace meas {
namespace algorithms {

/************************************************************************************************************/
/**
 * Constructor for a pcaPsf
 */
pcaPsf::pcaPsf(PTR(lsst::afw::math::Kernel) kernel ///< The desired Kernel
              ) : afwDetection::KernelPsf(kernel)
{
    //
    // Check that it's a LinearCombinationKernel
    //
    if (kernel.get() != NULL &&
        dynamic_cast<lsst::afw::math::LinearCombinationKernel *>(kernel.get()) == NULL) {
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                          "pcaPsf expects a LinearCombinationKernel");
    }
}

//
// We need to make an instance here so as to register it with createPSF
//
// \cond
namespace {
    volatile bool isInstance =
        lsst::afw::detection::Psf::registerMe<pcaPsf, PTR(lsst::afw::math::Kernel)>("PCA");
}

// \endcond
}}}
