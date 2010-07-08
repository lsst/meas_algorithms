// -*- LSST-C++ -*-
#if !defined(LSST_DETECTION_CR_H)
#define LSST_DETECTION_CR_H
//!
// Handle cosmic rays in a MaskedImage
//
#include <vector>
#include "lsst/base.h"
#include "lsst/afw/image/MaskedImage.h"

namespace lsst { namespace afw {
namespace detection {
    class Footprint;
    class Psf;
}}}

namespace lsst {
namespace meas {
namespace algorithms {

template <typename MaskedImageT>
std::vector<PTR(lsst::afw::detection::Footprint)>
findCosmicRays(MaskedImageT& image,
               lsst::afw::detection::Psf const &psf,
               double const bkgd,
               lsst::pex::policy::Policy const& policy,
               bool const keep = false
              );

}}}

#endif
