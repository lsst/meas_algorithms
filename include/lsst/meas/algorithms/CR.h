#if !defined(LSST_DETECTION_CR_H)
#define LSST_DETECTION_CR_H
//!
// Handle cosmic rays in a MaskedImage
//
#include <vector>
#include <lsst/afw/image/MaskedImage.h>
#include "lsst/afw/detection/Footprint.h"
#include "lsst/detection/PSF.h"

class lsst::pex::policy::Policy;

namespace lsst { namespace detection {

template <typename MaskedImageT>
std::vector<lsst::afw::detection::Footprint::Ptr>
findCosmicRays(MaskedImageT& image,
               lsst::detection::PSF const &psf,
               float const bkgd,
               lsst::pex::policy::Policy const& policy,
               bool const keep = false
              );

}}

#endif
