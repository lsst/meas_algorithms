#if !defined(LSST_DETECTION_MEASURE_H)
#define LSST_DETECTION_MEASURE_H
//!
// Measure properties of an image selected by a Footprint
//
#include <list>
#include <cmath>
#include "boost/cstdint.hpp"
#include "boost/shared_ptr.hpp"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/image/ImageUtils.h"
#include "lsst/afw/detection.h"
#include "lsst/meas/algorithms/PSF.h"

namespace lsst {
namespace meas {
namespace algorithms {

struct Flags {
    enum {
        EDGE                      = 0x1,    // object is in region labelled EDGE
        SHAPE_SHIFT               = 0x2,    // centroid shifted while estimating adaptive moments
        SHAPE_MAXITER             = 0x4,    // too many iterations for adaptive moments
        SHAPE_UNWEIGHTED          = 0x8,    // "adaptive" moments are unweighted
        SHAPE_UNWEIGHTED_PSF      = 0x10,   // the PSF's "adaptive" moments are unweighted
        SHAPE_UNWEIGHTED_BAD      = 0x20,   // even the unweighted moments were bad
        PEAKCENTER                = 0x40,   // given centre is position of peak pixel
        BINNED1                   = 0x80,   // object was found in 1x1 binned image
    };
};
            
/************************************************************************************************************/
//
// Actually measure an object
//
template<typename MaskedImageT>
void measureSource(lsst::afw::detection::Source::Ptr src, MaskedImageT& mimage,
                   lsst::afw::detection::Footprint const& fp,
                   lsst::pex::policy::Policy const& policy,
                   float background, PSF const* psf=NULL);
}}}
#endif // LSST_DETECTION_MEASURE_H
