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

namespace lsst { namespace meas { namespace algorithms {

/************************************************************************************************************/
//
// Actually measure an object
//
template<typename MaskedImageT>
void measureSource(lsst::afw::detection::Source::Ptr src, MaskedImageT& mimage,
                   lsst::afw::detection::Footprint const& fp, float background);
            
}}}
#endif // LSST_DETECTION_MEASURE_H
