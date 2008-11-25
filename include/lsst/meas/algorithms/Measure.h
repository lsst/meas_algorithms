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
#if 0
#include "lsst/afw/detection.h"         // requires afw > 3.2
#else
#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/detection/Source.h"
#endif

namespace lsst { namespace detection {

/*!
 * \brief Measure properties of an image selected by a Footprint
 *
 */
template<typename MaskedImageT>
class Measure : public lsst::daf::data::LsstBase {
public:
    Measure(MaskedImageT& img);
    void measureSource(lsst::afw::detection::Source::Ptr, lsst::afw::detection::Footprint const &fp, float background=0);
    void measureSource(lsst::afw::detection::Source::Ptr, lsst::afw::detection::Footprint::Ptr fpPtr, float background=0);
private:
    MaskedImageT _img;
};

}}
#endif // LSST_DETECTION_MEASURE_H
