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

template <typename MaskedImageT>
class FootprintFunctor {
public:
    FootprintFunctor(lsst::afw::detection::Footprint const& foot,
                     MaskedImageT const& mimage
                    ) : _foot(foot), _mimage(mimage) {}

    virtual ~FootprintFunctor() {}

    void apply() {
        if (_foot.getSpans().empty()) {
            return;
        }

        int ox1 = 0, oy = 0;            // Current position of the locator (in the SpanList loop)
        typename MaskedImageT::xy_locator loc = _mimage.xy_at(ox1, oy);

        for (lsst::afw::detection::Footprint::SpanList::const_iterator siter = _foot.getSpans().begin();
             siter != _foot.getSpans().end(); siter++) {
            lsst::afw::detection::Span::Ptr const span = *siter;

            int const y = span->getY();
            int const x0 = span->getX0();
            int const x1 = span->getX1();

            loc += lsst::afw::image::pair2I(x0 - ox1, y - oy);

            for (int x = x0; x <= x1; ++x, ++loc.x()) {
                operator()(loc, x, y);
            }

            ox1 = x1 + 1; oy = y;
        }
    }

    MaskedImageT const& getImage() const { return _mimage; }    

    virtual void operator()(typename MaskedImageT::xy_locator loc, int x, int y) = 0;
private:
    lsst::afw::detection::Footprint const& _foot;
    MaskedImageT const& _mimage;
};

/************************************************************************************************************/
//
// Actually measure an object
//
template<typename MaskedImageT>
void measureSource(lsst::afw::detection::Source::Ptr src, MaskedImageT& mimage,
                   lsst::afw::detection::Footprint const& fp, float background);
            
}}}
#endif // LSST_DETECTION_MEASURE_H
