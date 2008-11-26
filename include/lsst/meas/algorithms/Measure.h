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

        for (lsst::afw::detection::Footprint::SpanList::const_iterator siter = _foot.getSpans().begin();
             siter != _foot.getSpans().end(); siter++) {
            lsst::afw::detection::Span::Ptr const span = *siter;

            int const y = span->getY();
            int const x0 = span->getX0();
            int const x1 = span->getX1();

            typename MaskedImageT::xy_locator loc = _mimage.xy_at(x0, y);
            for (int x = x0; x <= x1; ++x, ++loc.x()) {
                operator()(loc, x, y);
            }
        }
    }

    virtual void operator()(typename MaskedImageT::xy_locator loc, int x, int y) = 0;
private:
    lsst::afw::detection::Footprint const& _foot;
    MaskedImageT const& _mimage;
};

template <typename MaskedImageT>
class FootprintCentroid : public FootprintFunctor<MaskedImageT> {
public:
    FootprintCentroid(lsst::afw::detection::Footprint const& foot,
                     MaskedImageT const& mimage
                     ) : FootprintFunctor<MaskedImageT>(foot, mimage), _n(0), _sum(0), _sumx(0), _sumy(0) {}

    void operator()(typename MaskedImageT::xy_locator loc, int x, int y) {
        typename MaskedImageT::Image::Pixel val = loc.image(0, 0);

        _n++;
        _sum += val;
        _sumx += lsst::afw::image::indexToPosition(x)*val;
        _sumy += lsst::afw::image::indexToPosition(y)*val;
    }

    int getN() const { return _n; }
    double getSum() const { return _sum; }
    double getX() const { return _sumx/_sum; }
    double getY() const { return _sumy/_sum; }
private:
    int _n;
    double _sum, _sumx, _sumy;
};

}}}
#endif // LSST_DETECTION_MEASURE_H
