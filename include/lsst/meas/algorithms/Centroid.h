#ifndef LSST_MEAS_ALGORITHMS_CENTROID_H
#define LSST_MEAS_ALGORITHMS_CENTROID_H
//!
// Measure centroid of a Footprint
//


#include "lsst/afw/detection.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/image/ImageUtils.h"

namespace lsst {
namespace meas {
namespace algorithms {



/******************************************************************************/
/**
 * \brief Calculate the footprint's flux weighted centroid
 */
template <typename MaskedImageT>
class FootprintCentroid : public lsst::afw::detection::FootprintFunctor<typename MaskedImageT> {
public:
    FootprintCentroid(MaskedImageT const& mimage)          ///< The image the source lives in
        : detection::FootprintFunctor<MaskedImageT>(mimage), 
          _n(0), _sum(0), _sumx(0), _sumy(0) {
    }

    /// \brief method called for each pixel by apply()
    void operator()(typename MaskedImageT::xy_locator loc, ///< locator pointing at the pixel
                    int x,                                 ///< column-position of pixel
                    int y                                  ///< row-position of pixel
    ) {
        typename MaskedImageT::Image::Pixel val = loc.image();

        _n++;
        _sum += val;
        _sumx += lsst::afw::image::indexToPosition(x)*val;
        _sumy += lsst::afw::image::indexToPosition(y)*val;
    }

    /// Return the number of pixels
    int getN() const { return _n; }
    /// Return the Footprint's flux
    double getSum() const { return _sum; }
    /// Return the Footprint's column centroid
    double getX() const { return _sumx/_sum; }
    /// Return the Footprint's row centroid
    double getY() const { return _sumy/_sum; }
private:
    int _n;
    double _sum, _sumx, _sumy;
};


}}} //namespace lsst::meas::algorithms

 
 #endif //LSST_MEAS_ALGORITHMS_CENTROID_H
