/// \file

#include "lsst/pex/logging/Trace.h"
#include "lsst/meas/algorithms/Measure.h"

namespace image = lsst::afw::image;
namespace detection = lsst::afw::detection;
namespace algorithms = lsst::meas::algorithms;

/************************************************************************************************************/
/**
 * \brief Calculate a detected source's moments
 */
template <typename MaskedImageT>
class FootprintCentroid : public detection::FootprintFunctor<MaskedImageT> {
public:
    FootprintCentroid(MaskedImageT const& mimage                    ///< The image the source lives in
                     ) : detection::FootprintFunctor<MaskedImageT>(mimage), _n(0), _sum(0), _sumx(0), _sumy(0) {}

    /// \brief method called for each pixel by apply()
    void operator()(typename MaskedImageT::xy_locator loc, ///< locator pointing at the pixel
                    int x,                                 ///< column-position of pixel
                    int y                                  ///< row-position of pixel
                   ) {
        typename MaskedImageT::Image::Pixel val = loc.image(0, 0);

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

/************************************************************************************************************/
/**
 * \brief Set some fields in a Source from foot (which was found in mimage)
 */
template<typename MaskedImageT>
void algorithms::measureSource(lsst::afw::detection::Source::Ptr src,    ///< the Source to receive results
                               MaskedImageT& mimage,      ///< image wherein Footprint dwells
                               lsst::afw::detection::Footprint const& foot, ///< Footprint to measure
                               float background                ///< background level to subtract
                              ) {
    FootprintCentroid<MaskedImageT> centroid(mimage);
    centroid.apply(foot);
    
    src->setColc(centroid.getX());
    src->setRowc(centroid.getY());
    src->setFlux(centroid.getSum());
}

//
// Explicit instantiations
//
template void algorithms::measureSource(detection::Source::Ptr src, image::MaskedImage<float>& mimage,
                                        detection::Footprint const &foot, float background);
//
// Why do we need double images?
//
#if 1
template void algorithms::measureSource(detection::Source::Ptr src, image::MaskedImage<double>& mimage,
                                        detection::Footprint const &foot, float background);
#endif
