#if !defined(LSST_MEAS_ALGORITHMS_NAIVECENTROID_H)
#define LSST_MEAS_ALGORITHMS_NAIVECENTROID_H 1
/**
 * @file
 */
#include "lsst/meas/algorithms/Centroid.h"
#include "lsst/meas/algorithms/CentroidImpl.h"

namespace lsst { namespace meas { namespace algorithms {

/**
 * @brief A class that knows how to calculate centroids as a simple unweighted first moment
 * of the 3x3 region around a pixel
 */
template<typename ImageT>
class NaiveMeasureCentroid : public MeasureCentroid<ImageT> {
public:
    /**
     * @brief Return the (unique) instance of NaiveMeasureCentroid
     */
    static MeasureCentroid<ImageT>* getInstance() {
        if (_instance == NULL) {
            _instance = new NaiveMeasureCentroid;
            MeasureCentroid<ImageT>::registerType("NAIVE", NAIVE);
        }
        return _instance;
    }
private:
    NaiveMeasureCentroid() {}
    Centroid doApply(ImageT const& image, int x, int y, PSF const* psf, double background) const;

    static NaiveMeasureCentroid* _instance;
};
}}}
#endif
