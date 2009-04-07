#if !defined(LSST_MEAS_ALGORITHMS_SDSSCENTROID_H)
#define LSST_MEAS_ALGORITHMS_SDSSCENTROID_H 1
/**
 * @file
 */
#include "lsst/meas/algorithms/Centroid.h"
#include "lsst/meas/algorithms/CentroidImpl.h"

namespace lsst { namespace meas { namespace algorithms {

/**
 * @brief A class that knows how to calculate centroids using the SDSS centroiding algorithm
 */
template<typename ImageT>
class SdssMeasureCentroid : public MeasureCentroid<ImageT> {
public:
    /**
     * @brief Return the (unique) instance of NaiveMeasureCentroid
     */
    static MeasureCentroid<ImageT>* getInstance() {
        if (_instance == NULL) {
            _instance = new SdssMeasureCentroid;
            MeasureCentroid<ImageT>::registerType("SDSS", SDSS);
        }
        return _instance;
    }
private:
    SdssMeasureCentroid() {}
    Centroid doApply(ImageT const& image, int x, int y, PSF const*, double background) const;

    static SdssMeasureCentroid* _instance;
};
}}}
#endif
