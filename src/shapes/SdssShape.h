#if !defined(LSST_MEAS_ALGORITHMS_SDSSSHAPE_H)
#define LSST_MEAS_ALGORITHMS_SDSSSHAPE_H 1
/**
 * @file
 */
#include "lsst/meas/algorithms/Shape.h"
#include "lsst/meas/algorithms/ShapeImpl.h"

namespace lsst { namespace meas { namespace algorithms {

/**
 * @brief A class that knows how to calculate centroids using the SDSS centroiding algorithm
 */
template<typename ImageT>
class SdssmeasureShape : public measureShape<ImageT> {
public:
    /**
     * @brief Return the (unique) instance of SdssmeasureShape
     */
    static measureShape<ImageT>* getInstance() {
        if (_instance == NULL) {
            _instance = new SdssmeasureShape;
            measureShape<ImageT>::registerType("SDSS", SDSS);
        }
        return _instance;
    }
private:
    SdssmeasureShape() {}
    Shape doApply(ImageT const& image, int x, int y, PSF const*, double background) const;

    static SdssmeasureShape* _instance;
};
}}}
#endif
