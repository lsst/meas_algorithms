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
template<typename MaskedImageT>
class SdssMeasureShape : public MeasureShape<MaskedImageT> {
public:
    SdssMeasureShape(int type=0) : MeasureShape<MaskedImageT>() {
        static bool _registered = false;

        if (!_registered) {
            SdssMeasureShape::declare("SDSS", new MeasureShapeFactory<SdssMeasureShape>());
            _registered = true;
        }        
    }
private:
    Shape doApply(MaskedImageT const& image, double xcen, double ycen, PSF const*, double background) const;
};
}}}
#endif
