#if !defined(LSST_MEAS_ALGORITHMS_NAIVEPHOTOMETRY_H)
#define LSST_MEAS_ALGORITHMS_NAIVEPHOTOMETRY_H 1
/**
 * @file NaivePhotometry.h
 * @brief Compute Aperture and PSF fluxes in a simple way
 * @ingroup meas/algorithms
 * @author Steve Bickerton (adapted from RHL's Shape class)
 */
#include "lsst/meas/algorithms/Photometry.h"
#include "lsst/meas/algorithms/PhotometryImpl.h"

namespace lsst { namespace meas { namespace algorithms {

/**
 * @class A class that knows how to calculate fluxes using the NAIVE photometry algorithm
 * @ingroup meas/algorithms
 */
template<typename MaskedImageT>
class measureNaivePhotometry : public measurePhotometry<MaskedImageT> {
public:
    /**
     * @brief Return the (unique) instance of measureNaivePhotometry
     */
    static measurePhotometry<MaskedImageT>* getInstance(float const radius) {
        if (_instance == NULL) {
            _instance = new measureNaivePhotometry(radius);
            measurePhotometry<MaskedImageT>::registerType("NAIVE", NAIVE);
        }
        _instance->setRadius(radius);
        return _instance;
    }
private:
    measureNaivePhotometry(float const radius) :  measurePhotometry<MaskedImageT>(radius) {}
    Photometry doApply(MaskedImageT const& image, double xcen, double ycen, PSF const*, double background) const;

    static measureNaivePhotometry* _instance;
};
}}}
#endif
