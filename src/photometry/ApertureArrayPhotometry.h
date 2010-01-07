#if !defined(LSST_MEAS_ALGORITHMS_APERTUREARRAYPHOTOMETRY_H)
#define LSST_MEAS_ALGORITHMS_APERTUREARRAYPHOTOMETRY_H 1
/**
 * @file ApertureArrayPhotometry.h
 * @brief Compute Aperture fluxes for an array of apertures
 * @ingroup meas/algorithms
 * @author Tim Axelrod (adapted from Steve Bickerton's NaivePhotometry class)
 */
#include "lsst/meas/algorithms/Photometry.h"
#include "lsst/meas/algorithms/PhotometryImpl.h"

namespace lsst { namespace meas { namespace algorithms {

/**
 * @brief A class that knows how to calculate fluxes in an array of circular apertures
 * @ingroup meas/algorithms
 */
template<typename MaskedImageT>
class measureApertureArrayPhotometry : public measurePhotometry<MaskedImageT> {
public:
    /**
     * @brief Return the (unique) instance of measureApertureArrayPhotometry
     */
    static measurePhotometry<MaskedImageT>* getInstance(float const radius) {
        if (_instance == NULL) {
            _instance = new measureApertureArrayPhotometry(radius);
            measurePhotometry<MaskedImageT>::registerType("APERTUREARRAY", APERTUREARRAY);
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
