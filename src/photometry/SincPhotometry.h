#if !defined(LSST_MEAS_ALGORITHMS_SINCPHOTOMETRY_H)
#define LSST_MEAS_ALGORITHMS_SINCPHOTOMETRY_H 1
/**
 * @file SincPhotometry.h
 * @brief Compute Aperture and PSF fluxes in a simple way
 * @ingroup meas/algorithms
 * @author Steve Bickerton (adapted from RHL's Shape class)
 */
#include "lsst/meas/algorithms/Photometry.h"
#include "lsst/meas/algorithms/PhotometryImpl.h"

// for debug only ... delete
#include "lsst/afw/image.h"

namespace lsst { namespace meas { namespace algorithms {

/// for debug only ... delete
template<typename PixelT>
typename lsst::afw::image::Image<PixelT>::Ptr getCoeffImage(double const xcen0, double const ycen0, double const radius);

            
/**
 * @brief A class that knows how to calculate fluxes using the SINC photometry algorithm
 * @ingroup meas/algorithms
 */
template<typename MaskedImageT>
class measureSincPhotometry : public measurePhotometry<MaskedImageT> {
public:
    /**
     * @brief Return the (unique) instance of measureSincPhotometry
     */
    static measurePhotometry<MaskedImageT>* getInstance(double const radius) {
        if (_instance == NULL) {
            _instance = new measureSincPhotometry(radius);
            measurePhotometry<MaskedImageT>::registerType("SINC", SINC);
        }
        _instance->setRadius(radius);
        return _instance;
    }
private:
    measureSincPhotometry(double const radius) :  measurePhotometry<MaskedImageT>(radius) {}
    Photometry doApply(MaskedImageT const& image, double xcen, double ycen, PSF const *, double background) const;

    static measureSincPhotometry* _instance;
};
            
}}}
#endif
