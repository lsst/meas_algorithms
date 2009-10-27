// -*- LSST-C++ -*-
#if !defined(LSST_MEAS_ALGORITHMS_SINCPHOTOMETRY_H)
#define LSST_MEAS_ALGORITHMS_SINCPHOTOMETRY_H 1
/**
 * @file SincPhotometry.h
 * @brief Compute Aperture and PSF fluxes in a simple way
 * @ingroup meas/algorithms
 * @author Steve Bickerton (adapted from RHL's Shape class)
 */
#include "lsst/meas/algorithms/Photometry.h"

// for debug only ... delete
#include "lsst/afw/image.h"

namespace lsst {
namespace meas {
namespace algorithms {

/// primarily for debug
template<typename PixelT>
typename lsst::afw::image::Image<PixelT>::Ptr getCoeffImage(double const xcen0,
                                                            double const ycen0,
                                                            double const radius);

            
/**
 * @brief A class that knows how to calculate fluxes using the SINC photometry algorithm
 * @ingroup meas/algorithms
 */
template<typename MaskedImageT>
class SincMeasurePhotometry : public MeasurePhotometry<MaskedImageT> {
public:
    static bool registerMe(std::string const& name);
protected:
    friend class MeasurePhotometryFactory<SincMeasurePhotometry>;
    SincMeasurePhotometry(float const radius) : MeasurePhotometry<MaskedImageT>(radius) {}
private:
    Photometry doApply(MaskedImageT const& image, double xcen, double ycen,
                       PSF const* psf, double background) const;
};

}}}
#endif
