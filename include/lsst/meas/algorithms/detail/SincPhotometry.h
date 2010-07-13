// -*- LSST-C++ -*-
#if !defined(LSST_MEAS_ALGORITHMS_SINCPHOTOMETRY_H)
#define LSST_MEAS_ALGORITHMS_SINCPHOTOMETRY_H 1
/**
 * @file SincPhotometry.h
 * @brief Compute Aperture using sinc interpolation
 * @ingroup meas/algorithms
 * @author Steve Bickerton (adapted from RHL's Shape class)
 */
namespace lsst {
namespace meas {
namespace algorithms {

/// primarily for debug
template<typename PixelT>
typename lsst::afw::image::Image<PixelT>::Ptr getCoeffImage(double const xcen0,
                                                            double const ycen0,
                                                            double const radius);

}}}
#endif
