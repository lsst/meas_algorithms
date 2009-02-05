#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"

#include "lsst/meas/algorithms/SdssCentroid.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;

namespace lsst { namespace meas { namespace algorithms {

/**
 * @brief the (unique) instance of SdssCentroider
 */
template<typename ImageT> SdssCentroider<ImageT>* SdssCentroider<ImageT>::_instance = 0;

/**
 * @brief Given an image and a pixel position, return a Centroid using the SDSS algorithm
 */
template<typename ImageT>
Centroid SdssCentroider<ImageT>::doApply(ImageT const& image, ///< The Image wherein dwells the object
                                         int x,               ///< object's column position
                                         int y,               ///< object's row position
                                         PSF const*,          ///< image's PSF
                                         double background    ///< image's background level
                                        ) const {
    typename ImageT::xy_locator im = image.xy_at(x, y);

    double const sum =
        (im(-1,  1) + im( 0,  1) + im( 1,  1) +
         im(-1,  0) + im( 0,  0) + im( 1,  0) +
         im(-1, -1) + im( 0, -1) + im( 1, -1)) - 9*background;

    if (sum == 0.0) {
        throw LSST_EXCEPT(pexExceptions::UnderflowErrorException,
                          (boost::format("Object at (%d, %d) has no counts") % x % y).str());
    }

    double const sum_x =
        -im(-1,  1) + im( 1,  1) +
        -im(-1,  0) + im( 1,  0) +
        -im(-1, -1) + im( 1, -1);
    double const sum_y =
        (im(-1,  1) + im( 0,  1) + im( 1,  1)) -
        (im(-1, -1) + im( 0, -1) + im( 1, -1));

    return Centroid(x + sum_x/sum, y + sum_y/sum);
}

//
// Explicit instantiations
//
// We need to make an instance here so as to register it with Centroider
//
// \cond
#define MAKE_CENTROIDERS(IMAGE_T) \
                namespace { \
                    Centroider<lsst::afw::image::Image<IMAGE_T> >* foo = \
                        SdssCentroider<lsst::afw::image::Image<IMAGE_T> >::getInstance(); \
                }
                
MAKE_CENTROIDERS(float)


// \endcond

}}}
