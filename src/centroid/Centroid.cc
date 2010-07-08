// -*- LSST-C++ -*-
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image/Image.h"

#include "lsst/meas/algorithms/Centroid.h"
#include "lsst/meas/algorithms/detail/MeasureFactory.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;

namespace lsst {
namespace meas {
namespace algorithms {

/**
 * @brief Call the concrete centroiding algorithm
 *
 * N.b. One purpose of this routine is to provide a place to specify default values for arguments
 */
template<typename ImageT>
Centroid MeasureCentroid<ImageT>::apply(ImageT const& image,
                                        int x,
                                        int y,
                                        lsst::afw::detection::Psf const* psf,
                                        double background
                                       ) const
{
    if (x - image.getX0() < 1 || x - image.getX0() > image.getWidth() - 2 ||
        y - image.getY0() < 1 || y - image.getY0() > image.getHeight() - 2) {
            throw LSST_EXCEPT(pexExceptions::RangeErrorException,
                              (boost::format("Object at (%d, %d) is too close "
                                             "to the edge of the frame") % x % y).str());
    }
    pexLogging::TTrace<8>("meas.algorithms.centroid", "Centroiding object at (%d, %d)", x, y);

    return doApply(image, x, y, psf, background);
}

/************************************************************************************************************/
//
// Explicit instantiations
// \cond
#define MAKE_CENTROIDERS(IMAGE_T)                                       \
    template class MeasureCentroid<IMAGE_T>;                            \
    template MeasureCentroid<IMAGE_T> *                                 \
    createMeasureProperty(std::string const&, IMAGE_T::ConstPtr, MeasureCentroid<IMAGE_T> const*);
                
MAKE_CENTROIDERS(lsst::afw::image::Image<int>)
MAKE_CENTROIDERS(lsst::afw::image::Image<float>)

// \endcond
                
}}}
