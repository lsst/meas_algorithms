/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
// -*- LSST-C++ -*-
/**
 * @file
 */

#include "boost/make_shared.hpp"

#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/meas/algorithms/CentroidControl.h"
#include "all.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace afwDet = lsst::afw::detection;
namespace afwImage = lsst::afw::image;

namespace lsst {
namespace meas {
namespace algorithms {

namespace {

/**
 * @brief A class that knows how to calculate centroids as a simple unweighted first moment
 * of the 3x3 region around a pixel
 */
class GaussianCentroid : public CentroidAlgorithm {
public:

    GaussianCentroid(
        GaussianCentroidControl const & ctrl,
        afw::table::Schema & schema
    ) : CentroidAlgorithm(ctrl, schema, "centroid measured with weighted Gaussian")
    {}

private:
    
    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(GaussianCentroid);
};

/**
 * @brief Given an image and a pixel position, calculate a position using a Gaussian fit
 */
template<typename PixelT>
void GaussianCentroid::_apply(
    afw::table::SourceRecord & source, 
    afw::image::Exposure<PixelT> const& exposure,
    afw::geom::Point2D const & center
) const {
    source.set(getKeys().meas, center); // better than NaN
    source.set(getKeys().flag, true); // say we've failed so that's the result if we throw
    typedef afw::image::Image<PixelT> ImageT;
    ImageT const& image = *exposure.getMaskedImage().getImage();

    int x = static_cast<int>(center.getX() + 0.5);
    int y = static_cast<int>(center.getY() + 0.5);

    x -= image.getX0();                 // work in image Pixel coordinates
    y -= image.getY0();

    FittedModel fit = twodg(image, x, y); // here's the fitter

    if (fit.params[FittedModel::PEAK] <= 0) {
        throw LSST_EXCEPT(pexExceptions::RuntimeErrorException,
                          (boost::format("Object at (%d, %d) has a peak of %g") %
                           x % y % fit.params[FittedModel::PEAK]).str());
    }

    source.set(getKeys().flag, false);
    source.set(
        getKeys().meas, 
        afw::geom::Point2D(
            lsst::afw::image::indexToPosition(image.getX0()) + fit.params[FittedModel::X0],
            lsst::afw::image::indexToPosition(image.getY0()) + fit.params[FittedModel::Y0]
        )
    );
    // FIXME: should report uncertainty
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(GaussianCentroid);

} // anonymous

PTR(AlgorithmControl) GaussianCentroidControl::_clone() const {
    return boost::make_shared<GaussianCentroidControl>(*this);
}

PTR(Algorithm) GaussianCentroidControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const &
) const {
    return boost::make_shared<GaussianCentroid>(*this, boost::ref(schema));
}

}}} // namespace lsst::meas::algorithms
