// -*- LSST-C++ -*-

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
 
/**
 * @file
 */

#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/meas/algorithms/CentroidControl.h"

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
class NaiveCentroid : public CentroidAlgorithm {
public:

    NaiveCentroid(NaiveCentroidControl const & ctrl, afw::table::Schema & schema) :
        CentroidAlgorithm(ctrl, schema, "unweighted 3x3 first moment centroid")
    {}

private:

    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(NaiveCentroid);
};
    
/**
 * @brief Given an image and a pixel position, return a Centroid using a naive 3x3 weighted moment
 */
template<typename PixelT>
void NaiveCentroid::_apply(
    afw::table::SourceRecord & source, 
    afw::image::Exposure<PixelT> const& exposure,
    afw::geom::Point2D const & center
) const {
    source.set(getKeys().flag, true); // say we've failed so that's the result if we throw
    typedef afw::image::Image<PixelT> ImageT;
    ImageT const& image = *exposure.getMaskedImage().getImage();

    int x = center.getX();  // FIXME: this is different from GaussianCentroid and SdssCentroid here,
    int y = center.getY();  //        and probably shouldn't be.

    x -= image.getX0();                 // work in image Pixel coordinates
    y -= image.getY0();

    if (x < 1 || x >= image.getWidth() - 1 || y < 1 || y >= image.getHeight() - 1) {
         throw LSST_EXCEPT(lsst::pex::exceptions::LengthErrorException,
                           (boost::format("Object at (%d, %d) is too close to the edge") % x % y).str());
    }

    typename ImageT::xy_locator im = image.xy_at(x, y);

    double const sum =
        (im(-1,  1) + im( 0,  1) + im( 1,  1) +
         im(-1,  0) + im( 0,  0) + im( 1,  0) +
         im(-1, -1) + im( 0, -1) + im( 1, -1))
        - 9 * static_cast<NaiveCentroidControl const &>(getControl()).background;

    if (sum == 0.0) {
        throw LSST_EXCEPT(pexExceptions::RuntimeErrorException,
                          (boost::format("Object at (%d, %d) has no counts") %
                           x % y).str());
    }

    double const sum_x =
        -im(-1,  1) + im( 1,  1) +
        -im(-1,  0) + im( 1,  0) +
        -im(-1, -1) + im( 1, -1);
    double const sum_y =
        (im(-1,  1) + im( 0,  1) + im( 1,  1)) -
        (im(-1, -1) + im( 0, -1) + im( 1, -1));

    source.set(getKeys().flag, false);
    source.set(
        getKeys().meas, 
        afw::geom::Point2D(
            lsst::afw::image::indexToPosition(x + image.getX0()) + sum_x / sum,
            lsst::afw::image::indexToPosition(y + image.getY0()) + sum_y / sum
        )
    );
    // FIXME: should report uncertainty
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(NaiveCentroid);

} // anonymous

PTR(AlgorithmControl) NaiveCentroidControl::_clone() const {
    return boost::make_shared<NaiveCentroidControl>(*this);
}

PTR(Algorithm) NaiveCentroidControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const &
) const {
    return boost::make_shared<NaiveCentroid>(*this, boost::ref(schema));
}

}}}  // namespace lsst::meas::algorithms
