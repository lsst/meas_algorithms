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
template<typename ExposureT>
class NaiveCentroid : public Algorithm<ExposureT>
{
public:
    typedef Algorithm<ExposureT> AlgorithmT;

    NaiveCentroid(NaiveCentroidControl const & ctrl, afw::table::Schema & schema) :
        AlgorithmT(ctrl), _background(ctrl.background),
        _keys(addCentroidFields(schema, ctrl.name, "unweighted 3x3 first moment centroid"))
    {}

    virtual void apply(afw::table::SourceRecord &, ExposurePatch<ExposureT> const&) const;

private:
    double _background;
    afw::table::KeyTuple<afw::table::Centroid> _keys;
};
    
/**
 * @brief Given an image and a pixel position, return a Centroid using a naive 3x3 weighted moment
 */
template<typename ExposureT>
void NaiveCentroid<ExposureT>::apply(
    afw::table::SourceRecord & source,
    ExposurePatch<ExposureT> const& patch
) const {
    CONST_PTR(ExposureT) exposure = patch.getExposure();
    typedef typename ExposureT::MaskedImageT::Image ImageT;
    ImageT const& image = *exposure->getMaskedImage().getImage();

    int x = patch.getCenter().getX();
    int y = patch.getCenter().getY();

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
         im(-1, -1) + im( 0, -1) + im( 1, -1)) - 9*_background;

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

    source.set(_keys.flag, true);
    source.set(
        _keys.meas, 
        afw::geom::Point2D(
            lsst::afw::image::indexToPosition(x + image.getX0()) + sum_x / sum,
            lsst::afw::image::indexToPosition(y + image.getY0()) + sum_y / sum
        )
    );
    // FIXME: should report uncertainty
}

} // anonymous

LSST_ALGORITHM_CONTROL_PRIVATE_IMPL(NaiveCentroidControl, NaiveCentroid)

}}}  // namespace lsst::meas::algorithms
