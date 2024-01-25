// -*- LSST-C++ -*-
/*
 * This file is part of meas_algorithms.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
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
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file
 *
 * @brief Utility functions for CloughTocher2DInterpolatorTask
 *
 */

#include "lsst/afw/geom/SpanSet.h"
#include "lsst/afw/image/Image.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/geom/Point.h"
#include "lsst/meas/algorithms/CloughTocher2DInterpolatorUtils.h"

namespace lsst {
namespace meas {
namespace algorithms {

std::pair<ndarray::Array<float, 2, 2>, ndarray::Array<float, 2, 2>>
CloughTocher2DInterpolatorUtils::findGoodPixelsAroundBadPixels(afw::image::MaskedImage<float> const &mimage,
                                                               std::vector<std::string> const &maskPlanes,
                                                               int const buffer) {
    auto bitMask = afw::image::Mask<>::getPlaneBitMask(maskPlanes);

    auto badSpanSet = afw::geom::SpanSet::fromMask(
            *mimage.getMask(), [bitMask](afw::image::MaskPixel pixel) { return bool(pixel & bitMask); });
    ndarray::Array<float, 2, 2> badPixelArray = ndarray::allocate(badSpanSet->getArea(), 3);
    badSpanSet->applyFunctor(
            [](geom::Point2I pt, float &x, float &y, float &valueOut, float valueIn) {
                x = pt.getX();
                y = pt.getY();
                valueOut = valueIn;
            },
            ndFlat(badPixelArray[ndarray::view()(0)].shallow()),
            ndFlat(badPixelArray[ndarray::view()(1)].shallow()),
            ndFlat(badPixelArray[ndarray::view()(2)].shallow()),
            *mimage.getImage());

    auto allSpanSet = badSpanSet->dilated(buffer, afw::geom::Stencil::BOX);
    auto goodSpanSet = allSpanSet->intersectNot(*badSpanSet);
    goodSpanSet = goodSpanSet->clippedTo(mimage.getBBox());

    badSpanSet.reset();
    allSpanSet.reset();

    ndarray::Array<float, 2, 2> goodPixelArray = ndarray::allocate(goodSpanSet->getArea(), 3);
    goodSpanSet->applyFunctor(
            [](geom::Point2I pt, float &x, float &y, float &valueOut, float valueIn) {
                x = pt.getX();
                y = pt.getY();
                valueOut = valueIn;
            },
            ndFlat(goodPixelArray[ndarray::view()(0)].shallow()),
            ndFlat(goodPixelArray[ndarray::view()(1)].shallow()),
            ndFlat(goodPixelArray[ndarray::view()(2)].shallow()),
            *mimage.getImage());

    goodSpanSet.reset();

    return std::make_pair(badPixelArray, goodPixelArray);
}

void CloughTocher2DInterpolatorUtils::updateArrayFromImage(ndarray::Array<float, 2, 2> const &pixelArray,
                                                           afw::image::Image<float> const &image) {
    afw::image::CheckIndices docheck(true);
    auto x0 = image.getX0();
    auto y0 = image.getY0();
    for (auto row : pixelArray) {
        int x = row[0] - x0;
        int y = row[1] - y0;
        row[2] = image(x, y, docheck);
    }
}

void CloughTocher2DInterpolatorUtils::updateImageFromArray(
        afw::image::Image<float> &image, ndarray::Array<float const, 2, 2> const &pixelArray) {
    afw::image::CheckIndices docheck(true);
    auto x0 = image.getX0();
    auto y0 = image.getY0();
    for (auto row : pixelArray) {
        int x = row[0] - x0;
        int y = row[1] - y0;
        image(x, y, docheck) = row[2];
    }
}

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst
