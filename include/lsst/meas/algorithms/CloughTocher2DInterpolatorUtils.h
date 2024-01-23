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

#ifndef LSST_MEAS_ALGORITHMS_CLOUGHTOCHER2DINTERPOLATORUTILS_H
#define LSST_MEAS_ALGORITHMS_CLOUGHTOCHER2DINTERPOLATORUTILS_H

/** \file
 * \brief Provide pixel-level utilities to support CloughTocher2DInterpolator.
 */

#include <vector>
#include "lsst/afw/image/Image.h"
#include "lsst/afw/image/MaskedImage.h"

namespace lsst {
namespace meas {
namespace algorithms {

/**
@brief This class contains static utility methods that are used by the CloughTocher2DInterpolatorTask
and exists solely to provide a namespace.
*/
class CloughTocher2DInterpolatorUtils {
    typedef float PixelT;

public:
    /**
     * @brief Find the positions of bad pixels as defined by the masks, and also
     * find the location and pixel value of good pixels around the bad pixels.
     *
     * @param[in] mimage MaskedImage to find the pixels from.
     * @param[in] maskPlanes List of mask planes to consider as bad pixels.
     * @param[in] buffer Number of pixels to dilate the bad pixels by.
     * @return A pair of arrays, the first containing the bad pixels and the second containing the good
     * pixels.
     *
     * @note The bad pixels array has shape (N, 3) where N is the number of bad pixels. The first column is
     * the x coordinate, the second column is the y coordinate, and the third column is the pixel value.
     * The good pixels array has shape (M, 3) where M is the number of good pixels and has the same format as
     * the previous one.
     */
    static std::pair<ndarray::Array<float, 2, 2>, ndarray::Array<float, 2, 2>> findGoodPixelsAroundBadPixels(
            afw::image::MaskedImage<PixelT, afw::image::MaskPixel, afw::image::VariancePixel> const &mimage,
            std::vector<std::string> const &maskPlanes, int const buffer);

    /**
     * @brief Update the values (third column) of the pixelArray from the image.
     *
     * @param[out] pixelArray The three-column array to update.
     * @param[in] image The image to update the pixel values from.
     *
     * @note The pixelArray is typically one of the arrays returned by findGoodPixelsAroundBadPixels.
     *
     * @see updateImageFromArray
     */
    static void updateArrayFromImage(ndarray::Array<float, 2, 2> const &pixelArray,
                                     afw::image::Image<PixelT> const &image);

    /**
     * @brief Update the pixel values of the image from the pixelArray
     *
     * @param[out] image The image to update.
     * @param[in] pixelArray The three-column array to update the pixel values from.
     *
     * @note The pixelArray is typically one of the arrays returned by findGoodPixelsAroundBadPixels.
     *
     * @see updateArrayFromImage
     */
    static void updateImageFromArray(afw::image::Image<PixelT> &image,
                                     ndarray::Array<float const, 2, 2> const &pixelArray);
};

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst

#endif
