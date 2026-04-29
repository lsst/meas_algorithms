// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
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

#ifndef LSST_MEAS_ALGORITHMS_SpanSetMoments_h_INCLUDED
#define LSST_MEAS_ALGORITHMS_SpanSetMoments_h_INCLUDED

#include "ndarray.h"
#include "lsst/geom/Point.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/geom/ellipses/Quadrupole.h"
#include "lsst/afw/geom/SpanSet.h"
#include "lsst/shapelet/ShapeletFunction.h"

namespace lsst {
namespace meas {
namespace algorithms {

/**
 *  A struct that computes the unweighted moments of the pixels in an
 * `lsst.afw.geom.SpanSet`.
 *
 * This class provides low-level pixel-processing code for
 * `ComputeRoughPsfShapeletsTask`.
 */
struct SpanSetMoments {

    /// Total flux within the SpanSet.
    double flux;

    /// Total variance within the SpanSet.
    double variance;

    /// Center derived from the unweighted first moment of the image.
    geom::Point2D center;

    /// Shape derived from the unweighted second moment of the image.
    afw::geom::ellipses::Quadrupole shape;

    /// The pixels actually used to compute the moments.
    std::shared_ptr<afw::geom::SpanSet> spans;

    /// Flag set if there were too many bad pixels to compute the moments.
    bool too_many_bad_pixels = false;

    /// Flag set if the center did not lie within the SpanSet.
    bool center_out_of_bounds = false;

    /// Flag set if there was a bad pixel too close to the center.
    bool bad_pixel_in_center = false;

    /// Flag set if the second moments did not resolve to an ellipse.
    bool singular_second_moments = false;

    /// Test whether any failure flag is set.
    bool any_flags_set() const {
        return too_many_bad_pixels || center_out_of_bounds || bad_pixel_in_center || singular_second_moments;
    }

    /// Return a flattened array of the x coordinates in `spans`.
    ndarray::Array<float, 1, 1> get_x_array() const;

    /// Return a flattened array of the y coordinates in `spans`.
    ndarray::Array<float, 1, 1> get_y_array() const;

    /**
     * Compute the unweighted moments of an image within a SpanSet.
     *
     * @param[in] spans          Pixel region to use.
     * @param[in] masked_image   Image to measure.
     * @param[in] bad_bitmask    Mask of bad pixels to remove from `spans`
     *                           before computing the moments.
     * @param[in] bad_pixel_mask_fraction     Maximum fraction of the pixels
     *                                        in `spans` that can be bad
     *                                        before giving up.
     * @param[in] bad_pixel_exlusion_radius   Radius around the estimated
     *                                        center where the present of a
     *                                        bad pixels will cause the
     *                                        algorithm to give up.
     */
    static std::shared_ptr<SpanSetMoments> compute(afw::geom::SpanSet const& spans,
                                                   afw::image::MaskedImage<float> const& masked_image,
                                                   afw::image::MaskPixel bad_bitmask,
                                                   double bad_pixel_max_fraction,
                                                   double bad_pixel_exclusion_radius);

    /**
     * Fit a common shapelet expansion to multiple sources whose moments have
     * already been computed.
     *
     * The expected use case is a sample of stars with similar but not
     * identical moments, which might plausibly have an approximately common
     * shapelet representation while benefitting from allowing the ellipse
     * used for the expansion to vary (i.e. to partially capture spatial
     * variation in the PSF).
     *
     * @param[in] masked_image    Image to measure.
     * @param[in] moments         Moments already measured on stars.
     * @param[in] order           Order of the shapelet expansion.
     * @param[in] scale           Factor to scale the moments ellipses by to
     *                            yield the ellipse for the shapelet fit.
     * @param[in] circular        If true, use a circular basis with the trace
     *                            radius of the moments, instead of an
     *                            elliptical basis.
     */
    static shapelet::ShapeletFunction fit_shapelets(
            afw::image::MaskedImage<float> const& masked_image,
            std::vector<std::shared_ptr<SpanSetMoments>> const& moments,
            int order, double scale, bool circular);
};

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst

#endif  // !LSST_MEAS_ALGORITHMS_SpanSetMoments_h_INCLUDED
