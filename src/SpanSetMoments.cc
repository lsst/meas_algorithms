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

#include <limits>
#include "Eigen/Dense"
#include "lsst/shapelet/constants.h"
#include "ndarray/eigen.h"
#include "lsst/meas/algorithms/SpanSetMoments.h"
#include <Eigen/src/Core/util/Constants.h>
#include "lsst/shapelet/MatrixBuilder.h"

namespace lsst {
namespace meas {
namespace algorithms {

ndarray::Array<float, 1, 1> SpanSetMoments::get_x_array() const {
    ndarray::Array<float, 1, 1> result = ndarray::allocate(spans->getArea());
    spans->applyFunctor([](geom::Point2I const& point, float& out) { out = point.getX(); },
                        ndarray::ndFlat(result));
    return result;
}

ndarray::Array<float, 1, 1> SpanSetMoments::get_y_array() const {
    ndarray::Array<float, 1, 1> result = ndarray::allocate(spans->getArea());
    spans->applyFunctor([](geom::Point2I const& point, float& out) { out = point.getY(); },
                        ndarray::ndFlat(result));
    return result;
}

std::shared_ptr<SpanSetMoments> SpanSetMoments::compute(afw::geom::SpanSet const& spans,
                                                        afw::image::MaskedImage<float> const& masked_image,
                                                        afw::image::MaskPixel bad_bitmask,
                                                        double bad_pixel_max_fraction,
                                                        double bad_pixel_exclusion_radius) {
    double const nan = std::numeric_limits<double>::quiet_NaN();
    std::shared_ptr<SpanSetMoments> result(
            new SpanSetMoments{nan, nan, geom::Point2D(), afw::geom::ellipses::Quadrupole()});
    auto area = spans.getArea();
    auto bbox = spans.getBBox();
    auto bad_spans = afw::geom::SpanSet::fromMask((*masked_image.getMask())[bbox], bad_bitmask);
    result->spans = spans.intersectNot(*bad_spans);
    if ((area - result->spans->getArea()) > (bad_pixel_max_fraction * area)) {
        result->too_many_bad_pixels = true;
        return result;
    }
    // Compute the variance of the flux.
    result->variance = 0.0;
    result->spans->applyFunctor(
            [&result](geom::Point2I const& point, float data) { result->variance += data; },
            *masked_image.getVariance());
    // Compute the flux (zeroth moment) and center (first moments).
    result->flux = 0.0;
    double qx = 0.0;
    double qy = 0.0;
    result->spans->applyFunctor(
            [&result, &qx, &qy](geom::Point2I const& point, float data) {
                result->flux += data;
                qx += point.getX() * data;
                qy += point.getY() * data;
            },
            *masked_image.getImage());
    result->center.setX(qx / result->flux);
    result->center.setY(qy / result->flux);
    // Check to see if the center is plausible, and that it is not too close
    // to a bad pixel.
    if (!spans.contains(geom::Point2I(result->center))) {
        result->center_out_of_bounds = true;
        return result;
    }
    double min_bad_pixel_center_distance_squared = std::numeric_limits<double>::infinity();
    for (auto const& span : *bad_spans) {
        for (auto point : span) {
            min_bad_pixel_center_distance_squared =
                    std::min(geom::Point2D(point).distanceSquared(result->center),
                             min_bad_pixel_center_distance_squared);
        }
    }
    if (min_bad_pixel_center_distance_squared < bad_pixel_exclusion_radius * bad_pixel_exclusion_radius) {
        result->bad_pixel_in_center = true;
        return result;
    }
    // Compute the second moments.
    double qxx = 0.0;
    double qyy = 0.0;
    double qxy = 0.0;
    result->spans->applyFunctor(
            [&result, &qxx, &qyy, &qxy](geom::Point2I const& point, float data) {
                double dx = point.getX() - result->center.getX();
                double dy = point.getY() - result->center.getY();
                qxx += dx * dx * data;
                qyy += dy * dy * data;
                qxy += dx * dy * data;
            },
            *masked_image.getImage());
    result->shape.setIxx(qxx / result->flux);
    result->shape.setIyy(qyy / result->flux);
    result->shape.setIxy(qxy / result->flux);
    if (result->shape.getDeterminant() <= 0.0 || !std::isfinite(result->shape.getTraceRadius())) {
        result->singular_second_moments = true;
        return result;
    }
    return result;
}

shapelet::ShapeletFunction SpanSetMoments::fit_shapelets(
        afw::image::MaskedImage<float> const& masked_image,
        std::vector<std::shared_ptr<SpanSetMoments>> const& sources, int order, double scale,
        bool circular) {
    std::size_t n_pixels = std::accumulate(sources.begin(), sources.end(), 0,
                                           [](std::size_t n, std::shared_ptr<SpanSetMoments> const& source) {
                                               return n + source->spans->getArea();
                                           });
    std::size_t n_coeff = shapelet::computeSize(order);
    ndarray::Array<float, 1, 1> data_array = ndarray::allocate(n_pixels);
    data_array.deep() = 0.0;
    // We allocate the design matrix via its transpose, because that yields
    // the memory order expected by the shapelet library and the better
    // memory order for the least-squares solver.
    ndarray::Array<float, 2, -2> design_matrix =
            ndarray::Array<float, 2, 2>(ndarray::allocate(n_coeff, n_pixels)).transpose();
    design_matrix.deep() = 0.0;
    std::size_t start = 0;
    for (auto const& source : sources) {
        shapelet::MatrixBuilder<float> builder(source->get_x_array(), source->get_y_array(), order);
        std::size_t stop = start + builder.getDataSize();
        auto source_data_array = data_array[ndarray::view(start, stop)].shallow();
        source->spans->flatten(source_data_array, masked_image.getImage()->getArray(), masked_image.getXY0());
        auto source_design_matrix = design_matrix[ndarray::view(start, stop)()].shallow();
        afw::geom::ellipses::Ellipse ellipse(source->shape, source->center);
        if (circular) {
            double radius = ellipse.getCore().getTraceRadius();
            ellipse.setCore(afw::geom::ellipses::Axes(radius, radius, 0.0));
        }
        ellipse.scale(scale);
        builder(source_design_matrix, ellipse);
        asEigenArray(source_design_matrix) *= source->flux;
        start = stop;
    }
    shapelet::ShapeletFunction result(order, shapelet::HERMITE);
    asEigenMatrix(result.getCoefficients()) = asEigenMatrix(design_matrix)
                                                      .bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
                                                      .solve(asEigenMatrix(data_array))
                                                      .cast<double>();
    return result;
}

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst
