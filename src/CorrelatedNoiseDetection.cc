// -*- LSST-C++ -*-
/*
 * LSST Data Management System
 * Copyright 2017 LSST/AURA.
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

#include "lsst/pex/exceptions.h"
#include "lsst/afw/math/LeastSquares.h"
#include "lsst/meas/algorithms/CorrelatedNoiseDetection.h"

namespace lsst { namespace meas { namespace algorithms {


afw::image::Image<float> measureCorrelationKernel(
    afw::image::MaskedImage<float> const & mi,
    int radius,
    afw::image::MaskPixel badBitMask
) {
    afw::image::Image<float> result(2*radius + 1, 2*radius + 1);
    afw::image::Image<int> count(result.getDimensions());
    int const width = mi.getWidth();
    int const height = mi.getHeight();
    // iterate over pixels in the MaskedImage, skipping any that meet our bad mask criteria
    for (int y1 = 0; y1 < height; ++y1) {
        int const y2a = std::max(0, y1 - radius);
        int const y2b = std::min(height, y1 + radius + 1);
        for (int x1 = 0; x1 < width; ++x1) {
            if ((*mi.getMask())(x1, y1) & badBitMask) {
                continue;
            }
            float z = (*mi.getImage())(x1, y1) / (*mi.getVariance())(x1, y1);
            // iterate over neighboring pixels, with bounds set to avoid image boundaries
            int const x2a = std::max(0, x1 - radius);
            int const x2b = std::min(height, x1 + radius + 1);
            for (int y2 = y2a; y2 < y2b; ++y2) {
                auto miIter = mi.row_begin(y2) + x2a;
                auto outIter = result.row_begin(radius + y2 - y1) + radius + x2a - x1;
                auto countIter = count.row_begin(radius + y2 - y1) + radius + x2a - x1;
                for (int x2 = x2a; x2 < x2b; ++x2, ++miIter, ++outIter, ++countIter) {
                    if (miIter.mask() & badBitMask) {
                        continue;
                    }
                    *outIter += z*miIter.image();
                    *countIter += 1;
                }
            }
        }
    }
    result.getArray().deep() /= count.getArray();
    result.setXY0(-radius, -radius);
    return result;
}

afw::image::Image<float> measureCorrelationKernel(
    afw::image::MaskedImage<float> const & image,
    int radius,
    std::vector<std::string> const & badMaskPlanes
) {
    afw::image::MaskPixel mask = 0x0;
    for (auto plane : badMaskPlanes) {
        mask |= image.getMask()->getPlaneBitMask(plane);
    }
    return measureCorrelationKernel(image, radius, mask);
}


namespace {

afw::geom::Extent2I getOddBoxHalfWidth(afw::geom::Box2I const & box, std::string const & name) {
    afw::geom::Extent2I r((box.getWidth() - 1) / 2, (box.getHeight() - 1) / 2);
    if (r.getX()*2 + 1 != box.getWidth()) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            name + " image width must be an odd integer"
        );
    }
    if (r.getY()*2 + 1 != box.getHeight()) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            name + " image height must be an odd integer"
        );
    }
    return r;
}

} // anonymous


afw::image::Image<double> fitGeneralDetectionKernel(
    afw::image::Image<double> const & psf,
    afw::image::Image<float> const & correlation,
    int radius
) {
    auto const psfR = getOddBoxHalfWidth(psf.getBBox(), "PSF");
    auto const corrR = getOddBoxHalfWidth(correlation.getBBox(), "Correlation kernel");
    afw::geom::Extent2I const outR(radius, radius);
    if (psfR.getX() <= corrR.getX()) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            "PSF image width must be greater than correlation kernel width"
        );
    }
    if (psfR.getY() <= corrR.getY()) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            "PSF image height must be greater than correlation kernel height"
        );
    }
    if (psfR.getX() <= outR.getX()) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            "PSF image width must be greater than output kernel width"
        );
    }
    if (psfR.getY() <= outR.getY()) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            "PSF image height must be greater than output kernel height"
        );
    }

    int const psfN = (psfR.getX()*2 + 1)*(psfR.getY()*2 + 1);
    int const outN = (outR.getX()*2 + 1)*(outR.getY()*2 + 1);

    // Get locators at the center of each image, so we can use indices with
    // the origin at the center and make the code easier to read.
    auto psfL = psf.xy_at(psfR.getX(), psfR.getY());
    auto corrL = correlation.xy_at(corrR.getX(), corrR.getY());

    // rhs is the PSF image, with rows and columns flattened into one dimension
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(psfN);

    // matrix represents convolution with the correlated noise kernel image
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(psfN, outN);
    int xyN = 0;
    for (int y = -psfR.getY(); y <= psfR.getY(); ++y) {
        for (int x = -psfR.getX(); x <= psfR.getX(); ++x) {
            int ijN = 0;
            for (int i = -outR.getY(); i <= outR.getY(); ++i) {
                for (int j = -outR.getX(); j <= outR.getX(); ++j) {
                    // Could move these checks into the inner loop bounds for performance,
                    // but this is easier to read and probably fast enough.
                    if (std::abs(y - i) <= corrR.getY() && std::abs(x - j) <= corrR.getX()) {
                        matrix(xyN, ijN) = corrL(x - j, y - i);
                    }
                    ++ijN;
                }
            }
            rhs[xyN] = psfL(x, y);
            ++xyN;
        }
    }

    // solve for the kernel that produces the PSF when convolved with the
    // noise correlation kernel
    auto lstsq = afw::math::LeastSquares::fromDesignMatrix(matrix, rhs);
    auto solution = lstsq.getSolution();

    // copy the result from the flattened solution vector into an image
    afw::image::Image<double> result(outR.getX()*2 + 1, outR.getY()*2 + 1);
    auto outL = result.xy_at(outR.getX(), outR.getY());
    xyN = 0;
    for (int y = -outR.getY(); y <= outR.getY(); ++y) {
        for (int x = -outR.getX(); x <= outR.getX(); ++x) {
            outL(x, y) = solution[xyN];
            ++xyN;
        }
    }

    result.setXY0(-outR.getX(), -outR.getY());
    return result;
}


}}} // namespace lsst::meas::algorithms
