// -*- lsst-c++ -*-
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

#ifndef LSST_MEAS_ALGORITHMS_CorrelatedNoiseDetection_h_INCLUDED
#define LSST_MEAS_ALGORITHMS_CorrelatedNoiseDetection_h_INCLUDED

#include <string>
#include <vector>

#include "lsst/afw/image/MaskedImage.h"

namespace lsst { namespace meas { namespace algorithms {

/**
 * Estimate the noise correlation kernel of an image, assuming it is stationary.
 *
 * @param[in]  image          The image to be measured.  The empirical covariance
 *                            will be divided by the values in image's variance plane,
 *                            and its mask will be used to reject pixels containing
 *                            objects or artifacts.
 * @param[in]  radius         Distance in pixels to which the correlation is measured
 *                            on either side; the returned image will have dimensions
 *                            (2*radius + 1, 2*radius + 1).
 * @param[in]  badBitMask     A bit mask indicating pixels that should not be included
 *                            in the measurement.  Should generally include at least
 *                            DETECTED.
 *
 * @return  the noise correlation kernel: an image in which the central pixel
 *          represents the fraction of the total variance/covariance in the variance,
 *          and neighboring pixels contain the correlation at different offsets.
 */
afw::image::Image<float> measureCorrelationKernel(
    afw::image::MaskedImage<float> const & image,
    int radius,
    afw::image::MaskPixel badBitMask
);

/**
 * Estimate the noise correlation kernel of an image, assuming it is stationary.
 *
 * @param[in]  image          The image to be measured.  The empirical covariance
 *                            will be divided by the values in image's variance plane,
 *                            and its mask will be used to reject pixels containing
 *                            objects or artifacts.
 * @param[in]  radius         Distance in pixels to which the correlation is measured
 *                            on either side; the returned image will have dimensions
 *                            (2*radius + 1, 2*radius + 1).
 * @param[in]  badMaskPlanes  A list of mask planenames indicating pixels that
 *                            should not be included in the measurement.
 *                            Should generally include at least DETECTED.
 *
 * @return  the noise correlation kernel: an image in which the central pixel
 *          represents the fraction of the total variance/covariance in the variance,
 *          and neighboring pixels contain the correlation at different offsets.
 */
afw::image::Image<float> measureCorrelationKernel(
    afw::image::MaskedImage<float> const & image,
    int radius,
    std::vector<std::string> const & badMaskPlanes
);


/**
 * Fit an optimal detection kernel for a PSF that corrects for correlated noise.
 *
 * This is the same as the kernel that yields the PSF when convolved with the
 * correlation kernel.
 *
 * The returned kernel cannot be used in detection in quite the same way as the
 * PSF is used in detection on images with noise correlated noise.
 * In the uncorrelated noise case with image @f$z@f$, PSF @f$\phi@f$, and
 * per-pixel variance @f$\sigma^2@f$, the significance image is
 * @f[
 *    \nu = \frac{\phi \ast z}{\sigma \sqrt{\phi \cdot \phi}}
 * @f]
 * In the correlated noise case, with @f$\psi@f$ the kernel computed by this
 * function, the significance image is
 * @f[
 *    \nu = \frac{\psi \ast z}{\sigma \sqrt{\phi \cdot \psi}}
 * @f]
 * (note the difference in the denominator).
 *
 * @param[in]  psf          Image of the PSF model.
 * @param[in]  correlation  Noise correlation kernel, of the sort estimated by
 *                          measureCorrelationKernel.
 * @param[in]  radius       Radius of the detection kernel image in pixels.
 *                          The returned image will have dimensions
 *                          (2*radius + 1, 2*radius + 1).
 *
 * @return  an image containing the optimal detection kernel.
 */
afw::image::Image<double> fitGeneralDetectionKernel(
    afw::image::Image<double> const & psf,
    afw::image::Image<float> const & correlation,
    int radius
);


}}} // namespace lsst::meas::algorithms

#endif // !LSST_MEAS_ALGORITHMS_CorrelatedNoiseDetection_h_INCLUDED
