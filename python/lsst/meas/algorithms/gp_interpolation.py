# This file is part of meas_algorithms.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
from lsst.meas.algorithms import CloughTocher2DInterpolatorUtils as ctUtils
from lsst.geom import Box2I, Point2I
from lsst.afw.geom import SpanSet
import copy
import treegp

import logging

__all__ = [
    "InterpolateOverDefectGaussianProcess",
    "GaussianProcessTreegp",
]


def updateMaskFromArray(mask, bad_pixel, interpBit):
    """
    Update the mask array with the given bad pixels.

    Parameters
    ----------
    mask : `lsst.afw.image.MaskedImage`
        The mask image to update.
    bad_pixel : `np.array`
        An array-like object containing the coordinates of the bad pixels.
        Each row should contain the x and y coordinates of a bad pixel.
    interpBit : `int`
        The bit value to set for the bad pixels in the mask.
    """
    x0 = mask.getX0()
    y0 = mask.getY0()
    for row in bad_pixel:
        x = int(row[0] - x0)
        y = int(row[1] - y0)
        mask.array[y, x] |= interpBit
    # TO DO --> might be better: mask.array[int(bad_pixel[:,1]-y0), int(bad_pixel[:,0]-x)] |= interpBit


def median_with_mad_clipping(data, mad_multiplier=2.0):
    """
    Calculate the median of the input data after applying Median Absolute Deviation (MAD) clipping.

    The MAD clipping method is used to remove outliers from the data. The median of the data is calculated,
    and then the MAD is calculated as the median absolute deviation from the median. The data is then clipped
    by removing values that are outside the range of median +/- mad_multiplier * MAD. Finally, the median of
    the clipped data is returned.

    Parameters:
    -----------
    data : `np.array`
        Input data array.
    mad_multiplier : `float`, optional
        Multiplier for the MAD value used for clipping. Default is 2.0.

    Returns:
    --------
    median_clipped : `float`
        Median value of the clipped data.

    Examples:
    ---------
    >>> data = [1, 2, 3, 4, 5, 100]
    >>> median_with_mad_clipping(data)
    3.5
    """
    median = np.median(data)
    mad = np.median(np.abs(data - median))
    clipping_range = mad_multiplier * mad
    clipped_data = np.clip(data, median - clipping_range, median + clipping_range)
    median_clipped = np.median(clipped_data)
    return median_clipped


class GaussianProcessTreegp:
    """
    Gaussian Process Treegp class for Gaussian Process interpolation.

    The basic GP regression, which uses Cholesky decomposition.

    Parameters:
    -----------
    std : `float`, optional
        Standard deviation of the Gaussian Process kernel. Default is 1.0.
    correlation_length : `float`, optional
        Correlation length of the Gaussian Process kernel. Default is 1.0.
    white_noise : `float`, optional
        White noise level of the Gaussian Process. Default is 0.0.
    mean : `float`, optional
        Mean value of the Gaussian Process. Default is 0.0.
    """

    def __init__(self, std=1.0, correlation_length=1.0, white_noise=0.0, mean=0.0):
        self.std = std
        self.correlation_length = correlation_length
        self.white_noise = white_noise
        self.mean = mean

        # Looks like weird to do that, but this is justified.
        # in GP if no noise is provided, even if matrix
        # can be inverted, it wont invert because of numerical
        # issue (det(K)~0). Add a little bit of noise allow
        # to compute a numerical solution in the case of no
        # external noise is added. Wont happened on real
        # image but help for unit test.
        if self.white_noise == 0.0:
            self.white_noise = 1e-5

    def fit(self, x_train, y_train):
        """
        Fit the Gaussian Process to the given training data.

        Parameters:
        -----------
        x_train : `np.array`
            Input features for the training data.
        y_train : `np.array`
            Target values for the training data.
        """
        kernel = f"{self.std}**2 * RBF({self.correlation_length})"
        self.gp = treegp.GPInterpolation(
            kernel=kernel,
            optimizer="none",
            normalize=False,
            white_noise=self.white_noise,
        )
        self.gp.initialize(x_train, y_train - self.mean)
        self.gp.solve()

    def predict(self, x_predict):
        """
        Predict the target values for the given input features.

        Parameters:
        -----------
        x_predict : `np.array`
            Input features for the prediction.

        Returns:
        --------
        y_pred : `np.array`
            Predicted target values.
        """
        y_pred = self.gp.predict(x_predict)
        return y_pred + self.mean


class InterpolateOverDefectGaussianProcess:
    """
    InterpolateOverDefectGaussianProcess class performs Gaussian Process
    (GP) interpolation over defects in an image.

    Parameters:
    -----------
    masked_image : `lsst.afw.image.MaskedImage`
        The masked image containing the defects to be interpolated.
    defects : `list`[`str`], optional
        The types of defects to be interpolated. Default is ["SAT"].
    fwhm : `float`, optional
        The full width at half maximum (FWHM) of the PSF. Default is 5.
    bin_spacing : `int`, optional
        The spacing between bins for good pixel binning. Default is 10.
    threshold_dynamic_binning : `int`, optional
        The threshold for dynamic binning. Default is 1000.
    threshold_subdivide : `int`, optional
        The threshold for sub-dividing the bad pixel array to avoid memory error. Default is 20000.
    correlation_length_cut : `int`, optional
        The factor by which to dilate the bounding box around defects. Default is 5.
    log : `lsst.log.Log`, `logging.Logger` or `None`, optional
        Logger object used to write out messages. If `None` a default
        logger will be used.
    """

    def __init__(
        self,
        masked_image,
        defects=["SAT"],
        fwhm=5,
        bin_image=True,
        bin_spacing=10,
        threshold_dynamic_binning=1000,
        threshold_subdivide=20000,
        correlation_length_cut=5,
        log=None,
    ):

        self.log = log or logging.getLogger(__name__)

        self.bin_image = bin_image
        self.bin_spacing = bin_spacing
        self.threshold_subdivide = threshold_subdivide
        self.threshold_dynamic_binning = threshold_dynamic_binning

        self.masked_image = masked_image
        self.defects = defects
        self.correlation_length = fwhm
        self.correlation_length_cut = correlation_length_cut

        self.interpBit = self.masked_image.mask.getPlaneBitMask("INTRP")

    def run(self):
        """
        Interpolate over the defects in the image.

        Change self.masked_image .
        """
        if self.defects == [] or self.defects is None:
            self.log.info("No defects found. No interpolation performed.")
        else:
            mask = self.masked_image.getMask()
            bad_pixel_mask = mask.getPlaneBitMask(self.defects)
            bad_mask_span_set = SpanSet.fromMask(mask, bad_pixel_mask).split()

            bbox = self.masked_image.getBBox()
            global_xmin, global_xmax = bbox.minX, bbox.maxX
            global_ymin, global_ymax = bbox.minY, bbox.maxY

            for spanset in bad_mask_span_set:
                bbox = spanset.getBBox()
                # Dilate the bbox to make sure we have enough good pixels around the defect
                # For now, we dilate by 5 times the correlation length
                # For GP with the isotropic kernel, points at the default value of
                # correlation_length_cut=5 have negligible effect on the prediction.
                bbox = bbox.dilatedBy(
                    int(self.correlation_length * self.correlation_length_cut)
                )  # need integer as input.
                xmin, xmax = max([global_xmin, bbox.minX]), min(global_xmax, bbox.maxX)
                ymin, ymax = max([global_ymin, bbox.minY]), min(global_ymax, bbox.maxY)
                localBox = Box2I(Point2I(xmin, ymin), Point2I(xmax - xmin, ymax - ymin))
                masked_sub_image = self.masked_image[localBox]

                masked_sub_image = self.interpolate_masked_sub_image(masked_sub_image)
                self.masked_image[localBox] = masked_sub_image

    def _good_pixel_binning(self, pixels):
        """
        Performs pixel binning using treegp.meanify

        Parameters:
        -----------
        pixels : `np.array`
            The array of pixels.

        Returns:
        --------
        `np.array`
            The binned array of pixels.
        """

        n_pixels = len(pixels[:, 0])
        dynamic_binning = int(np.sqrt(n_pixels / self.threshold_dynamic_binning))
        if n_pixels / self.bin_spacing**2 < n_pixels / dynamic_binning**2:
            bin_spacing = self.bin_spacing
        else:
            bin_spacing = dynamic_binning
        binning = treegp.meanify(bin_spacing=bin_spacing, statistics="mean")
        binning.add_field(
            pixels[:, :2],
            pixels[:, 2:].T,
        )
        binning.meanify()
        return np.array(
            [binning.coords0[:, 0], binning.coords0[:, 1], binning.params0]
        ).T

    def interpolate_masked_sub_image(self, masked_sub_image):
        """
        Interpolate the masked sub-image.

        Parameters:
        -----------
        masked_sub_image : `lsst.afw.image.MaskedImage`
            The sub-masked image to be interpolated.

        Returns:
        --------
        `lsst.afw.image.MaskedImage`
            The interpolated sub-masked image.
        """

        cut = int(
            self.correlation_length * self.correlation_length_cut
        )  # need integer as input.
        bad_pixel, good_pixel = ctUtils.findGoodPixelsAroundBadPixels(
            masked_sub_image, self.defects, buffer=cut
        )
        # Do nothing if bad pixel is None.
        if bad_pixel.size == 0 or good_pixel.size == 0:
            self.log.info("No bad or good pixels found. No interpolation performed.")
            return masked_sub_image
        # Do GP interpolation if bad pixel found.
        else:
            # gp interpolation
            sub_image_array = masked_sub_image.getVariance().array
            white_noise = np.sqrt(
                np.mean(sub_image_array[np.isfinite(sub_image_array)])
            )
            kernel_amplitude = np.max(good_pixel[:, 2:])
            if not np.isfinite(kernel_amplitude):
                filter_finite = np.isfinite(good_pixel[:, 2:]).T[0]
                good_pixel = good_pixel[filter_finite]
                if good_pixel.size == 0:
                    self.log.info(
                        "No bad or good pixels found. No interpolation performed."
                    )
                    return masked_sub_image
                # kernel amplitude might be better described by maximum value of good pixel given
                # the data and not really a random gaussian field.
                kernel_amplitude = np.max(good_pixel[:, 2:])

            if self.bin_image:
                try:
                    good_pixel = self._good_pixel_binning(copy.deepcopy(good_pixel))
                except Exception:
                    self.log.info(
                        "Binning failed, use original good pixel array in interpolation."
                    )

            # put this after binning as computing median is O(n*log(n))
            clipped_median = median_with_mad_clipping(good_pixel[:, 2:])

            gp = GaussianProcessTreegp(
                std=np.sqrt(kernel_amplitude),
                correlation_length=self.correlation_length,
                white_noise=white_noise,
                mean=clipped_median,
            )
            gp.fit(good_pixel[:, :2], np.squeeze(good_pixel[:, 2:]))
            if bad_pixel.size < self.threshold_subdivide:
                gp_predict = gp.predict(bad_pixel[:, :2])
                bad_pixel[:, 2:] = gp_predict.reshape(np.shape(bad_pixel[:, 2:]))
            else:
                self.log.info("sub-divide bad pixel array to avoid memory error.")
                for i in range(0, len(bad_pixel), self.threshold_subdivide):
                    end = min(i + self.threshold_subdivide, len(bad_pixel))
                    gp_predict = gp.predict(bad_pixel[i:end, :2])
                    bad_pixel[i:end, 2:] = gp_predict.reshape(
                        np.shape(bad_pixel[i:end, 2:])
                    )

            # Update values
            ctUtils.updateImageFromArray(masked_sub_image.image, bad_pixel)
            updateMaskFromArray(masked_sub_image.mask, bad_pixel, self.interpBit)
            return masked_sub_image
