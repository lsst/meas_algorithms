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

"""Gaussian Process interpolation over image defects.

The interpolation follows TURBO-GP (TUrbulence Removal By fourier-Optimized
Gaussian Process, RTN-129), the treegp implementation of the Gomes et al.
(2025, AJ 170:361) empirical kernel: the measured 2D 2-point correlation
function of the good pixels, cleaned by apodization and thresholding of its
Fourier power spectrum, is used directly as the kernel, so there are no
hyperparameters to set or fit. The linear algebra is done by
`treegp.GridConvolutionGP`: the covariance is never built, the solve is a
preconditioned conjugate gradient whose iterations cost one FFT pair on a
grid the size of the defect area, the kernel is positive semi-definite by
construction, and the prediction is a bilinear read-off of the posterior mean
field.

The 2-point correlation function is measured once per image on all the good
pixels (see `measure_2pcf_grid`), and the resulting kernel is shared by the
local solves done around each connected defect.
"""

import logging
import warnings

import numpy as np
import treecorr
import treegp
from scipy import fft

from lsst.afw.geom import SpanSet
from lsst.meas.algorithms import CloughTocher2DInterpolatorUtils as ctUtils

# We need to explicitly turn off multiprocessing in treecorr which is used
# by treegp.
treecorr.set_max_omp_threads(1)

__all__ = [
    "InterpolateOverDefectGaussianProcess",
    "GaussianProcessTreegp",
    "measure_2pcf_grid",
]


def updateMaskFromArray(mask, bad_pixel, interpBit):
    """Set a mask bit at the given pixel positions.

    Parameters
    ----------
    mask : `lsst.afw.image.Mask`
        The mask to update.
    bad_pixel : `numpy.ndarray`
        Array of shape (n, 3) whose first two columns are the x and y
        coordinates of the pixels (parent coordinates).
    interpBit : `int`
        The bit value to set for the given pixels.
    """
    x = np.rint(bad_pixel[:, 0]).astype(int) - mask.getX0()
    y = np.rint(bad_pixel[:, 1]).astype(int) - mask.getY0()
    mask.array[y, x] |= interpBit


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


def measure_2pcf_grid(image, good, max_sep):
    """Measure the 2D 2-point correlation function of a gridded image with
    gaps, by FFT.

    This is the pair-count estimator

        xi(dx, dy) = sum_{good pairs} z(x, y) z(x + dx, y + dy) / N_pairs(dx, dy)

    with unit weights, computed exactly at every integer lag as the ratio of
    the autocorrelation of the (zero-filled) image to the autocorrelation of
    the good-pixel mask, both evaluated by zero-padded FFTs. It includes the
    zero lag (the variance of the good pixels), which a pair counter such as
    treecorr excludes on gridded data since distinct pixels never coincide.
    The cost is O(N log N) in the number of pixels, independent of ``max_sep``.

    Parameters
    ----------
    image : `numpy.ndarray`
        Image array of shape (ny, nx), already mean-subtracted. Values at
        pixels that are not ``good`` are ignored (they may be NaN).
    good : `numpy.ndarray`
        Boolean array of shape (ny, nx), `True` for the pixels to use.
    max_sep : `int`
        Half width of the lag grid, in pixels.

    Returns
    -------
    xi : `numpy.ndarray`
        Correlation function on a (2 * max_sep, 2 * max_sep) grid of lags
        ``-max_sep .. max_sep - 1`` in each axis, indexed ``[iy, ix]`` with the
        zero lag at pixel ``max_sep`` (the treegp / treecorr TwoD layout, ready
        for `treegp.empirical_2pcf.clean` and `treegp.GridConvolutionGP`).
        Lags with no pair are set to zero.
    """
    image = np.asarray(image, dtype=float)
    good = np.asarray(good, dtype=bool)
    if image.ndim != 2 or image.shape != good.shape:
        raise ValueError(
            "image and good must be 2D arrays of the same shape. "
            f"Current shapes: {image.shape}, {good.shape}."
        )
    max_sep = int(max_sep)
    if max_sep < 1:
        raise ValueError(f"max_sep must be at least 1 pixel. Current value: {max_sep}.")

    z = np.where(good, image, 0.0)
    m = good.astype(float)
    ny, nx = image.shape
    # Pad by max_sep so that the circular autocorrelation equals the linear
    # one for every lag kept below.
    shape = (fft.next_fast_len(ny + max_sep, real=True), fft.next_fast_len(nx + max_sep, real=True))
    fz = fft.rfft2(z, s=shape)
    fm = fft.rfft2(m, s=shape)
    num = fft.irfft2(fz * np.conj(fz), s=shape)
    den = fft.irfft2(fm * np.conj(fm), s=shape)

    # Negative lags wrap around the end of the padded arrays.
    lags = np.arange(-max_sep, max_sep)
    rows = lags % shape[0]
    cols = lags % shape[1]
    num = num[np.ix_(rows, cols)]
    den = den[np.ix_(rows, cols)]

    xi = np.zeros_like(num)
    # den holds pair counts (integers up to FFT round-off).
    has_pairs = den > 0.5
    xi[has_pairs] = num[has_pairs] / den[has_pairs]
    return xi


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
    """Interpolate over the defects of an image with a Gaussian Process whose
    kernel is the measured 2-point correlation function of the image
    (TURBO-GP, RTN-129; Gomes et al. 2025, AJ 170:361).

    The kernel is measured once, on all the good pixels of the image
    (`measure_2pcf_grid`), and cleaned with `treegp.empirical_2pcf.clean`
    (apodization, thresholding of the Fourier power spectrum). Each
    connected defect is then interpolated from the good pixels within
    ``max_sep`` of it, with its own local mean (clipped median) and per-pixel
    noise (variance plane), using `treegp.GridConvolutionGP`: a matrix-free
    conjugate gradient solve on an FFT grid the size of the defect area, and
    a prediction that is a read-off of the posterior mean field.

    Parameters
    ----------
    masked_image : `lsst.afw.image.MaskedImage`
        The masked image containing the defects to be interpolated. Modified
        in place by `run`.
    defects : `list` [`str`], optional
        The mask planes to be interpolated over. Default is ["SAT"].
    fwhm : `float`, optional
        The full width at half maximum (FWHM) of the PSF, in pixels. Default
        is 5.
    kernel_half_width : `float`, optional
        Minimum half width, in pixels, of the 2-point correlation function
        grid, i.e. of the kernel support. Default is 30 (6 arcsec at 0.2
        arcsec/pixel, about 3 times the FWHM of the worst expected PSF).
    fwhm_factor : `float`, optional
        The kernel half width is at least ``fwhm_factor * fwhm``. Default is 3.
    power_threshold : `float`, optional
        Signal-to-noise threshold below which the Fourier modes of the
        measured correlation function are set to zero. Default is 2.5.
    apod_window : `str`, optional
        Apodization window applied to the measured correlation function
        before its Fourier transform, "hann" or "blackman-harris". Default is
        "hann".
    apod_radius : `float` or `None`, optional
        Radius, in pixels, where the apodization window reaches zero. If
        `None`, the window reaches zero at the edge of the kernel grid.
    white_noise : `float`, optional
        Additional noise (standard deviation, in image units) added in
        quadrature to the pixel errors from the variance plane. Default is 0.
    cg_rtol : `float`, optional
        Relative tolerance of the conjugate gradient solve. Default is 1e-4.
    cg_maxiter : `int`, optional
        Maximum number of conjugate gradient iterations per defect area.
        Default is 200.
    log : `lsst.log.Log`, `logging.Logger` or `None`, optional
        Logger object used to write out messages. If `None` a default
        logger will be used.
    bin_spacing, threshold_dynamic_binning, threshold_subdivide : optional
        Deprecated and ignored. They configured the binning of the good
        pixels and the chunking of the prediction that the former dense
        solver needed; `lsst.ip.isr` still passes them. A `FutureWarning`
        is emitted if any is given; they will be removed.

    Notes
    -----
    After `run`, the following diagnostics are available: ``max_sep`` (the
    kernel half width actually used, in pixels), ``mean_global`` (the clipped
    median subtracted before measuring the correlation function), ``xi`` and
    ``xi_clean`` (the measured and cleaned correlation functions, shape
    ``(2 * max_sep, 2 * max_sep)``, zero lag at pixel ``max_sep``), and
    ``n_iterations`` (the number of conjugate gradient iterations of each
    defect area).
    """

    def __init__(
        self,
        masked_image,
        defects=["SAT"],
        fwhm=5,
        kernel_half_width=30,
        fwhm_factor=3,
        power_threshold=2.5,
        apod_window="hann",
        apod_radius=None,
        white_noise=0.0,
        cg_rtol=1e-4,
        cg_maxiter=200,
        log=None,
        bin_spacing=None,
        threshold_dynamic_binning=None,
        threshold_subdivide=None,
    ):
        self.log = log or logging.getLogger(__name__)

        legacy = {
            "bin_spacing": bin_spacing,
            "threshold_dynamic_binning": threshold_dynamic_binning,
            "threshold_subdivide": threshold_subdivide,
        }
        given = sorted(name for name, value in legacy.items() if value is not None)
        if given:
            warnings.warn(
                f"InterpolateOverDefectGaussianProcess ignores the deprecated argument(s) {given}: "
                "the empirical-kernel Gaussian Process needs neither pixel binning nor prediction "
                "chunking. They will be removed.",
                FutureWarning,
                stacklevel=2,
            )

        self.masked_image = masked_image
        self.defects = defects
        self.fwhm = fwhm
        self.kernel_half_width = kernel_half_width
        self.fwhm_factor = fwhm_factor
        self.power_threshold = power_threshold
        self.apod_window = apod_window
        self.apod_radius = apod_radius
        self.white_noise = white_noise
        self.cg_rtol = cg_rtol
        self.cg_maxiter = cg_maxiter

        # Half width of the kernel grid (and support), in pixels. This is also
        # the distance around the defects within which good pixels are used
        # for the interpolation: pixels farther away have zero covariance with
        # every bad pixel.
        self.max_sep = int(np.ceil(max(kernel_half_width, fwhm_factor * fwhm)))
        if self.max_sep < 2:
            raise ValueError(
                "The kernel half width must span at least 2 pixels. "
                f"Current value: {self.max_sep} (kernel_half_width={kernel_half_width}, "
                f"fwhm_factor={fwhm_factor}, fwhm={fwhm})."
            )

        self.interpBit = self.masked_image.mask.getPlaneBitMask("INTRP")

        # Diagnostics, filled by run().
        self.mean_global = None
        self.xi = None
        self.xi_clean = None
        self.n_iterations = []
        self._engine = None

    def run(self):
        """
        Interpolate over the defects in the image.

        Change self.masked_image .
        """
        if self.defects == [] or self.defects is None:
            self.log.info("No defects found. No interpolation performed.")
            return

        mask = self.masked_image.getMask()
        bad_pixel_mask = mask.getPlaneBitMask(self.defects)
        if not np.any(mask.array & bad_pixel_mask):
            self.log.info("No bad pixels found. No interpolation performed.")
            return

        # Kernel measured once on the whole image, shared by all defects.
        engine = self._build_kernel(bad_pixel_mask)

        bad_mask_span_set = SpanSet.fromMask(mask, bad_pixel_mask).split()
        global_bbox = self.masked_image.getBBox()

        for spanset in bad_mask_span_set:
            # Dilate the bbox by the kernel half width to include all the good
            # pixels that have a non-zero covariance with the defect.
            localBox = spanset.getBBox().dilatedBy(self.max_sep)
            localBox.clip(global_bbox)
            masked_sub_image = self.masked_image[localBox]

            masked_sub_image = self.interpolate_masked_sub_image(masked_sub_image, engine)
            self.masked_image[localBox] = masked_sub_image

    def _build_kernel(self, bad_pixel_mask):
        """Measure and clean the 2-point correlation function of the good
        pixels of the whole image, and build the solve/predict engine.

        Parameters
        ----------
        bad_pixel_mask : `int`
            Bit mask of the defect planes.

        Returns
        -------
        engine : `treegp.GridConvolutionGP` or `None`
            The Gaussian Process engine holding the cleaned kernel, or `None`
            if no significant correlation was found (the defects are then
            filled with the local clipped median).
        """
        image = self.masked_image.image.array
        good = np.isfinite(image) & ((self.masked_image.mask.array & bad_pixel_mask) == 0)
        n_good = np.count_nonzero(good)
        if n_good == 0:
            self.log.warning("No good pixel in the image: the defects are filled with the local median.")
            return None

        self.mean_global = median_with_mad_clipping(image[good])
        self.xi = measure_2pcf_grid(image - self.mean_global, good, self.max_sep)

        # The treegp cleaning (apodization, Fourier thresholding) lives on the
        # empirical_2pcf solver. Its X, y, y_err arguments are only used by
        # the treecorr measurement, which is not done here (measure_2pcf_grid
        # is the exact equivalent on gridded data), so a minimal 2D field is
        # enough to set the kernel grid geometry.
        ny, nx = image.shape
        dummy_coords = np.array([[0.0, 0.0], [float(nx), float(ny)]])
        dummy_values = np.zeros(2)
        solver = treegp.empirical_2pcf(
            dummy_coords,
            dummy_values,
            dummy_values,
            max_sep=float(self.max_sep),
            pixel_size=1.0,
            power_threshold=self.power_threshold,
            apodize=True,
            apod_window=self.apod_window,
            apod_radius=self.apod_radius,
            apod_anisotropy=None,
        )
        if solver.npix != self.xi.shape[0]:
            raise RuntimeError(
                f"Inconsistent kernel grid: treegp expects {solver.npix} pixels, "
                f"the measured correlation function has {self.xi.shape[0]}."
            )
        try:
            self.xi_clean = solver.clean(self.xi)
        except RuntimeError as e:
            self.log.warning(
                "No significant correlation found in the image (%s): "
                "the defects are filled with the local median.",
                e,
            )
            self.xi_clean = None
            return None

        self.log.debug(
            "Empirical kernel measured on %d good pixels: half width %d pixels, "
            "variance %.3g (%.3g after cleaning).",
            n_good,
            self.max_sep,
            self.xi[self.max_sep, self.max_sep],
            self.xi_clean[self.max_sep, self.max_sep],
        )
        self._engine = treegp.GridConvolutionGP(
            self.xi_clean,
            pixel_size=1.0,
            upsample=1,
            cg_rtol=self.cg_rtol,
            cg_maxiter=self.cg_maxiter,
        )
        return self._engine

    def _pixel_errors(self, masked_sub_image, good_pixel):
        """Return the error of the given good pixels from the variance plane.

        Non-finite or non-positive variances are replaced by the median of
        the valid ones, and ``white_noise`` is added in quadrature.

        Parameters
        ----------
        masked_sub_image : `lsst.afw.image.MaskedImage`
            The sub-image the pixels belong to.
        good_pixel : `numpy.ndarray`
            Array of shape (n, 3) with the x and y (parent) coordinates of the
            pixels in its first two columns.

        Returns
        -------
        y_err : `numpy.ndarray`
            Standard deviation of each pixel, shape (n,), strictly positive.
        """
        x = np.rint(good_pixel[:, 0]).astype(int) - masked_sub_image.getX0()
        y = np.rint(good_pixel[:, 1]).astype(int) - masked_sub_image.getY0()
        variance = masked_sub_image.variance.array[y, x].astype(float)
        valid = np.isfinite(variance) & (variance > 0)
        if not np.all(valid):
            if np.any(valid):
                fill = np.median(variance[valid])
            else:
                fill = 1.0
                self.log.debug("No valid variance around the defect: unit pixel errors are used.")
            variance = np.where(valid, variance, fill)
        variance = variance + self.white_noise**2
        return np.sqrt(variance)

    def interpolate_masked_sub_image(self, masked_sub_image, engine):
        """
        Interpolate the masked sub-image.

        Parameters
        ----------
        masked_sub_image : `lsst.afw.image.MaskedImage`
            The sub-masked image to be interpolated.
        engine : `treegp.GridConvolutionGP` or `None`
            The Gaussian Process engine built by `_build_kernel`. If `None`,
            the bad pixels are filled with the local clipped median.

        Returns
        -------
        masked_sub_image : `lsst.afw.image.MaskedImage`
            The interpolated sub-masked image.
        """
        bad_pixel, good_pixel = ctUtils.findGoodPixelsAroundBadPixels(
            masked_sub_image, self.defects, buffer=self.max_sep
        )
        # Do nothing if there is nothing to interpolate or to interpolate from.
        if bad_pixel.size == 0 or good_pixel.size == 0:
            self.log.info("No bad or good pixels found. No interpolation performed.")
            return masked_sub_image

        finite = np.isfinite(good_pixel[:, 2])
        if not np.all(finite):
            good_pixel = good_pixel[finite]
            if good_pixel.size == 0:
                self.log.info("No finite good pixels found. No interpolation performed.")
                return masked_sub_image

        # Local mean: sky level around the defect.
        local_mean = median_with_mad_clipping(good_pixel[:, 2])

        if engine is None:
            bad_pixel[:, 2] = local_mean
        else:
            y_err = self._pixel_errors(masked_sub_image, good_pixel)
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                engine.solve(good_pixel[:, :2], good_pixel[:, 2] - local_mean, y_err)
            for w in caught:
                self.log.warning("Gaussian Process solve around %s: %s", bad_pixel[0, :2], w.message)
            self.n_iterations.append(engine.n_iterations)
            bad_pixel[:, 2] = engine.predict(bad_pixel[:, :2]) + local_mean

        # Update values
        ctUtils.updateImageFromArray(masked_sub_image.image, bad_pixel)
        updateMaskFromArray(masked_sub_image.mask, bad_pixel, self.interpBit)
        return masked_sub_image
