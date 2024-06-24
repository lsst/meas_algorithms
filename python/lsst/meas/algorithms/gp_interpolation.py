import numpy as np
from lsst.meas.algorithms import CloughTocher2DInterpolatorUtils as ctUtils
from lsst.geom import Box2I, Point2I, Extent2I
from lsst.afw.geom import SpanSet
import copy
import treegp

import jax
from jax import jit
import jax.numpy as jnp

import warnings

__all__ = ["interpolateOverDefectsGP"]


def updateMaskFromArray(mask, bad_pixel, interpBit):
    """
    Update the mask array with the given bad pixels.

    Parameters
    ----------
    mask : MaskedImage #TO DO, this is not the right format need to check.
        The mask image to update.
    bad_pixel : array-like
        An array-like object containing the coordinates of the bad pixels.
        Each row should contain the x and y coordinates of a bad pixel.
    interpBit : int
        The bit value to set for the bad pixels in the mask.

    Returns
    -------
    None
    """
    x0 = mask.getX0()
    y0 = mask.getY0()
    for row in bad_pixel:
        x = int(row[0] - x0)
        y = int(row[1] - y0)
        mask.array[y, x] |= interpBit


@jit
def median_with_mad_clipping(data, mad_multiplier=2.0):
    """
    Calculate the median of the input data after applying Median Absolute Deviation (MAD) clipping.

    Parameters:
    -----------
    data : array-like
        Input data array.
    mad_multiplier : float, optional
        Multiplier for the MAD value used for clipping. Default is 2.0.

    Returns:
    --------
    median_clipped : float
        Median value of the clipped data.

    Notes:
    ------
    The MAD clipping method is used to remove outliers from the data. The median of the data is calculated,
    and then the MAD is calculated as the median absolute deviation from the median. The data is then clipped
    by removing values that are outside the range of median +/- mad_multiplier * MAD. Finally, the median of
    the clipped data is returned.

    Examples:
    ---------
    >>> data = [1, 2, 3, 4, 5, 100]
    >>> median_with_mad_clipping(data)
    3.5
    """
    median = jnp.median(data)
    mad = jnp.median(jnp.abs(data - median))
    clipping_range = mad_multiplier * mad
    clipped_data = jnp.clip(data, median - clipping_range, median + clipping_range)
    median_clipped = jnp.median(clipped_data)
    return median_clipped


# Below are the jax functions for Gaussian Process regression.
# Kernel is assumed to be isotropic RBF kernel, and solved
# using exact Cholesky decomposition.
# The interpolation solution is obtained by solving the linear system:
# y_interp = kernel_rect @ (kernel + y_err**2 * I)^-1 @ y_training.
# See the Rasmussen and Williams book for more details.
# Each function is decorated with @jit to compile the function.
# Exist package like tinygp, that is implemented in jax also, but
# wanted to control the implementation here, as I am doing a lot
# of custom things (setting my own hyperparameters, fine tune mean function,
# dynamic binning, ...).


@jit
def jax_pdist_squared(x):
    """
    Calculate the squared pairwise Euclidean distances between points in x.

    Parameters:
    ----------
    x : array-like
        Input array of shape (n_samples, n_features).

    Returns:
    -------
    dist : ndarray
        Array of shape (n_samples, n_samples) containing the squared pairwise Euclidean distances.

    Notes:
    ------
    The pairwise Euclidean distance between two points `p` and `q` is defined as:
    dist(p, q) = sqrt(sum((p - q)^2))

    Examples:
    ---------
    >>> x = np.array([[0, 0], [1, 1], [2, 2]])
    >>> jax_pdist_squared(x)
    array([[0, 2, 8],
           [2, 0, 2],
           [8, 2, 0]])
    """
    return jnp.sum((x[:, None, :] - x[None, :, :]) ** 2, axis=-1)


@jit
def jax_cdist_squared(xa, xb):
    """
    Compute the squared Euclidean distance between two sets of points.

    Parameters
    ----------
    xa : array-like, shape (n_samples_a, n_features)
        The first set of points.
    xb : array-like, shape (n_samples_b, n_features)
        The second set of points.

    Returns
    -------
    dist : ndarray, shape (n_samples_a, n_samples_b)
        The squared Euclidean distance between each pair of points.

    Notes
    -----
    This function uses the broadcasting feature of NumPy to efficiently compute
    the squared Euclidean distance between two sets of points.

    Examples
    --------
    >>> xa = np.array([[0, 0], [1, 1], [2, 2]])
    >>> xb = np.array([[1, 0], [0, 1]])
    >>> jax_cdist_squared(xa, xb)
    array([[1, 2],
           [1, 2],
           [5, 8]])

    """
    return jnp.sum((xa[:, None, :] - xb[None, :, :]) ** 2, axis=-1)


@jit
def jax_rbf_kernel(x, sigma, correlation_length, y_err):
    """
    Computes the radial basis function (RBF) kernel matrix.

    Parameters:
    -----------
    x : array-like
        Input data points with shape (n_samples, n_features).
    sigma : float
        The scale parameter of the kernel.
    correlation_length : float
        The correlation length parameter of the kernel.
    y_err : float
        Measurement error for the input values.

    Returns:
    --------
    kernel : array-like
        RBF kernel matrix with shape (n_samples, n_samples).
    """
    distance_squared = jax_pdist_squared(x)
    kernel = (sigma**2) * jnp.exp(-0.5 * distance_squared / (correlation_length**2))
    y_err = jnp.ones(len(x[:, 0])) * y_err
    kernel += jnp.eye(len(y_err)) * (y_err**2)
    return kernel


@jit
def jax_rbf_kernel_rect(x1, x2, sigma, correlation_length):
    """
    Compute the radial basis function (RBF) kernel (rectangular matrix).
    Parameters:
    -----------
    x1 : array-like
        The first set of input points.
    x2 : array-like
        The second set of input points.
    sigma : float
        The scale parameter of the kernel.
    correlation_length : float
        The correlation length parameter of the kernel.

    Returns:
    --------
    kernel_rect : ndarray
        The computed RBF kernel (rectangular matrix).

    """
    l1 = jax_cdist_squared(x1, x2)
    kernel_rect = (sigma**2) * jnp.exp(-0.5 * l1 / (correlation_length**2))
    return kernel_rect


@jit
def jax_get_alpha(y, kernel):
    """
    Compute the alpha vector for Gaussian Process interpolation.

    Parameters:
    -----------
    y : array_like
        The target values of the Gaussian Process.
    kernel : array_like
        The kernel matrix of the Gaussian Process.

    Returns:
    --------
    alpha : ndarray
        The alpha vector computed using the Cholesky decomposition and solve.

    """
    factor = (jax.scipy.linalg.cholesky(kernel, overwrite_a=True, lower=False), False)
    alpha = jax.scipy.linalg.cho_solve(factor, y, overwrite_b=False)
    return alpha.reshape((len(alpha), 1))


@jit
def jax_get_y_predict(kernel_rect, alpha):
    """
    Compute the predicted values of y using the given kernel and alpha.

    Parameters:
    -----------
    kernel_rect : array-like
        The kernel matrix.
    alpha : array-like
        The alpha vector.

    Returns:
    --------
    array-like
        The predicted values of y.

    """
    return jnp.dot(kernel_rect, alpha).T[0]


class GaussianProcessJax:
    def __init__(self, std=1.0, correlation_length=1.0, white_noise=0.0, mean=0.0):
        self.std = std
        self.correlation_lenght = correlation_length
        self.white_noise = white_noise
        self.mean = mean
        self._alpha = None

    def fit(self, x_good, y_good):
        y = y_good - self.mean
        self._x = x_good
        kernel = jax_rbf_kernel(x_good, self.std, self.correlation_lenght, self.white_noise)
        self._alpha = jax_get_alpha(y, kernel)

    def predict(self, x_bad):
        kernel_rect = jax_rbf_kernel_rect(x_bad, self._x, self.std, self.correlation_lenght)
        y_pred = jax_get_y_predict(kernel_rect, self._alpha)
        return y_pred + self.mean


# Vanilla Gaussian Process regression using treegp package
# There is no fancy O(N*log(N)) solver here, just the basic GP regression (Cholesky).
class GaussianProcessTreegp:
    """
    Gaussian Process Treegp class for Gaussian Process interpolation.

    Parameters:
    -----------
    std : float, optional
        Standard deviation of the Gaussian Process kernel. Default is 1.0.
    correlation_length : float, optional
        Correlation length of the Gaussian Process kernel. Default is 1.0.
    white_noise : float, optional
        White noise level of the Gaussian Process. Default is 0.0.
    mean : float, optional
        Mean value of the Gaussian Process. Default is 0.0.

    Methods:
    --------
    fit(x_good, y_good):
        Fit the Gaussian Process to the given training data.

        Parameters:
        -----------
        x_good : array-like
            Input features for the training data.
        y_good : array-like
            Target values for the training data.

    predict(x_bad):
        Predict the target values for the given input features.

        Parameters:
        -----------
        x_bad : array-like
            Input features for the prediction.

        Returns:
        --------
        y_pred : array-like
            Predicted target values.

    """

    def __init__(self, std=1.0, correlation_length=1.0, white_noise=0.0, mean=0.0):
        self.std = std
        self.correlation_length = correlation_length
        self.white_noise = white_noise
        self.mean = mean

    def fit(self, x_good, y_good):
        """
        Fit the Gaussian Process to the given training data.

        Parameters:
        -----------
        x_good : array-like
            Input features for the training data.
        y_good : array-like
            Target values for the training data.
        """
        kernel = f"{self.std}**2 * RBF({self.correlation_length})"
        self.gp = treegp.GPInterpolation(
            kernel=kernel,
            optimizer="none",
            normalize=False,
            white_noise=self.white_noise,
        )
        self.gp.initialize(x_good, y_good - self.mean)
        self.gp.solve()

    def predict(self, x_bad):
        """
        Predict the target values for the given input features.

        Parameters:
        -----------
        x_bad : array-like
            Input features for the prediction.

        Returns:
        --------
        y_pred : array-like
            Predicted target values.
        """
        y_pred = self.gp.predict(x_bad)
        return y_pred + self.mean


class InterpolateOverDefectGaussianProcess:
    """
    InterpolateOverDefectGaussianProcess class performs Gaussian Process
    (GP) interpolation over defects in an image.

    Parameters:
    -----------
    maskedImage : `lsst.afw.image.MaskedImage`
        The masked image containing the defects to be interpolated.
    defects : list of str, optional
        The types of defects to be interpolated. Default is ["SAT"].
    method : str, optional
        The method to use for GP interpolation. Must be either "jax" or "treegp". Default is "treegp".
    fwhm : float, optional
        The full width at half maximum (FWHM) of the PSF. Default is 5.
    bin_spacing : int, optional
        The spacing between bins for good pixel binning. Default is 10.
    threshold_dynamic_binning : int, optional
        The threshold for dynamic binning. Default is 1000.
    threshold_subdivide : int, optional
        The threshold for sub-dividing the bad pixel array to avoid memory error. Default is 20000.
    correlation_length_cut : int, optional
        The factor by which to dilate the bounding box around defects. Default is 5.

    Raises:
    -------
    ValueError
        If an invalid method is provided.

    Attributes:
    -----------
    bin_spacing : int
        The spacing between bins for good pixel binning.
    threshold_subdivide : int
        The threshold for sub-dividing the bad pixel array to avoid memory error.
    threshold_dynamic_binning : int
        The threshold for dynamic binning.
    maskedImage : `lsst.afw.image.MaskedImage`
        The masked image containing the defects to be interpolated.
    defects : list of str
        The types of defects to be interpolated.
    correlation_length : float
        The correlation length (FWHM).
    correlation_length_cut : int
        The factor by which to dilate the bounding box around defects.
    interpBit : int
        The bit mask for the "INTRP" plane in the image mask.

    Methods:
    --------
    interpolate_over_defects()
        Interpolates over the defects in the image.
    _good_pixel_binning(good_pixel)
        Performs good pixel binning.
    interpolate_sub_masked_image(sub_masked_image)
        Interpolates the sub-masked image.

    """

    def __init__(
        self,
        maskedImage,
        defects=["SAT"],
        method="treegp",
        fwhm=5,
        bin_spacing=10,
        threshold_dynamic_binning=1000,
        threshold_subdivide=20000,
        correlation_length_cut=5,
    ):
        if method == "jax":
            self.GaussianProcess = GaussianProcessJax
        elif method == "treegp":
            self.GaussianProcess = GaussianProcessTreegp
        else:
            raise ValueError("Invalid method. Must be 'jax' or 'treegp'.")

        self.bin_spacing = bin_spacing
        self.threshold_subdivide = threshold_subdivide
        self.threshold_dynamic_binning = threshold_dynamic_binning

        self.maskedImage = maskedImage
        self.defects = defects
        self.correlation_length = fwhm
        self.correlation_length_cut = correlation_length_cut

        self.interpBit = self.maskedImage.mask.getPlaneBitMask("INTRP")

    def interpolate_over_defects(self):
        """
        Interpolates over the defects in the image.
        """

        mask = self.maskedImage.getMask()
        badPixelMask = mask.getPlaneBitMask(self.defects)
        badMaskSpanSet = SpanSet.fromMask(mask, badPixelMask).split()

        bbox = self.maskedImage.getBBox()
        glob_xmin, glob_xmax = bbox.minX, bbox.maxX
        glob_ymin, glob_ymax = bbox.minY, bbox.maxY

        for spanset in badMaskSpanSet:
            bbox = spanset.getBBox()
            # Dilate the bbox to make sure we have enough good pixels around the defect
            # For now, we dilate by 5 times the correlation length
            # For GP with isotropic kernel, points at 5 correlation lengths away have negligible
            # effect on the prediction.
            bbox = bbox.dilatedBy(
                int(self.correlation_length * self.correlation_length_cut)
            )  # need integer as input.
            xmin, xmax = max([glob_xmin, bbox.minX]), min(glob_xmax, bbox.maxX)
            ymin, ymax = max([glob_ymin, bbox.minY]), min(glob_ymax, bbox.maxY)
            localBox = Box2I(Point2I(xmin, ymin), Extent2I(xmax - xmin, ymax - ymin))
            try:
                sub_masked_image = self.maskedImage[localBox]
            except IndexError:
                raise ValueError("Sub-masked image not found.")

            sub_masked_image = self.interpolate_sub_masked_image(sub_masked_image)
            self.maskedImage[localBox] = sub_masked_image

    def _good_pixel_binning(self, good_pixel):
        """
        Performs good pixel binning.

        Parameters:
        -----------
        good_pixel : numpy.ndarray
            The array of good pixels.

        Returns:
        --------
        numpy.ndarray
            The binned array of good pixels.
        """

        n_pixels = len(good_pixel[:, 0])
        dynamic_binning = int(np.sqrt(n_pixels / self.threshold_dynamic_binning))
        if n_pixels / self.bin_spacing**2 < n_pixels / dynamic_binning**2:
            bin_spacing = self.bin_spacing
        else:
            bin_spacing = dynamic_binning
        binning = treegp.meanify(bin_spacing=bin_spacing, statistics="mean")
        binning.add_field(
            good_pixel[:, :2],
            good_pixel[:, 2:].T,
        )
        binning.meanify()
        return np.array(
            [binning.coords0[:, 0], binning.coords0[:, 1], binning.params0]
        ).T

    def interpolate_sub_masked_image(self, sub_masked_image):
        """
        Interpolates the sub-masked image.

        Parameters:
        -----------
        sub_masked_image : `lsst.afw.image.MaskedImage`
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
            sub_masked_image, self.defects, buffer=cut
        )
        # Do nothing if bad pixel is None.
        if bad_pixel.size == 0 or good_pixel.size == 0:
            warnings.warn("No bad or good pixels found. No interpolation performed.")
            return sub_masked_image
        # Do GP interpolation if bad pixel found.
        else:
            # gp interpolation
            sub_image_array = sub_masked_image.getVariance().array
            white_noise = np.sqrt(
                np.mean(sub_image_array[np.isfinite(sub_image_array)])
            )
            kernel_amplitude = np.std(good_pixel[:, 2:])
            if not np.isfinite(kernel_amplitude):
                filter_finite = np.isfinite(good_pixel[:, 2:]).T[0]
                good_pixel = good_pixel[filter_finite]
                if good_pixel.size == 0:
                    warnings.warn(
                        "No bad or good pixels found. No interpolation performed."
                    )
                    return sub_masked_image
                # kernel amplitude might be better described by maximum value of good pixel given
                # the data and not really a random gaussian field.
                kernel_amplitude = np.max(good_pixel[:, 2:])
            else:
                kernel_amplitude = np.max(good_pixel[:, 2:])
            try:
                good_pixel = self._good_pixel_binning(copy.deepcopy(good_pixel))
            except Exception:
                warnings.warn(
                    "Binning failed, use original good pixel array in interpolate over."
                )

            # put this after binning as comupting median is O(n*log(n))
            mean = median_with_mad_clipping(good_pixel[:, 2:])

            gp = self.GaussianProcess(
                std=np.sqrt(kernel_amplitude),
                correlation_length=self.correlation_length,
                white_noise=white_noise,
                mean=mean,
            )
            gp.fit(good_pixel[:, :2], np.squeeze(good_pixel[:, 2:]))
            if bad_pixel.size < self.threshold_subdivide:
                gp_predict = gp.predict(bad_pixel[:, :2])
                bad_pixel[:, 2:] = gp_predict.reshape(np.shape(bad_pixel[:, 2:]))
            else:
                warnings.warn("sub-divide bad pixel array to avoid memory error.")
                for i in range(0, len(bad_pixel), self.threshold_subdivide):
                    end = min(i + self.threshold_subdivide, len(bad_pixel))
                    gp_predict = gp.predict(bad_pixel[i:end, :2])
                    bad_pixel[i:end, 2:] = gp_predict.reshape(
                        np.shape(bad_pixel[i:end, 2:])
                    )

            # update_value
            ctUtils.updateImageFromArray(sub_masked_image.image, bad_pixel)
            updateMaskFromArray(sub_masked_image.mask, bad_pixel, self.interpBit)
            return sub_masked_image


def interpolateOverDefectsGP(
    image,
    fwhm,
    badList,
    method="treegp",
    bin_spacing=25,
    threshold_dynamic_binning=1000,
    threshold_subdivide=20000,
):
    """
    Interpolates over defects in an image using Gaussian Process interpolation.

    Parameters
    ----------
    image : ndarray
        The input image.
    fwhm : float
        The full width at half maximum (FWHM) of the Gaussian kernel used for interpolation.
    badList : list
        A list of defect coordinates in the image.
    method : str, optional
        The method used for interpolation. Default is "treegp".
    bin_spacing : int, optional
        The spacing between bins used for dynamic binning. Default is 25.
    threshold_dynamic_binning : int, optional
        The threshold for dynamic binning. Default is 1000.
    threshold_subdivide : int, optional
        The threshold for subdividing defects. Default is 20000.

    Returns
    -------
    None

    Raises
    ------
    UserWarning
        If no defects are found in the image.

    Notes
    -----
    This function performs Gaussian Process interpolation over defects in the input image.
    It uses the provided defect coordinates to identify and interpolate over the defects.
    The interpolated image is not returned, instead, the input image is modified in-place.

    """
    if badList == [] or badList is None:
        warnings.warn("WARNING: no defects found. No interpolation performed.")
        return
    gp = InterpolateOverDefectGaussianProcess(
        image,
        defects=badList,
        method=method,
        fwhm=fwhm,
        bin_spacing=bin_spacing,
        threshold_dynamic_binning=threshold_dynamic_binning,
        threshold_subdivide=threshold_subdivide,
    )
    gp.interpolate_over_defects()
