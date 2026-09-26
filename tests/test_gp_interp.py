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


import unittest

import numpy as np

import lsst.utils.tests
import lsst.geom
import lsst.afw.image as afwImage
from lsst.meas.algorithms import (
    InterpolateOverDefectGaussianProcess,
    GaussianProcessTreegp,
    measure_2pcf_grid,
)


def rbf_kernel(x1, x2, sigma, correlation_length):
    """
    Computes the radial basis function (RBF) kernel matrix.

    Parameters:
    -----------
    x1 : `np.array`
        Location of training data point with shape (n_samples, n_features).
    x2 : `np.array`
        Location of training/test data point with shape (n_samples, n_features).
    sigma : `float`
        The scale parameter of the kernel.
    correlation_length : `float`
        The correlation length parameter of the kernel.

    Returns:
    --------
    kernel : `np.array`
        RBF kernel matrix with shape (n_samples, n_samples).
    """
    distance_squared = np.sum((x1[:, None, :] - x2[None, :, :]) ** 2, axis=-1)
    kernel = (sigma**2) * np.exp(-0.5 * distance_squared / (correlation_length**2))
    return kernel


def brute_force_2pcf(image, good, max_sep):
    """Pair-count 2-point correlation function of a gridded image with gaps,
    by explicit loops over the lags (reference for `measure_2pcf_grid`).

    Parameters
    ----------
    image : `np.array`
        Mean-subtracted image, shape (ny, nx).
    good : `np.array`
        Boolean mask of the pixels to use, shape (ny, nx).
    max_sep : `int`
        Half width of the lag grid, in pixels.

    Returns
    -------
    xi : `np.array`
        Correlation function on a (2 * max_sep, 2 * max_sep) grid of lags,
        indexed [iy, ix], zero lag at pixel max_sep.
    """
    ny, nx = image.shape
    z = np.where(good, image, 0.0)
    xi = np.zeros((2 * max_sep, 2 * max_sep))
    for iy, dy in enumerate(range(-max_sep, max_sep)):
        for ix, dx in enumerate(range(-max_sep, max_sep)):
            num = 0.0
            den = 0
            for y in range(ny):
                y2 = y + dy
                if y2 < 0 or y2 >= ny:
                    continue
                for x in range(nx):
                    x2 = x + dx
                    if x2 < 0 or x2 >= nx:
                        continue
                    if good[y, x] and good[y2, x2]:
                        num += z[y, x] * z[y2, x2]
                        den += 1
            if den > 0:
                xi[iy, ix] = num / den
    return xi


class Measure2pcfGridTestCase(lsst.utils.tests.TestCase):
    """Test the FFT 2-point correlation function estimator."""

    def test_against_brute_force(self):
        rng = np.random.Generator(np.random.MT19937(7))
        ny, nx, max_sep = 13, 11, 3
        image = rng.normal(size=(ny, nx))
        # Correlate neighbouring pixels a bit so that xi is not a delta.
        image[:, 1:] += 0.5 * image[:, :-1]
        image[1:, :] += 0.3 * image[:-1, :]
        good = rng.uniform(size=(ny, nx)) > 0.3
        image[~good] = np.nan
        image_mean_sub = image - np.mean(image[good])

        xi = measure_2pcf_grid(image_mean_sub, good, max_sep)
        expected = brute_force_2pcf(image_mean_sub, good, max_sep)

        self.assertEqual(xi.shape, (2 * max_sep, 2 * max_sep))
        self.assertFloatsAlmostEqual(xi, expected, atol=1e-10, rtol=1e-10)
        # Zero lag is the variance of the good pixels, at pixel max_sep.
        self.assertFloatsAlmostEqual(
            xi[max_sep, max_sep], np.mean(image_mean_sub[good] ** 2), rtol=1e-10
        )
        # Point symmetry xi(-d) = xi(d) where both lags are on the grid.
        self.assertFloatsAlmostEqual(xi[1:, 1:], xi[1:, 1:][::-1, ::-1], atol=1e-10)

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            measure_2pcf_grid(np.zeros((4, 4)), np.ones((4, 5), dtype=bool), 2)
        with self.assertRaises(ValueError):
            measure_2pcf_grid(np.zeros((4, 4)), np.ones((4, 4), dtype=bool), 0)


class InterpolateOverDefectGaussianProcessTestCase(lsst.utils.tests.TestCase):
    """Test InterpolateOverDefectGaussianProcess."""

    def setUp(self):
        super().setUp()

        npoints = 1000
        self.std = 100
        self.correlation_length = 10.0
        self.white_noise = 1e-5

        rng = np.random.Generator(np.random.MT19937(5))

        x1 = rng.uniform(0, 99, npoints)
        x2 = rng.uniform(0, 120, npoints)
        coord1 = np.array([x1, x2]).T

        kernel = rbf_kernel(coord1, coord1, self.std, self.correlation_length)
        kernel += np.eye(npoints) * self.white_noise**2

        # Data augmentation. Create a gaussian random field
        # on a 100 * 100 is to slow. So generate 1e3 points
        # and then interpolate it with a GP to do data augmentation.

        z1 = rng.multivariate_normal(np.zeros(npoints), kernel)

        x1 = np.linspace(0, 99, 100)
        x2 = np.linspace(0, 120, 121)
        x2, x1 = np.meshgrid(x2, x1)
        coord2 = np.array([x1.reshape(-1), x2.reshape(-1)]).T

        tgp = GaussianProcessTreegp(
            std=self.std,
            correlation_length=self.correlation_length,
            white_noise=self.white_noise,
            mean=0.0,
        )
        tgp.fit(coord1, z1)
        z2 = tgp.predict(coord2)
        z2 = z2.reshape(100, 121)

        self.maskedimage = afwImage.MaskedImageF(100, 121)
        for x in range(100):
            for y in range(121):
                self.maskedimage[x, y] = (z2[x, y], 0, 1.0)

        # Clone the maskedimage so we can compare it after running the task.
        self.reference = self.maskedimage.clone()

        # Set some central pixels as SAT
        sliceX, sliceY = slice(30, 35), slice(40, 45)
        self.maskedimage.mask[sliceX, sliceY] = afwImage.Mask.getPlaneBitMask("SAT")
        self.maskedimage.image[sliceX, sliceY] = np.nan
        # Put nans here to make sure interp is done ok

        # Set an entire column as BAD
        self.maskedimage.mask[54:55, :] = afwImage.Mask.getPlaneBitMask("BAD")
        self.maskedimage.image[54:55, :] = np.nan

        # Set an entire row as BAD
        self.maskedimage.mask[:, 110:111] = afwImage.Mask.getPlaneBitMask("BAD")
        self.maskedimage.image[:, 110:111] = np.nan

        # Set a diagonal set of pixels as CR
        for i in range(74, 78):
            self.maskedimage.mask[i, i] = afwImage.Mask.getPlaneBitMask("CR")
            self.maskedimage.image[i, i] = np.nan

        # Set one of the edges as EDGE
        self.maskedimage.mask[0:1, :] = afwImage.Mask.getPlaneBitMask("EDGE")
        self.maskedimage.image[0:1, :] = np.nan

        # Set a smaller streak at the edge
        self.maskedimage.mask[25:28, 0:1] = afwImage.Mask.getPlaneBitMask("EDGE")
        self.maskedimage.image[25:28, 0:1] = np.nan

        # Update the reference image's mask alone, so we can compare them after
        # running the task.
        self.reference.mask.array[:, :] = self.maskedimage.mask.array

        # Create a noise image
        # self.noise = self.maskedimage.clone()
        # rng = np.random.Generator(np.random.MT19937(5))
        # self.noise.image.array[:, :] = rng.normal(size=self.noise.image.array.shape)

    def test_interpolation(self):
        """Test that the interpolation is done correctly, with both
        correlation function estimators."""
        for two_pcf_method in ("treecorr", "fft"):
            with self.subTest(two_pcf_method=two_pcf_method):
                maskedimage = self.maskedimage.clone()
                gp = InterpolateOverDefectGaussianProcess(
                    maskedimage,
                    defects=["BAD", "SAT", "CR", "EDGE"],
                    fwhm=self.correlation_length,
                    kernel_half_width=20,
                    fwhm_factor=3,
                    two_pcf_method=two_pcf_method,
                    log=None,
                )
                # The kernel half width is at least fwhm_factor * fwhm.
                self.assertEqual(gp.max_sep, 30)

                gp.run()

                # One Gaussian Process per connected defect, each either
                # solved or filled with the local median.
                self.assertGreater(gp.n_areas, 0)
                self.assertEqual(len(gp.n_iterations) + gp.n_fallback, gp.n_areas)
                if gp.n_fallback < gp.n_areas:
                    pixel_size = gp.two_pcf_pixel_size if two_pcf_method == "treecorr" else 1
                    npix = 2 * int(np.ceil(gp.max_sep / pixel_size))
                    self.assertEqual(gp.last_xi.shape, (npix, npix))
                    self.assertEqual(gp.last_xi_clean.shape, (npix, npix))
                    self.assertGreater(gp.last_xi_clean[npix // 2, npix // 2], 0.0)

                # Assert that the mask and the variance planes remain unchanged.
                self.assertImagesEqual(maskedimage.variance, self.reference.variance)

                # The interpolated pixels are flagged as such, and only them.
                interpBit = maskedimage.mask.getPlaneBitMask("INTRP")
                badBits = maskedimage.mask.getPlaneBitMask(["BAD", "SAT", "CR", "EDGE"])
                isInterp = (maskedimage.mask.array & interpBit) != 0
                isBad = (self.reference.mask.array & badBits) != 0
                np.testing.assert_array_equal(isInterp, isBad)

                # Check that interpolated pixels are close to the reference
                # (original), and that none of them is still NaN.
                self.assertTrue(np.isfinite(maskedimage.image.array).all())
                self.assertImagesAlmostEqual(
                    maskedimage.image[1:, :],
                    self.reference.image[1:, :],
                    atol=5,
                )

    def test_legacy_kwargs(self):
        """The kwargs of the former dense solver, still passed by ip_isr,
        are accepted, ignored, and warned about."""
        with self.assertWarns(FutureWarning):
            gp = InterpolateOverDefectGaussianProcess(
                self.maskedimage,
                defects=["SAT"],
                fwhm=15,
                bin_spacing=20,
                threshold_dynamic_binning=2000,
                threshold_subdivide=20000,
            )
        # The ignored arguments do not change the configuration.
        self.assertEqual(gp.max_sep, 45)
        self.assertFalse(hasattr(gp, "bin_spacing"))


def setup_module(module):
    lsst.utils.tests.init()


class MemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
