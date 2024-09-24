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


class InterpolateOverDefectGaussianProcessTestCase(lsst.utils.tests.TestCase):
    """Test InterpolateOverDefectGaussianProcess."""

    def setUp(self):
        super().setUp()

        npoints = 1000
        self.std = 100
        self.correlation_length = 10.0
        self.white_noise = 1e-5

        x1 = np.random.uniform(0, 99, npoints)
        x2 = np.random.uniform(0, 120, npoints)
        coord1 = np.array([x1, x2]).T

        kernel = rbf_kernel(coord1, coord1, self.std, self.correlation_length)
        kernel += np.eye(npoints) * self.white_noise**2

        # Data augmentation. Create a gaussian random field
        # on a 100 * 100 is to slow. So generate 1e3 points
        # and then interpolate it with a GP to do data augmentation.

        np.random.seed(42)
        z1 = np.random.multivariate_normal(np.zeros(npoints), kernel)

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
        # np.random.seed(12345)
        # self.noise.image.array[:, :] = np.random.normal(size=self.noise.image.array.shape)

    @lsst.utils.tests.methodParameters(method=("jax"))
    def test_interpolation(self, method: str):
        """Test that the interpolation is done correctly.

        Parameters
        ----------
        method : `str`
            Code used to solve gaussian process.
        """

        gp = InterpolateOverDefectGaussianProcess(
            self.maskedimage,
            defects=["BAD", "SAT", "CR", "EDGE"],
            method=method,
            fwhm=self.correlation_length,
            bin_image=False,
            bin_spacing=30,
            threshold_dynamic_binning=1000,
            threshold_subdivide=20000,
            correlation_length_cut=5,
            log=None,
        )

        gp.run()

        # Assert that the mask and the variance planes remain unchanged.
        self.assertImagesEqual(self.maskedimage.variance, self.reference.variance)

        # Check that interpolated pixels are close to the reference (original),
        # and that none of them is still NaN.
        self.assertTrue(np.isfinite(self.maskedimage.image.array).all())
        self.assertImagesAlmostEqual(
            self.maskedimage.image[1:, :],
            self.reference.image[1:, :],
            atol=2,
        )


def setup_module(module):
    lsst.utils.tests.init()


class MemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
