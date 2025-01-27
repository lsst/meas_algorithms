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

import itertools
import unittest
import warnings

import lsst.utils.tests
import numpy as np
from lsst.afw import image as afwImage
from lsst.meas.algorithms import (
    ComputeNoiseCorrelationConfig,
    ComputeNoiseCorrelationTask,
    ScaleVarianceTask,
)


class NoiseVarianceTestCase(lsst.utils.tests.TestCase):
    @classmethod
    def setUpClass(cls, seed: int = 12345, size: int = 512) -> None:
        """
        Set up a common noise field and a masked image for all test cases.

        Parameters
        ----------
        size : int, optional
            Size of the noise image to generate.
        seed : int, optional
            Seed for the random number generator.
        """

        super().setUpClass()

        np.random.seed(seed)
        cls.noise = np.random.randn(size, size).astype(np.float32)

        # We will clip the edges, so the variance plane will be smaller by 2.
        variance_array = np.ones((size - 2, size - 2), dtype=np.float32)

        # Randomly set some mask bits to be non-zero.
        mask_array = (np.random.geometric(0.85, size=(size - 2, size - 2)) - 1).astype(
            np.int32
        )
        # Set some masked variance values to zero.
        variance_array[mask_array > 1] = 0.0

        cls.mi = afwImage.makeMaskedImageFromArrays(
            image=cls.noise[1:-1, 1:-1].copy(), variance=variance_array, mask=mask_array
        )

    @classmethod
    def tearDownClass(cls) -> None:
        del cls.mi
        del cls.noise
        super().tearDownClass()

    def _prepareImage(self, rho1: float, rho2: float, background_value: float = 0.0) -> None:
        """Create a correlated Gaussian noise field using simple translations.

        Y[i,j] = ( X[i,j] + a1 X[i-1,j] + a2 X[i,j-1] )/sqrt(1 + a1**2 + a2**2)

        Var( X[i,j] ) = Var( Y[i,j] ) = 1
        Cov( Y[i,j], V[i-1,j] ) = a1
        Cov( Y[i,j], V[i,j-1] ) = a2

        rho_i = a_i/sqrt(1 + a1**2 + a2**2) for i = 1, 2

        Parameters
        ----------
        rho1 : float
            Correlation coefficient along the horizontal (x) direction.
        rho2 : float
            Correlation coefficient along the vertical (y) direction.
        background_value : float, optional
            A constant background to add to the image.
        """
        # Solve for the kernel parameters (a1, a2) & generate correlated noise.
        r2 = rho1**2 + rho2**2
        if r2 > 0:
            k = 0.5 * (1 + np.sqrt(1 - 4 * r2)) / r2
            a1, a2 = k * rho1, k * rho2
            self.noise += background_value
            try:
                corr_noise = (
                    self.noise
                    + a1 * np.roll(self.noise, 1, axis=0)
                    + a2 * np.roll(self.noise, 1, axis=1)
                ) / np.sqrt(1 + a1**2 + a2**2)
            finally:
                self.noise -= background_value
        else:
            a1, a2 = 0, 0
            corr_noise = self.noise + background_value

        self.mi.image.array = corr_noise[1:-1, 1:-1].astype(np.float32)

    @lsst.utils.tests.methodParameters(
        rho=((0.0, 0.0), (-0.2, 0.0), (0.0, 0.1), (0.15, 0.25), (0.25, -0.15))
    )
    def testScaleVariance(self, rho):
        """Test that the ScaleVarianceTask scales the variance plane correctly."""
        task = ScaleVarianceTask()
        rho1, rho2 = rho
        self._prepareImage(rho1, rho2)
        scaleFactors = task.computeScaleFactors(self.mi)

        # Check for consistency between pixelFactor and imageFactor
        self.assertFloatsAlmostEqual(
            scaleFactors.pixelFactor, scaleFactors.imageFactor, atol=1e-6
        )

        # Since the variance is expected to remain unity after introducing the
        # correlations, the scaleFactor should be 1.0 within statistical error.
        self.assertFloatsAlmostEqual(scaleFactors.pixelFactor, 1.0, rtol=2e-2)

    @lsst.utils.tests.methodParametersProduct(
        rho=((0.0, 0.0), (0.2, 0.0), (0.0, -0.1), (0.15, 0.25), (-0.25, 0.15)),
        scaleEmpircalVariance=(False, True),
        subtractEmpiricalMean=(False, True),
        background_value=(0.0, 100.0),
    )
    def testComputeCorrelation(
        self, rho, background_value, scaleEmpircalVariance, subtractEmpiricalMean
    ):
        """Test that the noise correlation coefficients are computed correctly."""
        corr_matrix_size = 5
        config = ComputeNoiseCorrelationConfig(size=corr_matrix_size)
        config.scaleEmpiricalVariance = scaleEmpircalVariance
        config.subtractEmpiricalMean = subtractEmpiricalMean
        task = ComputeNoiseCorrelationTask(config=config)

        rho1, rho2 = rho
        self._prepareImage(rho1, rho2, background_value=background_value)

        # Create a copy of the image before running the task.
        mi_copy = afwImage.MaskedImage(self.mi, dtype=self.mi.dtype)
        self.assertIsNot(mi_copy, self.mi)  # Check that it's a deepcopy.

        with warnings.catch_warnings(record=True) as warning_list:
            corr_matrix = task.run(self.mi)

        # Check that when dividing the background pixels by per-pixel variance
        # we did not divide by zero accidentally.
        self.assertEqual(sum(w.category is RuntimeWarning for w in warning_list), 0)

        # Check that the task did not modify the input in place accidentally.
        self.assertIsNot(mi_copy, self.mi)
        np.testing.assert_array_equal(mi_copy.image.array, self.mi.image.array)
        np.testing.assert_array_equal(mi_copy.mask.array, self.mi.mask.array)
        np.testing.assert_array_equal(mi_copy.variance.array, self.mi.variance.array)

        # corr_matrix elements should be zero except for (1,0), (0,1) & (0,0).
        # Use the other elements to get an estimate of the statistical
        # uncertainty in our estimates.
        err = np.std(
            [
                corr_matrix(i, j)
                for i, j in itertools.product(
                    range(corr_matrix_size), range(corr_matrix_size)
                )
                if (i + j > 1)
            ]
        )

        self.assertLess(abs(corr_matrix(1, 0) / corr_matrix(0, 0) - rho1), 3 * err)
        self.assertLess(abs(corr_matrix(0, 1) / corr_matrix(0, 0) - rho2), 3 * err)
        self.assertLess(err, 3e-3)  # Check that the err is much less than rho.


class MemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
