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
from typing import Optional
import numpy.typing as npt

import numpy as np

import lsst.utils.tests
from lsst.meas.algorithms import GaussianProcessTreegp
from lsst.meas.algorithms.computeExPsf import ComputeExPsfTask
import lsst.pipe.base as pipeBase


def rbf_kernel(
    x1: npt.NDArray, x2: npt.NDArray, sigma: float, correlation_length: float
) -> npt.NDArray:
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
    kernel : `np.ndarray`
        RBF kernel matrix with shape (n_samples, n_samples).
    """
    distance_squared = np.sum((x1[:, None, :] - x2[None, :, :]) ** 2, axis=-1)
    kernel = (sigma**2) * np.exp(-0.5 * distance_squared / (correlation_length**2))
    return kernel


def generate_gaussian_random_field(
    xmin: int = 0,
    xmax: int = 2000,
    ymin: int = 0,
    ymax: int = 2000,
    npoints: int = 10000,
    nmax: int = 1000,
    std: float = 1.0,
    correlation_length: float = 10.0,
    white_noise: float = 1.0,
    seed: int = 42,
    input_coord: Optional[npt.NDArray] = None,
) -> pipeBase.Struct:
    """Generate a Gaussian Random Field.

    Function to generate a Gaussian Random Field.
    Help for unit test and generate spatial correlated
    function and have an analytical 2-point correlation
    to compared with (Gaussian here).

    Parameters
    ----------
    xmin: `int`
        Min value in x direction.
        Default: ``0``
    xmax: `int`
        Max value in x direction.
        Default: ``2000``
    ymin: `int`
        Min value in y direction.
        Default: ``0``
    ymax: `int`
        Max value in y direction.
        Default: ``2000``
    npoints: `int`
        Number of data points generated.
        Default: ``10000``
    nmax: `int`
        Max number of data points generated using
        `np.random.Generator.multivariate_normal`. If
        npoints>nmax, a GP will be used in addition.
        Default: ``1000``
    std: `float`
        Amplitude of the gaussian random field.
        Default: ``1.0``
    correlation_length: `float`
        Correlation length of the gaussian random field.
        Default: ``10.0``
    white_noise: `float`
        Noise added to the gaussian random field.
        Default: ``1.0``
    seed: `int`
        Seed of the random generator.
        Default: ``42``
    input_coord: `np.ndarray`
        Take a input coord to generate the Gaussian Random field
        Default: ``None``

    Returns
    -------
    struct : `lsst.pipe.base.Struct`
            The struct contains the following data:
        ``coord``: `np.ndarray`
            2D coordinate of the gaussian random field
        ``z``: `np.ndarray`
            Scalar value of the gaussian random field
    """
    # Choose explicit RNG.
    rng = np.random.Generator(np.random.MT19937(seed))

    if input_coord is not None:
        npoints = len(input_coord[:, 0])

    if npoints > nmax:
        ngenerated = nmax
    else:
        ngenerated = npoints

    if input_coord is None or npoints > nmax:
        x1 = rng.uniform(xmin, xmax, ngenerated)
        x2 = rng.uniform(ymin, ymax, ngenerated)
        coord1 = np.array([x1, x2]).T
    else:
        coord1 = input_coord
    kernel = rbf_kernel(coord1, coord1, std, correlation_length)
    kernel += np.eye(ngenerated) * white_noise**2

    z1 = rng.multivariate_normal(np.zeros(ngenerated), kernel)

    # Data augmentation. Create a gaussian random field
    # with npoints>nmax is to slow. So generate nmax points
    # and then interpolate it with a GP to do data augmentation.

    if npoints > nmax:
        if input_coord is None:
            x1 = rng.uniform(xmin, xmax, npoints)
            x2 = rng.uniform(ymin, ymax, npoints)
            coord = np.array([x1, x2]).T
        else:
            coord = input_coord

        tgp = GaussianProcessTreegp(
            std=std,
            correlation_length=correlation_length,
            white_noise=white_noise,
            mean=0.0,
        )
        tgp.fit(coord1, z1)
        z = tgp.predict(coord)
    else:
        coord = coord1
        z = z1

    return pipeBase.Struct(coord=coord, z=z)


class ComputeExPsfTestCase(lsst.utils.tests.TestCase):
    """Test ComputeExPsfTask."""

    def setUp(self) -> None:
        super().setUp()

        output1 = generate_gaussian_random_field(
            xmin=0,
            xmax=2000,
            ymin=0,
            ymax=2000,
            npoints=10000,
            nmax=2000,
            std=1.0,
            correlation_length=200.0,
            white_noise=0.01,
            seed=1,
            input_coord=None,
        )

        output2 = generate_gaussian_random_field(
            xmin=0,
            xmax=2000,
            ymin=0,
            ymax=2000,
            npoints=10000,
            nmax=2000,
            std=1.0,
            correlation_length=100.0,
            white_noise=0.01,
            seed=44,
            input_coord=output1.coord,
        )

        self.coord1 = output1.coord
        self.coord2 = output2.coord
        self.de1 = output1.z
        self.de2 = output2.z

    def test_comp_ex_psf(self) -> None:
        """Test that ex metric are compute and make sense."""

        np.testing.assert_equal(self.coord1, self.coord2)

        ra = self.coord1[:, 0]
        dec = self.coord1[:, 1]

        config = ComputeExPsfTask.ConfigClass()

        config.treecorr.min_sep = 0.01
        config.treecorr.max_sep = 5.0
        config.treecorr.nbins = 1
        config.treecorr.bin_type = "Linear"
        config.treecorr.sep_units = "arcmin"

        task = ComputeExPsfTask(config)
        output1 = task.run(self.de1, self.de2, ra, dec, units="arcmin")

        # At small scale, expect the scalar two-point correlation function
        # to be close to the input variance for de1 and de2. Cross correlation
        # between de1 and de2 should be zeros are they are 2 indendant field.

        np.testing.assert_allclose(output1.metric_E1, 1.0, atol=2e-1)
        np.testing.assert_allclose(output1.metric_E2, 1.0, atol=2e-1)
        np.testing.assert_allclose(output1.metric_Ex, 0.0, atol=2e-1)

        config = ComputeExPsfTask.ConfigClass()

        config.treecorr.min_sep = 5.0
        config.treecorr.max_sep = 600.0
        config.treecorr.nbins = 1
        config.treecorr.bin_type = "Linear"
        config.treecorr.sep_units = "arcmin"

        # At intermediar scale, expect E1>E2>Ex.

        task = ComputeExPsfTask(config)
        output2 = task.run(self.de1, self.de2, ra, dec, units="arcmin")

        np.testing.assert_allclose(output2.metric_E1, 0.20, atol=2e-1)
        np.testing.assert_allclose(output2.metric_E2, 0.05, atol=2e-1)
        np.testing.assert_allclose(output2.metric_Ex, 0.0, atol=2e-1)

        config = ComputeExPsfTask.ConfigClass()

        config.treecorr.min_sep = 600.0
        config.treecorr.max_sep = 1000.0
        config.treecorr.nbins = 1
        config.treecorr.bin_type = "Linear"
        config.treecorr.sep_units = "arcmin"

        # At large scale, expect the scalar two-point correlation function to
        # be all close to 0.

        task = ComputeExPsfTask(config)
        output2 = task.run(self.de1, self.de2, ra, dec, units="arcmin")

        np.testing.assert_allclose(output2.metric_E1, 0.0, atol=2e-1)
        np.testing.assert_allclose(output2.metric_E2, 0.0, atol=2e-1)
        np.testing.assert_allclose(output2.metric_Ex, 0.0, atol=2e-1)


def setup_module(module):
    lsst.utils.tests.init()


class MemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
