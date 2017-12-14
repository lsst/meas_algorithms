#
# LSST Data Management System
#
# Copyright 2017 AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#
from __future__ import absolute_import, division, print_function

import unittest
import numpy as np
import scipy.ndimage

import lsst.afw.image
import lsst.meas.algorithms
import lsst.utils.tests

try:
    type(display)
except NameError:
    display = False


class CorrelatedNoiseDetectionTestCase(lsst.utils.tests.TestCase):
    """Test utilities for dealing with correlated noise when running detection."""

    def testMeasureCorrelatationKernel(self):
        random = np.random.RandomState(50)
        n = 300
        mi = lsst.afw.image.MaskedImage(n, n, dtype=np.float32)
        rms = 5.0
        mi.image.array[:, :] = rms*random.randn(n, n)
        mi.variance.array[:, :] = rms**2
        plane = "DETECTED"
        bitmask = mi.mask.getPlaneBitMask(plane)
        # Randomly mask some pixels and give them higher values so we'll know if
        # masking them out doesn't work.
        m = (random.randn(n, n) > 1)
        mi.mask.array[:, :] = bitmask*m
        mi.image.array[m] += 10.0
        radius = 5
        kernel1 = lsst.meas.algorithms.measureCorrelationKernel(mi, radius=radius, badBitMask=bitmask)
        kernel2 = lsst.meas.algorithms.measureCorrelationKernel(mi, radius=radius, badMaskPlanes=[plane])
        if display:
            from lsst.afw.display import Display
            d = Display(frame=0)
            d.mtv(kernel1)
        self.assertEqual(kernel1.array.shape, (2*radius + 1, 2*radius + 1))
        self.assertImagesEqual(kernel1, kernel2)
        self.assertFloatsAlmostEqual(kernel1.array[5, 5], 1.0, atol=1E-2)
        self.assertFloatsAlmostEqual(kernel1.array[:5, :], 0.0, atol=1E-2)
        self.assertFloatsAlmostEqual(kernel1.array[6:, :], 0.0, atol=1E-2)
        self.assertFloatsAlmostEqual(kernel1.array[:, :5], 0.0, atol=1E-2)
        self.assertFloatsAlmostEqual(kernel1.array[:, 6:], 0.0, atol=1E-2)

    def testFitGeneralDetectionKernel(self):
        psfImage = lsst.afw.detection.GaussianPsf(41, 41, sigma=5.0).computeKernelImage()
        corrImage = lsst.afw.detection.GaussianPsf(15, 15, sigma=2.0).computeKernelImage().convertF()
        detImage = lsst.meas.algorithms.fitGeneralDetectionKernel(psfImage, corrImage, radius=19)
        paddedArray = np.zeros((41, 41), dtype=np.float64)
        paddedArray[1:-1, 1:-1] = detImage.array
        checkArray = scipy.ndimage.convolve(paddedArray, corrImage.array, mode='constant')
        checkImage = lsst.afw.image.ImageD(checkArray.shape[1], checkArray.shape[0])
        checkImage.array[:, :] = checkArray
        if display:
            from lsst.afw.display import Display
            d1 = Display(frame=1)
            d1.mtv(psfImage)
            d2 = Display(frame=2)
            d2.mtv(corrImage)
            d3 = Display(frame=3)
            d3.mtv(detImage)
            d4 = Display(frame=4)
            d4.mtv(checkImage)
        self.assertImagesAlmostEqual(psfImage, checkImage, atol=1E-5)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
