#!/usr/bin/env python
#
# LSST Data Management System
#
# Copyright 2008-2016  AURA/LSST.
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
from builtins import range
import math
import unittest
import numpy as np

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetect
import lsst.afw.math as afwMath
import lsst.afw.display.ds9 as ds9
import lsst.afw.display.utils as displayUtils
import lsst.meas.algorithms as measAlg
import lsst.pex.exceptions as pexExceptions
import lsst.utils.tests


try:
    type(display)
except NameError:
    display = False


class GaussianPsfTestCase(lsst.utils.tests.TestCase):
    """Test SingleGaussianPsf and DoubleGaussianPsf.

    This test case may be extended to cover any new classes derived from KernelPsf.
    """
    def setUp(self):
        FWHM = 5
        self.ksize = 25                      # size of desired kernel
        sigma = FWHM/(2*math.sqrt(2*math.log(2)))
        self.psfDg = measAlg.DoubleGaussianPsf(self.ksize, self.ksize,
                                               sigma, 1, 0.1)
        self.psfSg = measAlg.SingleGaussianPsf(self.ksize, self.ksize, sigma)

    def tearDown(self):
        del self.psfDg
        del self.psfSg

    def testComputeImage(self):
        """Test the computation of the PSF's image at a point."""

        for psf in [self.psfDg, self.psfSg]:
            ccdXY = afwGeom.Point2D(0, 0)
            kIm = psf.computeImage(ccdXY)

            if False:
                ds9.mtv(kIm)

            self.assertEqual(kIm.getWidth(), self.ksize)
            kIm = psf.computeImage(ccdXY)
            self.assertAlmostEqual(afwMath.makeStatistics(kIm, afwMath.SUM).getValue(), 1.0)

    def testComputeImage2(self):
        """Test the computation of the PSF's image at a point."""

        ccdXY = afwGeom.Point2D(0, 0)
        for psf in [self.psfDg, self.psfSg]:
            kIm = psf.computeImage(ccdXY)
            self.assertEqual(kIm.getWidth(), self.ksize)
            self.assertAlmostEqual(afwMath.makeStatistics(kIm, afwMath.SUM).getValue(), 1.0)

    def testKernel(self):
        """Test the creation of the dgPsf's kernel."""
        for psf in [self.psfDg, self.psfSg]:
            kIm = afwImage.ImageD(psf.getKernel().getDimensions())
            psf.getKernel().computeImage(kIm, False)

            self.assertEqual(kIm.getWidth(), self.ksize)
            self.assertAlmostEqual(afwMath.makeStatistics(kIm, afwMath.SUM).getValue(), 1.0)

        if False:
            ds9.mtv(kIm)

    def testInvalidDgPsf(self):
        """Test parameters of dgPsfs, both valid and not."""
        sigma1, sigma2, b = 1, 0, 0                     # sigma2 may be 0 iff b == 0
        measAlg.DoubleGaussianPsf(self.ksize, self.ksize, sigma1, sigma2, b)

        def badSigma1():
            sigma1 = 0
            measAlg.DoubleGaussianPsf(self.ksize, self.ksize, sigma1, sigma2, b)

        with self.assertRaises(pexExceptions.DomainError):
            badSigma1()

        def badSigma2():
            sigma2, b = 0, 1
            measAlg.DoubleGaussianPsf(self.ksize, self.ksize, sigma1, sigma2, b)

        with self.assertRaises(pexExceptions.DomainError):
            badSigma2()

    def testInvalidSgPsf(self):
        """Test parameters of sgPsfs, both valid and not."""
        sigma = 1.
        measAlg.SingleGaussianPsf(self.ksize, self.ksize, sigma)

        def badSigma1():
            sigma = 0
            measAlg.SingleGaussianPsf(self.ksize, self.ksize, sigma)

        with self.assertRaises(pexExceptions.DomainError):
            badSigma1()

    def testGetImage(self):
        """Test returning a realisation of the dgPsf."""

        for psf in [self.psfSg, self.psfDg]:
            xcen = psf.getKernel().getWidth()//2
            ycen = psf.getKernel().getHeight()//2

            stamps = []
            trueCenters = []
            for x, y in ([10, 10], [9.4999, 10.4999], [10.5001, 10.5001]):
                fx, fy = x - int(x), y - int(y)
                if fx >= 0.5:
                    fx -= 1.0
                if fy >= 0.5:
                    fy -= 1.0

                im = psf.computeImage(afwGeom.Point2D(x, y)).convertF()

                stamps.append(im.Factory(im, True))
                trueCenters.append([xcen + fx, ycen + fy])

            if display:
                mos = displayUtils.Mosaic()     # control mosaics
                ds9.mtv(mos.makeMosaic(stamps))

                for i in range(len(trueCenters)):
                    bbox = mos.getBBox(i)

                    ds9.dot("+",
                            bbox.getMinX() + xcen, bbox.getMinY() + ycen, ctype=ds9.RED, size=1)
                    ds9.dot("+",
                            bbox.getMinX() + trueCenters[i][0], bbox.getMinY() + trueCenters[i][1])

                    ds9.dot("%.2f, %.2f" % (trueCenters[i][0], trueCenters[i][1]),
                            bbox.getMinX() + xcen, bbox.getMinY() + 2)

    def testKernelPsf(self):
        """Test creating a Psf from a Kernel."""

        x, y = 10.4999, 10.4999
        ksize = 15
        sigma1 = 1
        #
        # Make a PSF from that kernel
        #
        kPsf = measAlg.KernelPsf(afwMath.AnalyticKernel(ksize, ksize,
                                                        afwMath.GaussianFunction2D(sigma1, sigma1)))

        kIm = kPsf.computeImage(afwGeom.Point2D(x, y))
        #
        # And now via the dgPsf model
        #
        dgPsf = measAlg.DoubleGaussianPsf(ksize, ksize, sigma1)
        dgIm = dgPsf.computeImage(afwGeom.Point2D(x, y))
        #
        # Check that they're the same
        #
        diff = type(kIm)(kIm, True)
        diff -= dgIm
        stats = afwMath.makeStatistics(diff, afwMath.MAX | afwMath.MIN)
        self.assertAlmostEqual(stats.getValue(afwMath.MAX), 0.0, places=16)
        self.assertAlmostEqual(stats.getValue(afwMath.MIN), 0.0, places=16)

        for pad in [-2, 4, 0]:
            resizedKPsf = kPsf.resized(ksize + pad, ksize + pad)
            self.assertEqual(resizedKPsf.computeBBox().getDimensions(),
                             afwGeom.Extent2I(ksize + pad, ksize + pad))
            self.assertEqual(resizedKPsf.getKernel().getKernelParameters(),
                             kPsf.getKernel().getKernelParameters())
            self._compareKernelImages(kPsf, resizedKPsf)
        if display:
            mos = displayUtils.Mosaic()
            mos.setBackground(-0.1)
            ds9.mtv(mos.makeMosaic([kIm, dgIm, diff], mode="x"), frame=1)

    def testResize(self):
        """Test that resized Single and Double Gaussian PSFs have
        same model parameters, but new kernel dimensions."""

        for lengthNew in [1, 11, 99]:
            # Test Double Gaussian
            psfResized = self.psfDg.resized(lengthNew, lengthNew)
            self.assertEqual(psfResized.getSigma1(), self.psfDg.getSigma1())
            self.assertEqual(psfResized.getSigma2(), self.psfDg.getSigma2())
            self.assertEqual(psfResized.getB(), self.psfDg.getB())
            self._compareKernelImages(psfResized, self.psfDg)

            self.assertEqual(psfResized.getKernel().getWidth(), lengthNew)
            self.assertEqual(psfResized.getKernel().getHeight(), lengthNew)

            # Test Single Gaussian Parameters
            psfResized = self.psfSg.resized(lengthNew, lengthNew)
            self.assertEqual(psfResized.getSigma(), self.psfSg.getSigma())
            self._compareKernelImages(psfResized, self.psfSg)

            self.assertEqual(psfResized.getKernel().getWidth(), lengthNew)
            self.assertEqual(psfResized.getKernel().getHeight(), lengthNew)

    def _compareKernelImages(self, psf1, psf2):
        """Test that overlapping portions of kernel images are identical
        """
        im1 = psf1.computeKernelImage()
        im2 = psf2.computeKernelImage()
        bboxIntersection = im1.getBBox()
        bboxIntersection.clip(im2.getBBox())
        im1Intersection = afwImage.ImageD(im1, bboxIntersection)
        im2Intersection = afwImage.ImageD(im2, bboxIntersection)
        scale1 = im1.getArray().sum() / im1Intersection.getArray().sum()
        scale2 = im2.getArray().sum() / im2Intersection.getArray().sum()
        im1Arr = scale1 * im1Intersection.getArray()
        im2Arr = scale2 * im2Intersection.getArray()
        self.assertTrue(np.allclose(im1Arr, im2Arr),
                        "kernel images %s, %s do not match" % (im1Arr, im2Arr))

    def testComputeBBox(self):
        """Test that computeBBox returns same bbox as kernel
        """
        for psf in [self.psfDg, self.psfSg]:
            self.assertEqual(psf.computeBBox(), psf.getKernel().getBBox())

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
