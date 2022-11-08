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

import math
import unittest
import numpy as np

import lsst.geom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg
import lsst.pex.exceptions as pexExceptions
import lsst.utils.tests

try:
    display
except NameError:
    display = False
else:
    import lsst.afw.display as afwDisplay
    afwDisplay.setDefaultMaskTransparency(75)


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
            ccdXY = lsst.geom.Point2D(0, 0)
            kIm = psf.computeImage(ccdXY)

            if False:
                afwDisplay.Display(frame=1).mtv(kIm, title=self._testMethodName + ": kIm")

            self.assertEqual(kIm.getWidth(), self.ksize)
            kIm = psf.computeImage(ccdXY)
            self.assertAlmostEqual(afwMath.makeStatistics(kIm, afwMath.SUM).getValue(), 1.0)

    def testComputeImage2(self):
        """Test the computation of the PSF's image at a point.
        """
        ccdXY = lsst.geom.Point2D(0, 0)
        for psf in [self.psfDg, self.psfSg]:
            kIm = psf.computeImage(ccdXY)
            self.assertEqual(kIm.getWidth(), self.ksize)
            self.assertAlmostEqual(afwMath.makeStatistics(kIm, afwMath.SUM).getValue(), 1.0)

    def testKernel(self):
        """Test the creation of the dgPsf's kernel.
        """
        for psf in [self.psfDg, self.psfSg]:
            kIm = afwImage.ImageD(psf.getKernel().getDimensions())
            psf.getKernel().computeImage(kIm, False)

            self.assertEqual(kIm.getWidth(), self.ksize)
            self.assertAlmostEqual(afwMath.makeStatistics(kIm, afwMath.SUM).getValue(), 1.0)

        if False:
            afwDisplay.Display(frame=2).mtv(kIm, title=self._testMethodName + ": kIm")

    def testInvalidDgPsf(self):
        """Test parameters of dgPsfs, both valid and not.
        """
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
        """Test parameters of sgPsfs, both valid and not.
        """
        sigma = 1.
        measAlg.SingleGaussianPsf(self.ksize, self.ksize, sigma)

        def badSigma1():
            sigma = 0
            measAlg.SingleGaussianPsf(self.ksize, self.ksize, sigma)

        with self.assertRaises(pexExceptions.DomainError):
            badSigma1()

    def testGetImage(self):
        """Test returning a realisation of the dgPsf.
        """
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

                im = psf.computeImage(lsst.geom.Point2D(x, y)).convertF()

                stamps.append(im.Factory(im, True))
                trueCenters.append([xcen + fx, ycen + fy])

            if display:
                mos = afwDisplay.utils.Mosaic()     # control mosaics
                disp = afwDisplay.Display(frame=0)
                disp.mtv(mos.makeMosaic(stamps), title=self._testMethodName + ": mosaic")

                for i in range(len(trueCenters)):
                    bbox = mos.getBBox(i)

                    disp.dot("+",
                             bbox.getMinX() + xcen, bbox.getMinY() + ycen, ctype=afwDisplay.RED, size=1)
                    disp.dot("+",
                             bbox.getMinX() + trueCenters[i][0], bbox.getMinY() + trueCenters[i][1])

                    disp.dot("%.2f, %.2f" % (trueCenters[i][0], trueCenters[i][1]),
                             bbox.getMinX() + xcen, bbox.getMinY() + 2)

    def testKernelPsf(self):
        """Test creating a Psf from a Kernel.
        """
        x, y = 10.4999, 10.4999
        ksize = 15
        sigma1 = 1
        #
        # Make a PSF from that kernel
        #
        kPsf = measAlg.KernelPsf(afwMath.AnalyticKernel(ksize, ksize,
                                                        afwMath.GaussianFunction2D(sigma1, sigma1)))

        kIm = kPsf.computeImage(lsst.geom.Point2D(x, y))
        #
        # And now via the dgPsf model
        #
        dgPsf = measAlg.DoubleGaussianPsf(ksize, ksize, sigma1)
        dgIm = dgPsf.computeImage(lsst.geom.Point2D(x, y))
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
            self.assertEqual(
                resizedKPsf.computeBBox(resizedKPsf.getAveragePosition()).getDimensions(),
                lsst.geom.Extent2I(ksize + pad, ksize + pad)
            )
            self.assertEqual(resizedKPsf.getKernel().getKernelParameters(),
                             kPsf.getKernel().getKernelParameters())
            self._compareKernelImages(kPsf, resizedKPsf)
        if display:
            mos = afwDisplay.utils.Mosaic()
            mos.setBackground(-0.1)
            afwDisplay.Display(frame=1).mtv(mos.makeMosaic([kIm, dgIm, diff], mode="x"),
                                            title=self._testMethodName + ": mosaic")

    def testResize(self):
        """Test that resized Single and Double Gaussian PSFs have
        same model parameters, but new kernel dimensions.
        """

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
        """Test that overlapping portions of kernel images are identical.
        """
        im1 = psf1.computeKernelImage(psf1.getAveragePosition())
        im2 = psf2.computeKernelImage(psf2.getAveragePosition())
        bboxIntersection = im1.getBBox()
        bboxIntersection.clip(im2.getBBox())
        im1Intersection = afwImage.ImageD(im1, bboxIntersection)
        im2Intersection = afwImage.ImageD(im2, bboxIntersection)
        scale1 = im1.getArray().sum()/im1Intersection.getArray().sum()
        scale2 = im2.getArray().sum()/im2Intersection.getArray().sum()
        im1Arr = scale1*im1Intersection.getArray()
        im2Arr = scale2*im2Intersection.getArray()
        self.assertTrue(np.allclose(im1Arr, im2Arr),
                        "kernel images %s, %s do not match" % (im1Arr, im2Arr))

    def testComputeBBox(self):
        """Test that computeBBox returns same bbox as kernel.
        """
        for psf in [self.psfDg, self.psfSg]:
            self.assertEqual(
                psf.computeBBox(psf.getAveragePosition()),
                psf.getKernel().getBBox()
            )


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
