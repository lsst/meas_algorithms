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
import unittest

from lsst.afw.image import ExposureF
from lsst.meas.algorithms import SingleGaussianPsf, DoubleGaussianPsf
from lsst.meas.algorithms.installGaussianPsf import InstallGaussianPsfTask, FwhmPerSigma
import lsst.utils.tests


class CandidateMaskingTestCase(lsst.utils.tests.TestCase):
    """Test InstallGaussianPsfTask."""

    def testNoPsf(self):
        """Test InstallGaussianPsfTask when the input exposure has no PSF."""
        for width in (21, 25):
            for fwhm in (2.8, 7.1):
                config = InstallGaussianPsfTask.ConfigClass()
                config.width = width
                config.fwhm = fwhm
                task = InstallGaussianPsfTask(config=config)
                exposure = ExposureF(100, 100)
                task.run(exposure=exposure)
                self.assertTrue(exposure.hasPsf())
                psf = exposure.getPsf()
                psfIm = psf.computeImage()
                self.assertEqual(psfIm.getWidth(), width)
                self.assertEqual(psfIm.getHeight(), width)
                measFwhm = psf.computeShape().getDeterminantRadius()*FwhmPerSigma
                self.assertAlmostEqual(measFwhm, fwhm, delta=1e-3)

    def testMatchDoubleGaussianPsf(self):
        """Test InstallGaussianPsfTask when the input exposure has a DoubleGaussian PSF."""
        config = InstallGaussianPsfTask.ConfigClass()
        task = InstallGaussianPsfTask(config=config)

        for doubleGaussParms in (
            # width, height, inner sigma, outer sigma, outer/inner peak amplitude
            (21, 23, 1.2, 3.5, 0.02),
            (23, 25, 3.5, 9.0, 0.02),
        ):
            exposure = ExposureF(100, 100)
            inPsf = DoubleGaussianPsf(*doubleGaussParms)
            exposure.setPsf(inPsf)
            desWidth, desHeight, innerSigma = doubleGaussParms[0:3]
            task.run(exposure=exposure)
            self.assertTrue(exposure.hasPsf())
            psf = exposure.getPsf()
            psfIm = psf.computeImage()
            self.assertEqual(psfIm.getWidth(), desWidth)
            self.assertEqual(psfIm.getHeight(), desHeight)
            self.assertAlmostEqual(psf.computeShape().getDeterminantRadius(), innerSigma, delta=0.1)

    def testMatchSingleGaussianPsf(self):
        """Test InstallGaussianPsfTask when the input exposure has a single Gaussian PSF."""
        config = InstallGaussianPsfTask.ConfigClass()
        task = InstallGaussianPsfTask(config=config)

        for desWidth, desHeight, desSigma in (
            (21, 23, 1.2),
            (23, 25, 3.5),
        ):
            exposure = ExposureF(100, 100)
            inPsf = SingleGaussianPsf(desWidth, desHeight, desSigma)
            exposure.setPsf(inPsf)
            task.run(exposure=exposure)
            self.assertTrue(exposure.hasPsf())
            psf = exposure.getPsf()
            psfIm = psf.computeImage()
            self.assertEqual(psfIm.getWidth(), desWidth)
            self.assertEqual(psfIm.getHeight(), desHeight)
            self.assertAlmostEqual(psf.computeShape().getDeterminantRadius(), desSigma, delta=1e-3)

    def testBadDim(self):
        for width in (8, 10, 20):
            config = InstallGaussianPsfTask.ConfigClass()
            config.width = width
            with self.assertRaises(Exception):
                config.validate()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
