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

from lsst.meas.algorithms import GaussianPsfFactory, SigmaPerFwhm, SingleGaussianPsf, DoubleGaussianPsf
from lsst.pex.config import Config, FieldValidationError
import lsst.utils.tests


class GaussianPsfFactoryTestCase(lsst.utils.tests.TestCase):

    def testApply(self):
        """Test apply and computeSizeAndSigma methods."""
        factory = GaussianPsfFactory()
        for fixedSize in (None, 5, 6):
            factory.size = fixedSize
            for fwhm in (None, 1, 5, 10):
                if fwhm is None:
                    predFwhm = factory.defaultFwhm
                else:
                    predFwhm = fwhm
                for sizeFactor in (1, 3, 5):
                    factory.sizeFactor = sizeFactor
                    if fixedSize is None:
                        # predicted size without applying limits = fwhm * sizeFactor, increased to be odd
                        naiveSize = int(0.5 + (predFwhm * sizeFactor))
                        if naiveSize % 2 == 0:
                            naiveSize += 1
                    else:
                        naiveSize = fixedSize
                    for minSize in (max(1, naiveSize-1), naiveSize, naiveSize + 2):
                        if minSize <= 0:
                            continue
                        factory.minSize = minSize
                        for maxSize in (naiveSize-1, naiveSize):
                            if maxSize <= 0:
                                continue
                            if maxSize < minSize:
                                continue
                            factory.maxSize = maxSize
                            for addWing in (False, True):
                                factory.addWing = addWing
                                size, sigma = factory.computeSizeAndSigma(predFwhm)
                                self.assertAlmostEqual(sigma, SigmaPerFwhm * predFwhm)

                                if fixedSize is not None:
                                    self.assertEqual(fixedSize, size)
                                else:
                                    if minSize is not None:
                                        self.assertLessEqual(minSize, size)
                                    if maxSize is not None:
                                        self.assertGreaterEqual(maxSize, size)
                                    if None not in (minSize, maxSize) and minSize < size < maxSize:
                                        self.assertEqual(naiveSize, size)
                                        self.assertEqual(size % 2, 1)

                                psfModel = factory.apply(fwhm)
                                if addWing:
                                    self.assertIsInstance(psfModel, DoubleGaussianPsf)
                                    self.assertAlmostEqual(psfModel.getSigma1(), sigma)
                                    self.assertAlmostEqual(
                                        psfModel.getSigma2(), sigma * factory.wingFwhmFactor)
                                    self.assertAlmostEqual(psfModel.getB(), factory.wingAmplitude)
                                else:
                                    self.assertIsInstance(psfModel, SingleGaussianPsf)
                                    self.assertAlmostEqual(psfModel.getSigma(), sigma)

                                psfKernel = psfModel.getKernel()
                                self.assertEqual(size, psfKernel.getHeight())
                                self.assertEqual(size, psfKernel.getWidth())

    def testValidate(self):
        """Test the validate method and field-by-field validation."""
        # test field-by-field validation
        for fieldName in (
            "size",
            "sizeFactor",
            "minSize",
            "maxSize",
            "defaultFwhm",
            "wingFwhmFactor",
            "wingAmplitude",
        ):
            for value in (-1, 0):
                factory = GaussianPsfFactory()
                with self.assertRaises(FieldValidationError):
                    setattr(factory, fieldName, value)

        # test the validate method
        for fieldName in ("sizeFactor", "defaultFwhm", "addWing", "wingFwhmFactor", "wingAmplitude"):
            factory = GaussianPsfFactory()
            setattr(factory, fieldName, None)
            with self.assertRaises(Exception):
                factory.validate()

        for minSize in (None, 5, 9):
            for maxSize in (None, 3, 7):
                factory = GaussianPsfFactory()
                factory.minSize = minSize
                factory.maxSize = maxSize
                if None not in (minSize, maxSize) and maxSize < minSize:
                    with self.assertRaises(Exception):
                        factory.validate()
                else:
                    factory.validate()  # should not raise

    def testMakeField(self):
        """Test the makeField method."""
        for addWing in (False, True):
            psfConfig = PsfConfig()
            psfConfig.psfModel.defaultFwhm = 2.7
            psfConfig.psfModel.sizeFactor = 3.0
            psfConfig.psfModel.minSize = 1
            psfConfig.psfModel.maxSize = None
            psfConfig.psfModel.addWing = addWing
            psfModel = psfConfig.psfModel.apply()
            if addWing:
                self.assertIsInstance(psfModel, DoubleGaussianPsf)
                self.assertAlmostEqual(psfModel.getSigma1(), SigmaPerFwhm * psfConfig.psfModel.defaultFwhm)
            else:
                self.assertIsInstance(psfModel, SingleGaussianPsf)
                self.assertAlmostEqual(psfModel.getSigma(), SigmaPerFwhm * psfConfig.psfModel.defaultFwhm)

            for fwhm in (None, 2.5):
                predFwhm = psfConfig.psfModel.defaultFwhm if fwhm is None else fwhm
                psfModel = psfConfig.psfModel.apply(fwhm=fwhm)
                if addWing:
                    self.assertIsInstance(psfModel, DoubleGaussianPsf)
                    self.assertAlmostEqual(psfModel.getSigma1(), SigmaPerFwhm * predFwhm)
                else:
                    self.assertIsInstance(psfModel, SingleGaussianPsf)
                    self.assertAlmostEqual(psfModel.getSigma(), SigmaPerFwhm * predFwhm)


class PsfConfig(Config):
    psfModel = GaussianPsfFactory.makeField("PSF model factory")


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
