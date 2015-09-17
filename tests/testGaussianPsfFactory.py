#!/usr/bin/env python
from __future__ import division, absolute_import

import unittest
import lsst.utils.tests as utilsTests
from lsst.pex.config import Config, FieldValidationError
from lsst.meas.algorithms import GaussianPsfFactory, SigmaPerFwhm, SingleGaussianPsf, DoubleGaussianPsf

class GaussianPsfFactoryTestCase(unittest.TestCase):
    def testApply(self):
        """Test apply and computeSizeAndSigma methods
        """
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
                    for minSize in (None, max(1, naiveSize-1), naiveSize, naiveSize + 2):
                        if minSize <= 0:
                            continue
                        factory.minSize = minSize
                        for maxSize in (None, naiveSize-1, naiveSize):
                            if maxSize <= 0:
                                continue
                            if None not in (minSize, maxSize) and maxSize < minSize:
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
        """Test the validate method and field-by-field validation
        """
        # test field-by-field validation
        for fieldName in \
            ("size", "sizeFactor", "minSize", "maxSize", "defaultFwhm", "wingFwhmFactor", "wingAmplitude"):
            for value in (-1, 0):
                factory = GaussianPsfFactory()
                self.assertRaises(FieldValidationError, setattr, factory, fieldName, value)

        # test the validate method
        for fieldName in ("sizeFactor", "defaultFwhm", "addWing", "wingFwhmFactor", "wingAmplitude"):
            factory = GaussianPsfFactory()
            setattr(factory, fieldName, None)
            self.assertRaises(Exception, factory.validate)

        for minSize in (None, 5, 9):
            for maxSize in (None, 3, 7):
                factory = GaussianPsfFactory()
                factory.minSize = minSize
                factory.maxSize = maxSize
                if None not in (minSize, maxSize) and maxSize < minSize:
                    self.assertRaises(Exception, factory.validate)
                else:
                    factory.validate() # should not raise

    def testMakeField(self):
        """Test the makeField method
        """
        for addWing in (False, True):
            testConfig = TestConfig()
            testConfig.psfModel.defaultFwhm = 2.7
            testConfig.psfModel.sizeFactor = 3.0
            testConfig.psfModel.minSize = 1
            testConfig.psfModel.maxSize = None
            testConfig.psfModel.addWing = addWing
            psfModel = testConfig.psfModel.apply()
            if addWing:
                self.assertIsInstance(psfModel, DoubleGaussianPsf)
                self.assertAlmostEqual(psfModel.getSigma1(), SigmaPerFwhm * testConfig.psfModel.defaultFwhm)
            else:
                self.assertIsInstance(psfModel, SingleGaussianPsf)
                self.assertAlmostEqual(psfModel.getSigma(), SigmaPerFwhm * testConfig.psfModel.defaultFwhm)

            for fwhm in (None, 2.5):
                predFwhm = testConfig.psfModel.defaultFwhm if fwhm is None else fwhm
                psfModel = testConfig.psfModel.apply(fwhm=fwhm)
                if addWing:
                    self.assertIsInstance(psfModel, DoubleGaussianPsf)
                    self.assertAlmostEqual(psfModel.getSigma1(), SigmaPerFwhm * predFwhm)
                else:
                    self.assertIsInstance(psfModel, SingleGaussianPsf)
                    self.assertAlmostEqual(psfModel.getSigma(), SigmaPerFwhm * predFwhm)


class TestConfig(Config):
    psfModel = GaussianPsfFactory.makeField("PSF model factory")


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(GaussianPsfFactoryTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the utilsTests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
