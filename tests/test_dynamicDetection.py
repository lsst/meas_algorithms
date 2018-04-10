#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

from builtins import zip

import unittest
import numpy as np

import lsst.utils.tests

from lsst.afw.geom import Box2I, Point2I, Point2D, Extent2I, SpherePoint, degrees, makeCdMatrix, makeSkyWcs
from lsst.afw.table import SourceTable
from lsst.meas.algorithms import DynamicDetectionTask
from lsst.meas.algorithms.testUtils import plantSources


class DynamicDetectionTest(lsst.utils.tests.TestCase):
    def setUp(self):
        xy0 = Point2I(12345, 67890)  # xy0 for image
        dims = Extent2I(1234, 567)  # Dimensions of image
        box = Box2I(xy0, dims)  # Bounding box of image
        sigma = 3.21  # PSF sigma
        buffer = 4.0  # Buffer for star centers around edge
        nSigmaForKernel = 5.0  # Number of PSF sigmas for kernel
        sky = 12345.6  # Sky level
        numStars = 100  # Number of stars
        noise = np.sqrt(sky)*np.pi*sigma**2  # Poisson noise per PSF
        faint = 1.0*noise  # Faintest level for star fluxes
        bright = 10.0*noise  # Brightest level for star fluxes
        starBox = Box2I(box)  # Area on image in which we can put star centers
        starBox.grow(-int(buffer*sigma))
        scale = 1.0e-5*degrees  # Pixel scale

        np.random.seed(12345)
        stars = [(xx, yy, ff, sigma) for xx, yy, ff in
                 zip(np.random.uniform(starBox.getMinX(), starBox.getMaxX(), numStars),
                     np.random.uniform(starBox.getMinY(), starBox.getMaxY(), numStars),
                     np.linspace(faint, bright, numStars))]
        self.exposure = plantSources(box, 2*int(nSigmaForKernel*sigma) + 1, sky, stars, True)
        self.exposure.setWcs(makeSkyWcs(crpix=Point2D(0, 0),
                                        crval=SpherePoint(0, 0, degrees),
                                        cdMatrix=makeCdMatrix(scale=scale)))

        self.config = DynamicDetectionTask.ConfigClass()
        # Real simple background subtraction
        self.config.background.useApprox = False
        self.config.background.binSize = max(*self.exposure.getDimensions())
        self.config.skyObjects.nSources = 300

        # Relative tolerance for tweak factor
        # Not sure why this isn't smaller; maybe due to use of Poisson instead of Gaussian noise?
        self.rtol = 0.1

    def tearDown(self):
        del self.exposure

    def check(self, expectFactor):
        schema = SourceTable.makeMinimalSchema()
        task = DynamicDetectionTask(config=self.config, schema=schema)
        table = SourceTable.make(schema)

        results = task.run(table, self.exposure, expId=12345)
        self.assertFloatsAlmostEqual(results.factor, expectFactor, rtol=self.rtol)

    def testVanilla(self):
        """Dynamic detection used as normal detection"""
        self.check(1.0)

    def testDynamic(self):
        """Modify the variance plane, and see if the task is able to determine what we did"""
        factor = 2.0
        self.exposure.maskedImage.variance /= factor
        self.check(1.0/np.sqrt(factor))

    def testNoSources(self):
        self.config.skyObjects.nSources = self.config.minNumSources - 1
        self.check(1.0)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    import sys
    setup_module(sys.modules['__main__'])
    unittest.main()
