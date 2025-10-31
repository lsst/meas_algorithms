import unittest

import lsst.utils.tests
import numpy as np
from lsst.afw.geom import makeCdMatrix, makeSkyWcs
from lsst.afw.table import SourceTable
from lsst.geom import Box2I, Extent2I, Point2D, Point2I, SpherePoint, degrees
from lsst.meas.algorithms.adaptive_thresholds import AdaptiveThresholdDetectionTask
from lsst.meas.algorithms.testUtils import plantSources


class AdaptiveThresholdDetectionTest(lsst.utils.tests.TestCase):
    def setUp(self):
        xy0 = Point2I(12345, 67890)  # xy0 for image
        dims = Extent2I(2345, 2345)  # Dimensions of image
        self.box = Box2I(xy0, dims)  # Bounding box of image
        self.sigma = 3.21  # PSF sigma
        self.nSigmaForKernel = 5.0  # Number of PSF sigmas for kernel
        self.sky = 12345.6  # Sky level
        noise = np.sqrt(self.sky)*np.pi*self.sigma**2  # Poisson noise per PSF
        self.faint = 1.0*noise  # Faintest level for star fluxes
        self.bright = 100.0*noise  # Brightest level for star fluxes
        self.starBox = Box2I(self.box)  # Area on image in which we can put star centers
        starBoxBuffer = 4.0  # Buffer for star centers around edge
        self.starBox.grow(-int(starBoxBuffer*self.sigma))
        self.scale = 1.0e-5*degrees  # Pixel scale
        self.config = AdaptiveThresholdDetectionTask.ConfigClass()
        self.rtol = 1e-4

    def makeMockExposure(self, numStars=300):
        rng = np.random.Generator(np.random.MT19937(12345))
        stars = [(xx, yy, ff, self.sigma) for xx, yy, ff in
                 zip(rng.uniform(self.starBox.getMinX(), self.starBox.getMaxX(), numStars),
                     rng.uniform(self.starBox.getMinY(), self.starBox.getMaxY(), numStars),
                     np.linspace(self.faint, self.bright, numStars))]
        exposure = plantSources(self.box, 2*int(self.nSigmaForKernel*self.sigma) + 1, self.sky, stars, True)
        exposure.setWcs(makeSkyWcs(crpix=Point2D(0, 0),
                                   crval=SpherePoint(0, 0, degrees),
                                   cdMatrix=makeCdMatrix(scale=self.scale)))
        return exposure

    def tearDown(self):
        # del self.exposure
        del self.config

    def check(self, expectFactor, initialThreshold=None):
        schema = SourceTable.makeMinimalSchema()
        table = SourceTable.make(schema)
        task = AdaptiveThresholdDetectionTask(config=self.config)
        results = task.run(table, self.exposure, initialThreshold=initialThreshold)
        self.assertFloatsAlmostEqual(results.thresholdValue, expectFactor, rtol=self.rtol)

    def testUncrowded(self):
        """Add sparse-ish numbers of stars.
        """
        self.exposure = self.makeMockExposure()
        self.check(5.0)
        self.exposure = self.makeMockExposure(numStars=100)
        self.check(5.0)

    def testCrowded(self):
        """Add enough stars to be considered crowded and test adjusting
        initial detection threshold.
        """
        self.exposure = self.makeMockExposure(numStars=8000)
        self.check(23.8)
        self.exposure = self.makeMockExposure(numStars=8000)
        self.check(50.0, initialThreshold=50.0)
        self.exposure = self.makeMockExposure(numStars=25000)
        self.check(173.4, initialThreshold=20.0)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    import sys
    setup_module(sys.modules['__main__'])
    unittest.main()
