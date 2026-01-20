import unittest

import lsst.utils.tests
import numpy as np
from lsst.afw.geom import makeCdMatrix, makeSkyWcs
from lsst.afw.table import SourceTable
from lsst.geom import Box2I, Extent2I, Point2D, Point2I, SpherePoint, degrees
from lsst.meas.algorithms import DynamicDetectionTask, InsufficientSourcesError
from lsst.meas.algorithms.testUtils import plantSources
from lsst.pex.config import FieldValidationError


class DynamicDetectionTest(lsst.utils.tests.TestCase):
    def setUp(self):
        xy0 = Point2I(12345, 67890)  # xy0 for image
        dims = Extent2I(2345, 2345)  # Dimensions of image
        box = Box2I(xy0, dims)  # Bounding box of image
        sigma = 3.21  # PSF sigma
        buffer = 4.0  # Buffer for star centers around edge
        nSigmaForKernel = 5.0  # Number of PSF sigmas for kernel
        sky = 12345.6  # Sky level
        numStars = 100  # Number of stars
        noise = np.sqrt(sky)*np.pi*sigma**2  # Poisson noise per PSF
        faint = 1.0*noise  # Faintest level for star fluxes
        bright = 100.0*noise  # Brightest level for star fluxes
        starBox = Box2I(box)  # Area on image in which we can put star centers
        starBox.grow(-int(buffer*sigma))
        scale = 1.0e-5*degrees  # Pixel scale

        rng = np.random.Generator(np.random.MT19937(12345))
        stars = [(xx, yy, ff, sigma) for xx, yy, ff in
                 zip(rng.uniform(starBox.getMinX(), starBox.getMaxX(), numStars),
                     rng.uniform(starBox.getMinY(), starBox.getMaxY(), numStars),
                     np.linspace(faint, bright, numStars))]
        self.exposure = plantSources(box, 2*int(nSigmaForKernel*sigma) + 1, sky, stars, True)
        self.exposure.setWcs(makeSkyWcs(crpix=Point2D(0, 0),
                                        crval=SpherePoint(0, 0, degrees),
                                        cdMatrix=makeCdMatrix(scale=scale)))

        self.config = DynamicDetectionTask.ConfigClass()
        self.rtol = 0.1

    def tearDown(self):
        del self.exposure

    def check(self, expectFactor, config):
        schema = SourceTable.makeMinimalSchema()
        task = DynamicDetectionTask(config=config, schema=schema)
        table = SourceTable.make(schema)

        results = task.run(table, self.exposure, expId=12345)
        self.assertFloatsAlmostEqual(results.factor, expectFactor, rtol=self.rtol)

    def testVanilla(self):
        """Dynamic detection used as normal detection."""
        self.check(1.0, self.config)

    def testDynamic(self):
        """Modify the variance plane, and see if the task is able to determine
        what we did.
        """
        factor = 2.0
        self.exposure.maskedImage.variance /= factor
        self.check(1.0/np.sqrt(factor), self.config)

    def testNoWcs(self):
        """Check that dynamic detection runs when the exposure wcs is None."""
        self.exposure.setWcs(None)
        self.check(1.0, self.config)

    def testMinimalSkyObjects(self):
        """Check that dynamic detection runs when there are a relatively small
        number of sky objects.
        """
        config = DynamicDetectionTask.ConfigClass()
        config.skyObjects.nSources = int(0.1 * self.config.skyObjects.nSources)
        self.check(1.0, config)

    def testDynamicNoLimits(self):
        """Test setting the threshold scaling and background tweak limits to
        None (i.e. no limits imposed).
        """
        config = DynamicDetectionTask.ConfigClass()
        config.minThresholdScaleFactor = None
        config.maxThresholdScaleFactor = None
        config.minBackgroundTweak = None
        config.maxBackgroundTweak = None
        self.check(1.0, config)

    def testNoThresholdScaling(self):
        """Check that dynamic detection runs when doThresholdScaling is False.
        """
        config = DynamicDetectionTask.ConfigClass()
        config.doThresholdScaling = False
        self.check(1.0, config)

    def testNoBackgroundTweak(self):
        """Check that dynamic detection runs when doBackgroundTweak is False.
        """
        config = DynamicDetectionTask.ConfigClass()
        config.doBackgroundTweak = False
        self.check(1.0, config)

    def testThresholdScalingAndNoBackgroundTweak(self):
        """Check that dynamic detection runs when both doThresholdScaling and
        doBackgroundTweak are False.
        """
        config = DynamicDetectionTask.ConfigClass()
        config.doThresholdScaling = False
        config.doBackgroundTweak = False
        self.check(1.0, config)

    def testBrightDetectionPass(self):
        """Check that the maximum number of bright detection iterations is
        observed.

        The bright detection loop is forced by setting config.brightMultiplier
        so low such that the entire image is marked as DETECTED.  As such, the
        task is doomed to fail where it tries to lay down sky sources."
        """
        config = DynamicDetectionTask.ConfigClass()
        config.brightMultiplier = 0.05
        config.brightDetectionIterMax = 2
        with self.assertRaisesRegex(
                InsufficientSourcesError, "Insufficient good sky source flux measurements"):
            self.check(1.0, config)

    def testThresholdsOutsideBounds(self):
        """Check that dynamic detection properly sets threshold limits.
        """
        schema = SourceTable.makeMinimalSchema()
        config = DynamicDetectionTask.ConfigClass()
        config.minThresholdScaleFactor = 1.05
        config.maxBackgroundTweak = -0.1
        table = SourceTable.make(schema)
        task = DynamicDetectionTask(config=config, schema=schema)
        task.run(table, self.exposure, expId=12345)

        config = DynamicDetectionTask.ConfigClass()
        config.maxThresholdScaleFactor = 0.99
        config.minBackgroundTweak = 0.1
        table = SourceTable.make(schema)
        task = DynamicDetectionTask(config=config, schema=schema)
        task.run(table, self.exposure, expId=12345)

    def testConfigValidation(self):
        """Check that the field validation is working correctly.
        """
        schema = SourceTable.makeMinimalSchema()
        config = DynamicDetectionTask.ConfigClass()
        config.minThresholdScaleFactor = 1.05
        config.maxThresholdScaleFactor = 1.01
        with self.assertRaisesRegex(
                FieldValidationError, "minThresholdScaleFactor must be <= maxThresholdScaleFactor"):
            DynamicDetectionTask(config=config, schema=schema)

        config = DynamicDetectionTask.ConfigClass()
        config.minBackgroundTweak = 2.0
        config.maxBackgroundTweak = 1.0
        with self.assertRaisesRegex(
                FieldValidationError, "minBackgroundTweak must be <= maxBackgroundTweak"):
            DynamicDetectionTask(config=config, schema=schema)

    def testNoSkyObjects(self):
        """Check that dynamic detection runs when there are no sky objects.

        Notes
        -----
        This test, originally named 'testNoSources', was deactivated on
        DM-39809 because the dynamic detection task is not able to function
        without any sky objects. This test should be reactivated once DM-39826
        has been completed.
        """
        pass


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    import sys
    setup_module(sys.modules['__main__'])
    unittest.main()
