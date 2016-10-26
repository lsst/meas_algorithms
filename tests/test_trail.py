from __future__ import print_function, division

import sys
from contextlib import contextmanager
import numpy as np

import unittest
import lsst.utils.tests

import lsst.afw.display
import lsst.afw.image
import lsst.afw.detection
import lsst.afw.geom

import lsst.meas.algorithms

from lsst.pipe.base import Struct
from lsst.pex.config import Config, Field

import lsst.log


from contextlib import contextmanager
from timeit import default_timer
@contextmanager
def timer():
    start = default_timer()
    elapser = lambda: default_timer() - start
    yield lambda: elapser()
    end = default_timer()
    elapser = lambda: end - start


class Trail(Struct):
    def __init__(self, peak, width, slope, intercept):
        Struct.__init__(self, peak=peak, width=width, slope=slope, intercept=intercept)

    def func(self, distance):
        raise NotImplementedError("Subclasses should define")

    def makeImage(self, width, height):
        xx, yy = np.meshgrid(np.arange(width), np.arange(height))
        distance = (self.slope*xx - yy + self.intercept)/np.hypot(self.slope, -1.0)
        return self.peak*self.func(distance)

class Satellite(Trail):
    def func(self, distance):
        return np.exp(-0.5*distance**2/self.width**2)

class WingedSatellite(Satellite):
    def func(self, distance):
        return Satellite.func(distance) + 0.1*self.peak*np.exp(-0.5*distance**2/(2.0*self.width)**2)

class Aircraft(Trail):
    def func(self, distance):
        return np.where(np.abs(distance) < self.width, 1.0, 0.0)


class CheckConfig(Config):
    width = Field(dtype=int, default=2048, doc="Width of image")
    height = Field(dtype=int, default=2048, doc="Height of image")
    x0 = Field(dtype=int, default=12345, doc="x offset for image")
    y0 = Field(dtype=int, default=6789, doc="y offset for image")
    psfSigma = Field(dtype=float, default=4.321, doc="Gaussian sigma for PSF")
    noise = Field(dtype=float, default=10.0, doc="Noise to add to image")
    seed = Field(dtype=int, default=987654321, doc="Seed for RNG")
    meanFactor = Field(dtype=float, default=0.01, doc="Factor of noise for threshold on cleaned image mean")
    stdFactor = Field(dtype=float, default=1.1, doc="Factor of noise for threshold on cleaned image stdev")
    numStars = Field(dtype=int, default=100, doc="Number of stars in image")
    starFlux = Field(dtype=float, default=1000.0, doc="Flux of stars")
    starSize = Field(dtype=float, default=5.0, doc="Size of star in units of psfSigma")
    thetaTolerance = Field(dtype=float, default=0.05, doc="Absolute tolerance of recovered angle (radians)")
    radiusTolerance = Field(dtype=float, default=10.0, doc="Absolute tolerance of recovered radius (pixels)")
    detectionThreshold = Field(dtype=float, default=3.0, doc="Threshold for detection")


class TrailTestCase(lsst.utils.tests.TestCase):
    """A test case for trail detection"""

#    @lsst.utils.tests.debugger(Exception)
    def checkTrails(self, trailList, config=None, check=None):
        if config is None:
            config = lsst.meas.algorithms.TrailConfig()
        if check is None:
            check = CheckConfig()

        exposure = lsst.afw.image.ExposureF(check.width, check.height)
        exposure.setXY0(lsst.afw.geom.PointI(check.x0, check.y0))
        psfSize = 2*int(5.0*check.psfSigma) + 1
        exposure.setPsf(lsst.afw.detection.GaussianPsf(psfSize, psfSize, check.psfSigma))
        image = exposure.getMaskedImage().getImage()
        array = image.getArray()
        rng = np.random.RandomState(check.seed)
        noise = image.Factory(image.getBBox())
        noise.getArray()[:] = check.noise*rng.normal(size=array.shape)

        stars = image.Factory(image.getBBox())
        size = check.starSize*check.psfSigma
        for ii in range(check.numStars):
            xCenter = rng.uniform(0, check.width)
            yCenter = rng.uniform(0, check.height)
            box = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(int(xCenter), int(yCenter)),
                                      lsst.afw.geom.Extent2I(1, 1))
            box.grow(int(size + 0.5))
            box.clip(exposure.getBBox(lsst.afw.image.LOCAL))
            xx, yy = np.meshgrid(np.arange(box.getMinX(), box.getMaxX() + 1),
                                 np.arange(box.getMinY(), box.getMaxY() + 1))

            distance2 = (xx - xCenter)**2 + (yy - yCenter)**2
            stars[box].getArray()[:] += check.starFlux*np.exp(-0.5*distance2/check.psfSigma**2)

        mask = exposure.getMaskedImage().getMask()
        mask.getArray()[:] = np.where(stars.getArray() > check.detectionThreshold*check.noise,
                                      mask.getPlaneBitMask("DETECTED"), 0)

        exposure.getMaskedImage().getVariance().set(check.noise**2)

        for trail in trailList:
            array += trail.makeImage(check.width, check.height)

        image += stars
        image += noise

        lsst.log.Log.getLogger("trail").setLevel(lsst.log.DEBUG)
        task = lsst.meas.algorithms.TrailTask(config=config, name="trail")

        with timer() as elapsed:
            found = task.run(exposure)
        task.log.info("Elapsed time: %f", elapsed())

        def trailToHesse(trail):
            """Convert our trail to Hesse normal form

            The Hesse normal form is: x.cos(theta) + y.sin(theta) = r
            Our trails have been specified in the form: y = slope.x + intercept

            It's more convenient to compare the trails in Hesse form because
            the radius (r) and angle (theta) can be compared directly with
            a simple threshold, while the slope and intercept are more
            difficult to compare because the thresholds would change depending
            upon the slope.
            """
            theta = np.arctan(-1.0/trail.slope)
            radius = trail.intercept*np.sin(theta)
            if radius < 0:
                radius *= -1
                theta += np.pi
            if theta < 0:
                theta += 2*np.pi
            return Struct(r=radius, theta=theta*lsst.afw.geom.radians)

        self.assertEqual(len(found), len(trailList))
        # Sort the trail lists by angle to ensure we're comparing the same trail in each.
        for inTrail, outTrail in zip(sorted([trailToHesse(x) for x in trailList], key=lambda x: x.theta),
                                     sorted(found, key=lambda x: x.theta)):
            self.assertClose(outTrail.theta.asRadians(), inTrail.theta.asRadians(), rtol=None,
                             atol=check.thetaTolerance)
            self.assertClose(outTrail.r, inTrail.r, rtol=None, atol=check.radiusTolerance)

        mask = exposure.getMaskedImage().getMask()
        maskVal = mask.getPlaneBitMask("TRAIL")
        clean = mask.getArray() & maskVal == 0
        image -= stars
        self.assertLess(np.abs(array[clean].mean()), check.meanFactor*check.noise)
        self.assertLess(array[clean].std(), check.stdFactor*check.noise)
        image -= noise
        self.assertLess(array[clean].max(), config.maskFluxFraction*check.noise)

    def testSatellite(self):
        check = CheckConfig()
        trails = [
            Satellite(100.0, check.psfSigma, 1.2345, -123.45),
            Satellite(30.0, check.psfSigma, -1.2345, 2345.67),
            Satellite(10.0, check.psfSigma, 0.12345, 432.1),
            Satellite(5.0, check.psfSigma, 5.4321, -2345.67),
        ]

        self.checkTrails(trails)


    def testAircraft(self):
        trails = [
            Aircraft(300.0, 32.1, 1.02030, 123.45),
            Aircraft(100.0, 54.3, -1.2345, 2345.67),
        ]

        config = lsst.meas.algorithms.TrailConfig()
        config.widths = [40.0, 70.0, 100.0]
        config.sigmaSmooth = 2.0
        config.kernelSigma = 9.0
        config.kernelWidth = 15
        config.bins = 8

        self.checkTrails(trails, config)

    def testNone(self):
        self.checkTrails([])


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module, backend="virtualDevice"):
    lsst.utils.tests.init()
    try:
        lsst.afw.display.setDefaultBackend(backend)
    except:
        print("Unable to configure display backend: %s" % backend)


@contextmanager
def profiling(restriction=30, stream=sys.stdout):
    import cProfile, pstats
    prof = cProfile.Profile()
    prof.enable()
    try:
        yield prof
    finally:
        prof.disable()
        pstats.Stats(prof, stream=stream).sort_stats("cumulative").print_stats(restriction)


@contextmanager
def nullContext():
    yield


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--backend', type=str, default="virtualDevice",
                        help="The backend to use, e.g. 'ds9'. Be sure to 'setup display_<backend>'")
    parser.add_argument('--profile', default=False, action="store_true", help="Activate profiling")
    parser.add_argument('unittest', nargs='*')
    args = parser.parse_args()
    sys.argv[1:] = args.unittest

    setup_module(sys.modules[__name__], backend=args.backend)

    profiler = profiling if args.profile else nullContext
    with profiler():
        unittest.main(exit=False)
