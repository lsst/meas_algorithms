from __future__ import print_function, division

import numpy as np

import unittest
import lsst.utils.tests

import lsst.afw.display
import lsst.afw.image
import lsst.afw.detection
import lsst.afw.geom

import lsst.meas.algorithms


class TrailTestCase(lsst.utils.tests.TestCase):
    """A test case for trail detection"""

    @lsst.utils.tests.debugger(Exception)
    def testTrail(self):
        width, height = 2048, 2048
        x0, y0 = 12345, 6789
        psfSigma = 4.321
        psfSize = 2*int(5.0*psfSigma) + 1
        slope, intercept = 1.23, -123.45
        trailWidth = psfSigma # 5.4321
        trailPeak = 100.0
        noise = 10.0
        seed = 987654321

        exposure = lsst.afw.image.ExposureF(width, height)
        exposure.setXY0(lsst.afw.geom.PointI(x0, y0))
        exposure.setPsf(lsst.afw.detection.GaussianPsf(psfSize, psfSize, psfSigma))
        image = exposure.getMaskedImage().getImage().getArray()

        xx, yy = np.meshgrid(np.arange(width), np.arange(height))

        """
        y = slope*x + intercept
        0 = slope*x + -1*y + intercept
        distance to line = (slope*x - y + intercept)/sqrt(slope^2 + 1)
        """
        distance = (slope*xx - yy + intercept)/np.hypot(slope, -1.0)
        image[:] = trailPeak*np.exp(-0.5*distance**2/trailWidth**2)

        rng = np.random.RandomState(seed)
        image += noise*rng.normal(size=image.shape)
        exposure.getMaskedImage().getVariance().set(noise**2)

        lsst.afw.display.getDisplay(1).mtv(exposure, title="Input image")

        config = lsst.meas.algorithms.TrailConfig()

        task = lsst.meas.algorithms.TrailTask(config=config, name="trail")
        trails = task.run(exposure)

        lsst.afw.display.getDisplay(2).mtv(exposure, title="Output image")
        import pdb;pdb.set_trace()
        trails.display(exposure.getDimensions(), frame=2)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module, backend="virtualDevice"):
    lsst.utils.tests.init()
    try:
        lsst.afw.display.setDefaultBackend(backend)
    except:
        print("Unable to configure display backend: %s" % backend)


if __name__ == "__main__":
    import sys

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--backend', type=str, default="virtualDevice",
                        help="The backend to use, e.g. 'ds9'. Be sure to 'setup display_<backend>'")
    parser.add_argument('unittest', nargs='*')
    args = parser.parse_args()
    sys.argv[1:] = args.unittest

    setup_module(sys.modules[__name__], backend=args.backend)
    unittest.main()
