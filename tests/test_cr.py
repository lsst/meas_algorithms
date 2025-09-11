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
import os
import sys
import unittest

import lsst.geom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.meas.algorithms as algorithms
import lsst.meas.algorithms.testUtils as testUtils
import lsst.pex.config as pexConfig
import lsst.utils
import lsst.utils.logging
import lsst.utils.tests

# Increase the number for more verbose messages
lsst.utils.logging.trace_set_at("lsst.algorithms.CR", 3)

try:
    display
except NameError:
    display = False
else:
    import lsst.afw.display as afwDisplay
    afwDisplay.setDefaultMaskTransparency(75)

try:
    afwdataDir = lsst.utils.getPackageDir('afwdata')
    imageFile0 = os.path.join(afwdataDir, "CFHT", "D4", "cal-53535-i-797722_1.fits")
except Exception:
    imageFile0 = None


class CosmicRayTestCase(lsst.utils.tests.TestCase):
    """A test case for Cosmic Ray detection."""

    def setUp(self):
        self.FWHM = 5                   # pixels
        self.psf = algorithms.DoubleGaussianPsf(29, 29, self.FWHM/(2*math.sqrt(2*math.log(2))))

        self.mi = afwImage.MaskedImageF(imageFile0)
        self.XY0 = lsst.geom.PointI(0, 0)  # origin of the subimage we use

        if imageFile0:
            if True:                        # use full image, trimmed to data section
                self.XY0 = lsst.geom.PointI(32, 2)
                self.mi = self.mi.Factory(self.mi, lsst.geom.BoxI(self.XY0, lsst.geom.PointI(2079, 4609)),
                                          afwImage.LOCAL)
                self.nCR = 1076                 # number of CRs we should detect
            else:                               # use sub-image
                if True:
                    self.XY0 = lsst.geom.PointI(824, 140)
                    self.nCR = 10
                else:
                    self.XY0 = lsst.geom.PointI(280, 2750)
                    self.nCR = 2
                self.mi = self.mi.Factory(self.mi, lsst.geom.BoxI(self.XY0, lsst.geom.ExtentI(256, 256),
                                                                  afwImage.LOCAL))
                self.mi.setXY0(lsst.geom.PointI(0, 0))
        else:
            self.nCR = None

        self.mi.getMask().addMaskPlane("DETECTED")

    def tearDown(self):
        del self.mi
        del self.psf

    @unittest.skipUnless(imageFile0, "afwdata not available")
    def testDetection(self):
        """Test CR detection."""
        #
        # Subtract background
        #
        bctrl = afwMath.BackgroundControl(afwMath.Interpolate.NATURAL_SPLINE)
        bctrl.setNxSample(int(self.mi.getWidth()/256) + 1)
        bctrl.setNySample(int(self.mi.getHeight()/256) + 1)
        bctrl.getStatisticsControl().setNumSigmaClip(3.0)
        bctrl.getStatisticsControl().setNumIter(2)

        im = self.mi.getImage()
        try:
            backobj = afwMath.makeBackground(im, bctrl)
        except Exception as e:
            print(e, file=sys.stderr)

            bctrl.setInterpStyle(afwMath.Interpolate.CONSTANT)
            backobj = afwMath.makeBackground(im, bctrl)

        im -= backobj.getImageF()

        if display:
            frame = 0
            disp = afwDisplay.Display(frame=frame)
            disp.mtv(self.mi, title=self._testMethodName + ": Raw")  # raw frame
            if self.mi.getWidth() > 256:
                disp.pan(944 - self.mi.getX0(), 260 - self.mi.getY0())
        #
        # Mask known bad pixels
        #
        badPixels = testUtils.makeDefectList()
        algorithms.interpolateOverDefects(self.mi, self.psf, badPixels)

        stats = afwMath.makeStatistics(self.mi.getImage(), afwMath.MEANCLIP | afwMath.STDEVCLIP)
        background = stats.getValue(afwMath.MEANCLIP)

        crConfig = algorithms.FindCosmicRaysConfig()
        # These tests were written with the original 0.6 default.
        crConfig.cond3_fac2 = 0.6
        crConfig.cond3_fac = 2.5
        crConfig.niteration = 3
        crs = algorithms.findCosmicRays(self.mi, self.psf, background, pexConfig.makePropertySet(crConfig))

        if display:
            frame += 1
            disp = afwDisplay.Display(frame=frame)
            disp.mtv(self.mi, title=self._testMethodName + ": CRs removed")
            if self.mi.getWidth() > 256:
                disp.pan(944 - self.mi.getX0(), 260 - self.mi.getY0())

        print("Detected %d CRs" % len(crs))
        if display and False:
            for cr in crs:
                bbox = cr.getBBox()
                bbox.shift(lsst.geom.ExtentI(-self.mi.getX0(), -self.mi.getY0()))
                disp.line([(bbox.getMinX() - 0.5, bbox.getMinY() - 0.5),
                           (bbox.getMaxX() + 0.5, bbox.getMinY() - 0.5),
                           (bbox.getMaxX() + 0.5, bbox.getMaxY() + 0.5),
                           (bbox.getMinX() - 0.5, bbox.getMaxY() + 0.5),
                           (bbox.getMinX() - 0.5, bbox.getMinY() - 0.5)])

        if self.nCR is not None:
            self.assertEqual(len(crs), self.nCR)


class CosmicRayNullTestCase(unittest.TestCase):
    """A test case for no Cosmic Ray detection."""

    def setUp(self):
        self.FWHM = 5                   # pixels
        self.size = 128

        self.psf = algorithms.DoubleGaussianPsf(29, 29, self.FWHM/(2*math.sqrt(2*math.log(2))))
        self.mi = afwImage.MaskedImageF(128, 128)
        self.mi.set((0, 0, 1))

    def tearDown(self):
        del self.psf
        del self.mi

    def testDetection(self):
        crConfig = algorithms.FindCosmicRaysConfig()
        crs = algorithms.findCosmicRays(self.mi, self.psf, 0.0, pexConfig.makePropertySet(crConfig))
        self.assertEqual(len(crs), 0, "Found %d CRs in empty image" % len(crs))


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
