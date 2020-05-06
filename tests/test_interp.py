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

import os
import unittest
import math
import numpy as np

import lsst.geom
import lsst.afw.image as afwImage
import lsst.meas.algorithms as algorithms
import lsst.utils.tests
from lsst.daf.base import PropertyList

try:
    type(display)
except NameError:
    display = False
else:
    import lsst.afw.display as afwDisplay
    afwDisplay.setDefaultMaskTransparency(75)

# Determine if we have afwdata
try:
    afwdataDir = lsst.utils.getPackageDir('afwdata')
except Exception:
    afwdataDir = None

TESTDIR = os.path.abspath(os.path.dirname(__file__))


class DefectsTestCase(lsst.utils.tests.TestCase):
    """Tests for collections of Defect."""

    def assertMetadata(self, first, second):
        """Compare the metadata associated with Defects"""

        # Must strip out DATE metadata before comparison
        meta1 = first.getMetadata()
        meta2 = second.getMetadata()
        for d in (meta1, meta2):
            for k in ("DATE", "CALIB_CREATION_DATE", "CALIB_CREATION_TIME"):
                if k in d:
                    del d[k]

        self.assertEqual(meta1, meta2)
        meta1["NEW"] = "additional header"
        self.assertNotEqual(first.getMetadata(), second.getMetadata())
        del meta1["NEW"]

    def test_defects(self):
        defects = algorithms.Defects()

        defects.append(algorithms.Defect(lsst.geom.Box2I(lsst.geom.Point2I(5, 6),
                                         lsst.geom.Point2I(41, 50))))

        defects.append(lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                                       lsst.geom.Point2I(4, 5)))
        defects.append(lsst.geom.Point2I(10, 12))
        defects.append(afwImage.DefectBase(lsst.geom.Box2I(lsst.geom.Point2I(100, 200),
                                                           lsst.geom.Extent2I(5, 5))))
        self.assertEqual(len(defects), 4)

        for d in defects:
            self.assertIsInstance(d, algorithms.Defect)

        # Transposition
        transposed = defects.transpose()
        self.assertEqual(len(transposed), len(defects))
        self.assertEqual(transposed[0].getBBox(),
                         lsst.geom.Box2I(lsst.geom.Point2I(6, 5),
                                         lsst.geom.Extent2I(45, 37)))

        # Serialization round trip
        meta = PropertyList()
        meta["TESTHDR"] = "testing"
        defects.setMetadata(meta)

        table = defects.toFitsRegionTable()
        defects2 = algorithms.Defects.fromTable(table)
        self.assertEqual(defects2, defects)

        # via FITS
        with lsst.utils.tests.getTempFilePath(".fits") as tmpFile:
            defects.writeFits(tmpFile)
            defects2 = algorithms.Defects.readFits(tmpFile)

        # Equality tests the bounding boxes so metadata is tested separately.
        self.assertEqual(defects2, defects)
        self.assertMetadata(defects2, defects)

        # via text file
        with lsst.utils.tests.getTempFilePath(".ecsv") as tmpFile:
            defects.writeText(tmpFile)
            defects2 = algorithms.Defects.readText(tmpFile)

        # Equality tests the bounding boxes so metadata is tested separately.
        self.assertEqual(defects2, defects)
        self.assertMetadata(defects2, defects)

        # Check bad values
        with self.assertRaises(ValueError):
            defects.append(lsst.geom.Box2D(lsst.geom.Point2D(0., 0.),
                                           lsst.geom.Point2D(3.1, 3.1)))
        with self.assertRaises(ValueError):
            defects.append("defect")

    def testAstropyRegion(self):
        """Read a FITS region file created by Astropy regions."""

        with self.assertLogs():
            defects = algorithms.Defects.readFits(os.path.join(TESTDIR, "data", "fits_region.fits"))

        # Should be able to read 3 regions from the file
        self.assertEqual(len(defects), 3)

    def testLsstTextfile(self):
        """Read legacy LSST text file format"""
        with lsst.utils.tests.getTempFilePath(".txt") as tmpFile:
            with open(tmpFile, "w") as fh:
                print("""# X0  Y0  width height
     996        0       56       24
       0     4156     2048       20
       0        0       17     4176
    1998     4035       50      141
    1023        0        2     4176
    2027        0       21     4176
       0     4047       37      129
# Some rows without fixed column widths
14 20 2000 50
10 10 10 10
""", file=fh)

            defects = algorithms.Defects.readLsstDefectsFile(tmpFile)

            self.assertEqual(len(defects), 9)
            self.assertEqual(defects[3].getBBox(), lsst.geom.Box2I(lsst.geom.Point2I(1998, 4035),
                                                                   lsst.geom.Extent2I(50, 141)))

    def test_normalize_defects(self):
        """A test for the lsst.meas.algorithms.Defect.normalize() method.
        """
        defects = algorithms.Defects()
        #  First series of 1-pixel contiguous defects
        defects.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(15, 1), dimensions=lsst.geom.Extent2I(1, 1)))
        defects.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(15, 2), dimensions=lsst.geom.Extent2I(1, 1)))
        defects.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(15, 3), dimensions=lsst.geom.Extent2I(1, 1)))
        defects.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(15, 4), dimensions=lsst.geom.Extent2I(1, 1)))
        defects.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(15, 5), dimensions=lsst.geom.Extent2I(1, 1)))

        # Second series of 1-pixel contiguos defects
        defects.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(20, 11), dimensions=lsst.geom.Extent2I(1, 1)))
        defects.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(20, 12), dimensions=lsst.geom.Extent2I(1, 1)))
        defects.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(20, 13), dimensions=lsst.geom.Extent2I(1, 1)))
        defects.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(20, 14), dimensions=lsst.geom.Extent2I(1, 1)))
        defects.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(20, 15), dimensions=lsst.geom.Extent2I(1, 1)))

        # the normalize() method should replace the 10 defects above by just
        # two defect boxes with adjusted dimensions
        newDefects = defects.normalize()

        boxesMeasured = []
        for defect in newDefects:
            boxesMeasured.append(defect.getBBox())

        # The normalizing function should have created the following two boxes out
        # of the individual 1-pixel defects from above
        expectedDefects = [lsst.geom.Box2I(corner=lsst.geom.Point2I(15, 1),
                                           dimensions=lsst.geom.Extent2I(1, 5)),
                           lsst.geom.Box2I(corner=lsst.geom.Point2I(20, 11),
                                           dimensions=lsst.geom.Extent2I(1, 5))]

        self.assertEqual(len(expectedDefects), len(boxesMeasured))
        for expDef, measDef in zip(expectedDefects, boxesMeasured):
            self.assertEqual(expDef, measDef)

        # Normalize two distinct sets of Defects and ensure they compare to the same thing
        defects = algorithms.Defects()
        # Set 1
        defects.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(25, 1), dimensions=lsst.geom.Extent2I(1, 1)))
        defects.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(25, 2), dimensions=lsst.geom.Extent2I(1, 1)))
        defects.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(25, 3), dimensions=lsst.geom.Extent2I(1, 1)))
        defects.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(25, 4), dimensions=lsst.geom.Extent2I(1, 1)))
        defects.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(25, 5), dimensions=lsst.geom.Extent2I(1, 1)))
        defects.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(25, 6), dimensions=lsst.geom.Extent2I(1, 1)))
        defects.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(25, 7), dimensions=lsst.geom.Extent2I(1, 1)))
        defects.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(25, 8), dimensions=lsst.geom.Extent2I(1, 1)))

        # Set 2
        defects2 = algorithms.Defects()
        defects2.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(25, 1), dimensions=lsst.geom.Extent2I(1, 5)))
        defects2.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(25, 5), dimensions=lsst.geom.Extent2I(1, 4)))

        # normalize both sets
        newDefects = defects.normalize()
        newDefects2 = defects2.normalize()

        self.assertEqual(len(newDefects), len(newDefects2))

        boxesMeasured, boxesMeasured2 = [], []
        for defect, defect2 in zip(newDefects, newDefects2):
            boxesMeasured.append(defect.getBBox())
            boxesMeasured2.append(defect2.getBBox())

        expectedDefects = [lsst.geom.Box2I(corner=lsst.geom.Point2I(25, 1),
                                           dimensions=lsst.geom.Extent2I(1, 8))]

        self.assertEqual(len(expectedDefects), len(boxesMeasured))
        for expDef, measDef in zip(expectedDefects, boxesMeasured):
            self.assertEqual(expDef, measDef)

        self.assertEqual(len(expectedDefects), len(boxesMeasured2))
        for expDef, measDef in zip(expectedDefects, boxesMeasured2):
            self.assertEqual(expDef, measDef)


class InterpolationTestCase(lsst.utils.tests.TestCase):
    """A test case for interpolation."""

    def setUp(self):
        self.FWHM = 5
        self.psf = algorithms.DoubleGaussianPsf(15, 15, self.FWHM/(2*math.sqrt(2*math.log(2))))
        maskedImageFile = os.path.join(afwdataDir, "CFHT", "D4", "cal-53535-i-797722_1.fits")

        self.mi = afwImage.MaskedImageF(maskedImageFile)
        if False:                       # use sub-image?
            self.mi = self.mi.Factory(self.mi, afwImage.BBox(afwImage.PointI(760, 20), 256, 256))
        self.mi.getMask().addMaskPlane("INTERP")

        measAlgorithmsDir = lsst.utils.getPackageDir('meas_algorithms')
        self.badPixels = algorithms.Defects.readText(os.path.join(measAlgorithmsDir,
                                                                  "policy", "BadPixels.ecsv"))

    def tearDown(self):
        del self.mi
        del self.psf
        del self.badPixels

    @unittest.skipUnless(afwdataDir, "afwdata not available")
    def testDetection(self):
        """Test Interp algorithms."""

        if display:
            frame = 0
            afwDisplay.Display(frame=frame).mtv(self.mi, title=self._testMethodName + ": Original")

        algorithms.interpolateOverDefects(self.mi, self.psf, self.badPixels)

        if display:
            frame += 1
            afwDisplay.Display(frame=frame).mtv(self.mi, title=self._testMethodName + ": Interpolated")
            frame += 1
            afwDisplay.Display(frame=frame).mtv(self.mi.getVariance(),
                                                title=self._testMethodName + ": Variance")

    @unittest.skipUnless(afwdataDir, "afwdata not available")
    def test818(self):
        """A test case for #818; the full test is in /lsst/DC3root/ticketFiles/818"""

        badPixels = []
        defects = [((82, 663), 6, 8),
                   ((83, 659), 9, 6),
                   ((85, 660), 10, 11),
                   ((87, 669), 3, 3),
                   ]

        for xy0, width, height in defects:
            x0, y0 = xy0
            bbox = lsst.geom.BoxI(lsst.geom.PointI(x0, y0), lsst.geom.ExtentI(width, height))
            badPixels.append(algorithms.Defect(bbox))

        mi = afwImage.MaskedImageF(517, 800)

        algorithms.interpolateOverDefects(mi, self.psf, badPixels)

    @unittest.skipUnless(afwdataDir, "afwdata not available")
    def test1295(self):
        """A test case for #1295 (failure to interpolate over groups of defects."""
        im = afwImage.ImageF(lsst.geom.ExtentI(100, 100))
        mi = afwImage.makeMaskedImage(im)
        mi.set(100)
        flat = afwImage.ImageF(im.getDimensions())
        flat.set(1)
        flat[50:51, :, afwImage.LOCAL] = 0.0
        flat[55:56, :, afwImage.LOCAL] = 0.0
        flat[58:59, :, afwImage.LOCAL] = 0.0
        flat[51:60, 51:, afwImage.LOCAL] = 0.0

        mi /= flat

        if display:
            afwDisplay.Display(frame=0).mtv(mi, title=self._testMethodName + ": Raw")

        defectList = algorithms.Defects()
        bbox = lsst.geom.BoxI(lsst.geom.PointI(50, 0), lsst.geom.ExtentI(1, 100))
        defectList.append(algorithms.Defect(bbox))
        bbox = lsst.geom.BoxI(lsst.geom.PointI(55, 0), lsst.geom.ExtentI(1, 100))
        defectList.append(algorithms.Defect(bbox))
        bbox = lsst.geom.BoxI(lsst.geom.PointI(58, 0), lsst.geom.ExtentI(1, 100))
        defectList.append(algorithms.Defect(bbox))
        bbox = lsst.geom.BoxI(lsst.geom.PointI(51, 51), lsst.geom.ExtentI(9, 49))
        defectList.append(algorithms.Defect(bbox))

        psf = algorithms.DoubleGaussianPsf(15, 15, 1./(2*math.sqrt(2*math.log(2))))
        algorithms.interpolateOverDefects(mi, psf, defectList, 50.)

        if display:
            afwDisplay.Display(frame=1).mtv(mi, title=self._testMethodName + ": Interpolated")

        self.assertTrue(np.isfinite(mi.image[56, 51, afwImage.LOCAL]))

    @unittest.skipUnless(afwdataDir, "afwdata not available")
    def testEdge(self):
        """Test that we can interpolate to the edge"""
        mi = afwImage.MaskedImageF(80, 30)

        ima = mi.getImage().getArray()
        #
        # Loop over number of bad columns at left or right edge of image
        #
        for nBadCol in range(0, 20):
            mi.set((0, 0x0, 0))

            np.random.seed(666)
            ima[:] = np.random.uniform(-1, 1, ima.shape)

            defects = []

            if nBadCol > 0:
                #
                # Bad left edge
                #
                ima[:, 0:nBadCol] = 10
                defects.append(lsst.geom.BoxI(lsst.geom.PointI(0, 0),
                                              lsst.geom.ExtentI(nBadCol, mi.getHeight())))
                #
                # With another bad set of columns next to bad left edge
                #
                ima[:, -nBadCol:] = 10
                defects.append(lsst.geom.BoxI(lsst.geom.PointI(mi.getWidth() - nBadCol, 0),
                                              lsst.geom.ExtentI(nBadCol, mi.getHeight())))
                #
                # Bad right edge
                #
                ima[0:10, nBadCol+1:nBadCol+4] = 100
                defects.append(lsst.geom.BoxI(lsst.geom.PointI(nBadCol+1, 0),
                                              lsst.geom.ExtentI(3, 10)))
                #
                # With another bad set of columns next to bad right edge
                #
                ima[0:10, -nBadCol-4:-nBadCol-1] = 100
                defects.append((lsst.geom.BoxI(lsst.geom.PointI(mi.getWidth() - nBadCol - 4, 0),
                                               lsst.geom.ExtentI(3, 10))))
            #
            # Test cases that left and right bad patches nearly (or do) coalesce
            #
            ima[-3:, 0:mi.getWidth()//2-1] = 100
            defects.append(lsst.geom.BoxI(lsst.geom.PointI(0, mi.getHeight() - 3),
                                          lsst.geom.ExtentI(mi.getWidth()//2-1, 1)))

            ima[-3:, mi.getWidth()//2+1:] = 100
            defects.append(lsst.geom.BoxI(lsst.geom.PointI(mi.getWidth()//2 + 1, mi.getHeight() - 3),
                                          lsst.geom.ExtentI(mi.getWidth()//2 - 1, 1)))

            ima[-2:, 0:mi.getWidth()//2] = 100
            defects.append(lsst.geom.BoxI(lsst.geom.PointI(0, mi.getHeight() - 2),
                                          lsst.geom.ExtentI(mi.getWidth()//2, 1)))

            ima[-2:, mi.getWidth()//2+1:] = 100
            defects.append(lsst.geom.BoxI(lsst.geom.PointI(mi.getWidth()//2 + 1, mi.getHeight() - 2),
                                          lsst.geom.ExtentI(mi.getWidth()//2 - 1, 1)))

            ima[-1:, :] = 100
            defects.append(lsst.geom.BoxI(lsst.geom.PointI(0, mi.getHeight() - 1),
                                          lsst.geom.ExtentI(mi.getWidth(), 1)))

            # Test fix for HSC-978: long defect stops one pixel shy of the edge (when nBadCol == 0)
            ima[13, :-1] = 100
            defects.append(lsst.geom.BoxI(lsst.geom.PointI(0, 13), lsst.geom.ExtentI(mi.getWidth() - 1, 1)))
            ima[14, 1:] = 100
            defects.append(lsst.geom.BoxI(lsst.geom.PointI(1, 14), lsst.geom.ExtentI(mi.getWidth() - 1, 1)))

            #
            # Build list of defects to interpolate over
            #
            defectList = algorithms.Defects()

            for bbox in defects:
                defectList.append(algorithms.Defect(bbox))
            #
            # Guess a PSF and do the work
            #
            if display:
                afwDisplay.Display(frame=2).mtv(mi, title=self._testMethodName + ": image")

            psf = algorithms.DoubleGaussianPsf(15, 15, 1./(2*math.sqrt(2*math.log(2))))
            algorithms.interpolateOverDefects(mi, psf, defectList, 0, True)

            if display:
                afwDisplay.Display(frame=3).mtv(mi, title=self._testMethodName + ": image")

            self.assertGreater(np.min(ima), -2)
            self.assertGreater(2, np.max(ima))


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
