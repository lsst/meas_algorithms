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
from builtins import range
import os
import unittest
import math
import numpy as np

import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.display.ds9 as ds9
import lsst.meas.algorithms as algorithms
import lsst.meas.algorithms.defects as defects
import lsst.utils.tests


try:
    type(display)
except NameError:
    display = False

# Determine if we have afwdata
try:
    afwdataDir = lsst.utils.getPackageDir('afwdata')
except Exception:
    afwdataDir = None


class interpolationTestCase(lsst.utils.tests.TestCase):
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
        self.badPixels = defects.policyToBadRegionList(
            os.path.join(measAlgorithmsDir, "policy", "BadPixels.paf"))

    def tearDown(self):
        del self.mi
        del self.psf
        del self.badPixels

    @unittest.skipUnless(afwdataDir, "afwdata not available")
    def testDetection(self):
        """Test Interp algorithms."""

        if display:
            frame = 0
            ds9.mtv(self.mi, frame=frame, title="Original")

        algorithms.interpolateOverDefects(self.mi, self.psf, self.badPixels)

        if display:
            ds9.mtv(self.mi, frame=frame + 1, title="Interpolated")
            ds9.mtv(self.mi.getVariance(), frame=frame + 2, title="Variance")

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
            bbox = afwGeom.BoxI(afwGeom.PointI(x0, y0), afwGeom.ExtentI(width, height))
            badPixels.append(algorithms.Defect(bbox))

        mi = afwImage.MaskedImageF(517, 800)

        algorithms.interpolateOverDefects(mi, self.psf, badPixels)

    @unittest.skipUnless(afwdataDir, "afwdata not available")
    def test1295(self):
        """A test case for #1295 (failure to interpolate over groups of defects."""
        im = afwImage.ImageF(afwGeom.ExtentI(100, 100))
        mi = afwImage.makeMaskedImage(im)
        mi.set(100)
        flat = afwImage.ImageF(im.getDimensions())
        flat.set(1)
        for i in range(100):
            for j in range(100):
                if i == 50 or i == 55 or i == 58:
                    flat.set(i, j, 0)
                if i < 60 and i > 50 and j > 50:
                    flat.set(i, j, 0)

        mi /= flat

        if display:
            ds9.mtv(mi, frame=0, title="Raw")

        defectList = []
        bbox = afwGeom.BoxI(afwGeom.PointI(50, 0), afwGeom.ExtentI(1, 100))
        defectList.append(algorithms.Defect(bbox))
        bbox = afwGeom.BoxI(afwGeom.PointI(55, 0), afwGeom.ExtentI(1, 100))
        defectList.append(algorithms.Defect(bbox))
        bbox = afwGeom.BoxI(afwGeom.PointI(58, 0), afwGeom.ExtentI(1, 100))
        defectList.append(algorithms.Defect(bbox))
        bbox = afwGeom.BoxI(afwGeom.PointI(51, 51), afwGeom.ExtentI(9, 49))
        defectList.append(algorithms.Defect(bbox))

        psf = algorithms.DoubleGaussianPsf(15, 15, 1./(2*math.sqrt(2*math.log(2))))
        algorithms.interpolateOverDefects(mi, psf, defectList, 50.)

        if display:
            ds9.mtv(mi, frame=1, title="Interpolated")

        self.assertTrue(np.isfinite(mi.getImage().get(56, 51)))

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
                defects.append(afwGeom.BoxI(afwGeom.PointI(0, 0),
                                            afwGeom.ExtentI(nBadCol, mi.getHeight())))
                #
                # With another bad set of columns next to bad left edge
                #
                ima[:, -nBadCol:] = 10
                defects.append(afwGeom.BoxI(afwGeom.PointI(mi.getWidth() - nBadCol, 0),
                                            afwGeom.ExtentI(nBadCol, mi.getHeight())))
                #
                # Bad right edge
                #
                ima[0:10, nBadCol+1:nBadCol+4] = 100
                defects.append(afwGeom.BoxI(afwGeom.PointI(nBadCol+1, 0),
                                            afwGeom.ExtentI(3, 10)))
                #
                # With another bad set of columns next to bad right edge
                #
                ima[0:10, -nBadCol-4:-nBadCol-1] = 100
                defects.append((afwGeom.BoxI(afwGeom.PointI(mi.getWidth() - nBadCol - 4, 0),
                                             afwGeom.ExtentI(3, 10))))
            #
            # Test cases that left and right bad patches nearly (or do) coalesce
            #
            ima[-3:, 0:mi.getWidth()//2-1] = 100
            defects.append(afwGeom.BoxI(afwGeom.PointI(0, mi.getHeight() - 3),
                                        afwGeom.ExtentI(mi.getWidth()//2-1, 1)))

            ima[-3:, mi.getWidth()//2+1:] = 100
            defects.append(afwGeom.BoxI(afwGeom.PointI(mi.getWidth()//2 + 1, mi.getHeight() - 3),
                                        afwGeom.ExtentI(mi.getWidth()//2 - 1, 1)))

            ima[-2:, 0:mi.getWidth()//2] = 100
            defects.append(afwGeom.BoxI(afwGeom.PointI(0, mi.getHeight() - 2),
                                        afwGeom.ExtentI(mi.getWidth()//2, 1)))

            ima[-2:, mi.getWidth()//2+1:] = 100
            defects.append(afwGeom.BoxI(afwGeom.PointI(mi.getWidth()//2 + 1, mi.getHeight() - 2),
                                        afwGeom.ExtentI(mi.getWidth()//2 - 1, 1)))

            ima[-1:, :] = 100
            defects.append(afwGeom.BoxI(afwGeom.PointI(0, mi.getHeight() - 1),
                                        afwGeom.ExtentI(mi.getWidth(), 1)))

            # Test fix for HSC-978: long defect stops one pixel shy of the edge (when nBadCol == 0)
            ima[13, :-1] = 100
            defects.append(afwGeom.BoxI(afwGeom.PointI(0, 13), afwGeom.ExtentI(mi.getWidth() - 1, 1)))
            ima[14, 1:] = 100
            defects.append(afwGeom.BoxI(afwGeom.PointI(1, 14), afwGeom.ExtentI(mi.getWidth() - 1, 1)))

            #
            # Build list of defects to interpolate over
            #
            defectList = []

            for bbox in defects:
                defectList.append(algorithms.Defect(bbox))
            #
            # Guess a PSF and do the work
            #
            if display:
                ds9.mtv(mi, frame=0)

            psf = algorithms.DoubleGaussianPsf(15, 15, 1./(2*math.sqrt(2*math.log(2))))
            algorithms.interpolateOverDefects(mi, psf, defectList, 0, True)

            if display:
                ds9.mtv(mi, frame=1)

            self.assertGreater(np.min(ima), -2)
            self.assertGreater(2, np.max(ima))

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
