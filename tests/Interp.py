#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008-2015 LSST Corporation.
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
# see <http://www.lsstcorp.org/LegalNotices/>.
#

"""
Tests for bad pixel interpolation code

Run with:
   python Interp.py
or
   python
   >>> import Interp; Interp.run()
"""

import os
from math import *
import unittest
import lsst.utils
import math, numpy
import lsst.utils.tests as tests
import lsst.pex.logging as logging
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.display.ds9 as ds9
import lsst.meas.algorithms as algorithms
import lsst.meas.algorithms.defects as defects

try:
    type(verbose)
except NameError:
    verbose = 0
    logging.Trace_setVerbosity("algorithms.Interp", verbose)

try:
    type(display)
except NameError:
    display = False

class interpolationTestCase(unittest.TestCase):
    """A test case for interpolation"""
    def setUp(self):
        self.FWHM = 5
        afwdataDir = lsst.utils.getPackageDir('afwdata')
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
        del self.badPixels

    def testDetection(self):
        """Test Interp algorithms"""

        if display:
            frame = 0
            ds9.mtv(self.mi, frame=frame, title="Original")

        algorithms.interpolateOverDefects(self.mi, self.badPixels)

        if display:
            ds9.mtv(self.mi, frame = frame + 1, title="Interpolated")
            ds9.mtv(self.mi.getVariance(), frame = frame + 2, title="Variance")

    def test818(self):
        """A test case for #818; the full test is in /lsst/DC3root/ticketFiles/818"""

        badPixels = algorithms.DefectListT()
        defects = [((82, 663), 6, 8),
                   ((83, 659), 9, 6),
                   ((85, 660), 10, 11),
                   ((87, 669), 3, 3),
                   ]

        for xy0, width, height in defects:
            x0, y0 = xy0
            bbox = afwGeom.BoxI(afwGeom.PointI(x0, y0), afwGeom.ExtentI(width, height))
            badPixels.push_back(algorithms.Defect(bbox))

        mi = afwImage.MaskedImageF(517, 800)

        algorithms.interpolateOverDefects(mi, badPixels)

    def test1295(self):
        """A test case for #1295 (failure to interpolate over groups of defects"""

        im = afwImage.ImageF(afwGeom.ExtentI(100, 100))
        mi = afwImage.makeMaskedImage(im)
        mi.set(100)
        flat = afwImage.ImageF(im.getDimensions())
        flat.set(1)
        for i in range(100):
            for j in range(100):
                if i == 50 or i == 55 or i == 58:
                    flat.set(i,j,0)
                if i < 60 and i > 50 and j > 50:
                    flat.set(i,j,0)

        mi /= flat

        if display:
            ds9.mtv(mi, frame=0, title="Raw")

        defectList = algorithms.DefectListT()
        bbox = afwGeom.BoxI(afwGeom.PointI(50,0), afwGeom.ExtentI(1,100))
        defectList.append(algorithms.Defect(bbox))
        bbox = afwGeom.BoxI(afwGeom.PointI(55,0), afwGeom.ExtentI(1,100))
        defectList.append(algorithms.Defect(bbox))
        bbox = afwGeom.BoxI(afwGeom.PointI(58,0), afwGeom.ExtentI(1,100))
        defectList.append(algorithms.Defect(bbox))
        bbox = afwGeom.BoxI(afwGeom.PointI(51,51), afwGeom.ExtentI(9,49))
        defectList.append(algorithms.Defect(bbox))

        algorithms.interpolateOverDefects(mi, defectList, 50.)

        if display:
            ds9.mtv(mi, frame=1, title="Interpolated")

        self.assertTrue(numpy.isfinite(mi.getImage().get(56, 51)))

    def testEdge(self):
        """Test that we can interpolate to the edge"""
        mi = afwImage.MaskedImageF(80, 30)

        ima = mi.getImage().getArray()
        #
        # Loop over number of bad columns at left or right edge of image
        #
        for nBadCol in range(0, 20):
            mi.set((0, 0x0, 0))

            numpy.random.seed(666)
            ima[:] = numpy.random.uniform(-1, 1, ima.shape)

            defects = []

            if nBadCol > 0:
                #
                # Bad left edge
                #
                ima[:, 0:nBadCol] = 10
                defects.append(afwGeom.BoxI(afwGeom.PointI(0,0),
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
                defects.append(afwGeom.BoxI(afwGeom.PointI(nBadCol+1,0),
                                            afwGeom.ExtentI(3, 10)))
                #
                # With another bad set of columns next to bad right edge
                #
                ima[0:10, -nBadCol-4:-nBadCol-1] = 100
                defects.append((afwGeom.BoxI(afwGeom.PointI(mi.getWidth() - nBadCol - 4,0),
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
            defectList = algorithms.DefectListT()

            for bbox in defects:
                defectList.append(algorithms.Defect(bbox))
            #
            # Do the work
            #
            if display:
                ds9.mtv(mi, frame=0)

            algorithms.interpolateOverDefects(mi, defectList, 0, True)

            if display:
                ds9.mtv(mi, frame=1)

            self.assertGreater(numpy.min(ima), -2)
            self.assertGreater(2, numpy.max(ima))

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    try:
        lsst.utils.getPackageDir('afwdata')
        suites += unittest.makeSuite(interpolationTestCase)
    except Exception:
        print "Skipping interpolation test case as afwdata isn't set up"
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)


def run(exit = False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
