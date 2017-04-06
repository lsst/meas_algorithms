#
# LSST Data Management System
#
# Copyright 2008-2017  AURA/LSST.
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
import math
import unittest
import itertools

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.coord as afwCoord
from lsst.meas.algorithms import BinnedWcs
import lsst.utils.tests


class BinnedWcsTest(lsst.utils.tests.TestCase):

    def setUp(self):
        self.scale = (1.0*afwGeom.arcseconds).asDegrees()  # degrees/pixel
        self.wcs = afwImage.makeWcs(afwCoord.Coord(0.0*afwGeom.degrees, 0.0*afwGeom.degrees),
                                    afwGeom.Point2D(0.0, 0.0), self.scale, 0.0, 0.0, self.scale)

    def tearDown(self):
        del self.wcs

    def assertPointEqual(self, p, q, places=6):
        self.assertAlmostEqual(p.getX(), q.getX(), places=places)
        self.assertAlmostEqual(p.getY(), q.getY(), places=places)

    def assertWcs(self, origWcs, binnedWcs, xBin, yBin, x0, y0, size=1000):
        self.assertAnglesAlmostEqual(binnedWcs.pixelScale(), origWcs.pixelScale()*math.sqrt(xBin*yBin))
        for xOrig, yOrig in itertools.product((x0, x0 + size), (y0, y0 + size)):
            orig = afwGeom.Point2D(xOrig, yOrig)
            binned = afwGeom.Point2D(float(xOrig - x0)/xBin, float(yOrig - y0)/yBin)
            coord = origWcs.pixelToSky(orig)

            actual = binnedWcs.skyToPixel(origWcs.pixelToSky(orig))
            self.assertPointEqual(actual, binned)

            actual = origWcs.skyToPixel(binnedWcs.pixelToSky(binned))
            self.assertPointEqual(actual, orig)

            self.assertAlmostEqual(binnedWcs.pixArea(binned), origWcs.pixArea(orig)*xBin*yBin)

            units = afwGeom.arcseconds
            sky = afwGeom.Point2D(coord.getLongitude().asAngularUnits(units),
                                  coord.getLatitude().asAngularUnits(units))
            skyToPixel = binnedWcs.linearizeSkyToPixel(coord, units)
            self.assertPointEqual(skyToPixel(sky), binned)
            pixelToSky = binnedWcs.linearizePixelToSky(coord, units)
            self.assertPointEqual(pixelToSky(binned), sky)

    def testCases(self):
        for xBin, yBin, x0, y0 in [(1, 1, 0, 0),  # Pass-through
                                   (1, 1, 12345, 6789),  # Offset only
                                   (100, 100, 0, 0),  # Binning only
                                   (8, 3, 0, 0),     # Different binnings
                                   (100, 100, 12345, 6789),  # Binning and offset
                                   (4, 7, 9876, 54321),  # Different binnings and offset
                                   ]:
            print("Testing:", xBin, yBin, x0, y0)
            binnedWcs = BinnedWcs(self.wcs, xBin, yBin, afwGeom.Point2I(x0, y0))
            self.assertWcs(self.wcs, binnedWcs, xBin, yBin, x0, y0)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
