#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2018 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as
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
from __future__ import absolute_import, division, print_function

import numpy as np
import unittest

import lsst.utils.tests
from lsst.pex.exceptions import InvalidParameterError
from lsst.afw.geom import Point2D, Extent2D, Point2I, Box2D, Box2I, degrees
from lsst.afw.coord import IcrsCoord
from lsst.afw.image import TransmissionCurve, makeWcs
from lsst.afw.geom.polygon import Polygon
from lsst.afw.table import ExposureTable, ExposureCatalog
from lsst.meas.algorithms import makeCoaddTransmissionCurve
from lsst.meas.algorithms.testUtils import makeRandomTransmissionCurve


class CoaddBoundedFieldTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        # Test geometry:
        #
        # -100,99                99,99
        #     +--------------------+
        #     |AAAAAAAAAACCCCCDDDDD|    A == only in epoch A
        #     |AAAAAAAAAACCCCCDDDDD|    B == only in epoch B
        #     |AAAAAAAAAACCCCCDDDDD|    C == in both epoch A and epoch B
        #     |AAAAAAAAAACCCCCDDDDD|    D == in epoch A; in B's bbox but outside its ValidPolygon
        #     |AAAAAAAAAACCCCCDDDDD|
        #     |          BBBBBBBBBB|    All WCSs have the same CRVAL and CD.
        #     |          BBBBBBBBBB|
        #     |          BBBBBBBBBB|    Coadd has CRPIX=(0, 0)
        #     |          BBBBBBBBBB|    Epoch A has CRPIX=(0, -50)
        #     |          BBBBBBBBBB|    Epoch B has CRPIX=(-50, 0)
        #     +--------------------+
        # -100,-100             99,-100
        #
        self.rng = np.random.RandomState(50)
        crval = IcrsCoord(45.0*degrees, 45.0*degrees)
        cd = (5E-5, 0.0, 0.0, 5E-5)
        self.wcsCoadd = makeWcs(crval, Point2D(0.0, 0.0), *cd)
        self.wcsA = makeWcs(crval, Point2D(0.0, -50.0), *cd)
        self.wcsB = makeWcs(crval, Point2D(-50.0, 0.0), *cd)
        self.bboxCoadd = Box2I(Point2I(-100, -100), Point2I(99, 99))
        self.bboxA = Box2I(Point2I(-100, -50), Point2I(99, 49))
        self.bboxB = Box2I(Point2I(-50, -100), Point2I(49, 99))
        self.polygonA = None
        polygonD = Polygon(Box2D(Box2I(Point2I(0, 0), Point2I(49, 99))))
        self.polygonB, = polygonD.symDifference(Polygon(Box2D(self.bboxB)))
        self.curveA = makeRandomTransmissionCurve(self.rng)
        self.curveB = makeRandomTransmissionCurve(self.rng)
        self.weightA = 0.6
        self.weightB = 0.2
        schema = ExposureTable.makeMinimalSchema()
        weightKey = schema.addField("weight", type=float, doc="relative weight of image in Coadd")
        catalog = ExposureCatalog(schema)
        recordA = catalog.addNew()
        recordA[weightKey] = self.weightA
        recordA.setWcs(self.wcsA)
        recordA.setValidPolygon(self.polygonA)
        recordA.setBBox(self.bboxA)
        recordA.setTransmissionCurve(self.curveA)
        recordB = catalog.addNew()
        recordB[weightKey] = self.weightB
        recordB.setWcs(self.wcsB)
        recordB.setValidPolygon(self.polygonB)
        recordB.setBBox(self.bboxB)
        recordB.setTransmissionCurve(self.curveB)
        self.curveCoadd = makeCoaddTransmissionCurve(self.wcsCoadd, catalog)

    def tearDown(self):
        del self.wcsCoadd
        del self.wcsA
        del self.wcsB
        del self.curveCoadd

    def makeRandomPoint(self, *args, **kwds):
        """Draw a random Point2D within a Box2I.

        All arguments are forwarded directly to the Box2I constructor, allowing
        the caller to pass a fully-constructed Box2I, a (Point2I, Point2I) pair,
        or a (Point2I, Extent2I) pair.
        """
        bboxD = Box2D(Box2I(*args, **kwds))
        return bboxD.getMin() + Extent2D(bboxD.getWidth()*self.rng.rand(),
                                         bboxD.getHeight()*self.rng.rand())

    def testSampleAt(self):
        """Test the behavior of TransmissionCurve.sampleAt on the subclass
        returned by makeCoaddTransmissionCurve.
        """
        wavelengths = np.linspace(4000, 7000, 200)

        # Points in coadd coordinates in each of the distinct regions
        point0 = self.makeRandomPoint(Point2I(-100, -100), Point2I(-1, -1))
        pointA = self.makeRandomPoint(Point2I(-100, 0), Point2I(-1, 99))
        pointB = self.makeRandomPoint(Point2I(0, -100), Point2I(99, -1))
        pointC = self.makeRandomPoint(Point2I(0, 0), Point2I(49, 99))
        pointD = self.makeRandomPoint(Point2I(50, 0), Point2I(99, 99))
        points = [point0, pointA, pointB, pointC, pointD]

        # The same points, in sky coordinates
        coords = [self.wcsCoadd.pixelToSky(point) for point in points]

        # The same points, in Epoch A's coordinates
        point0A, pointAA, pointBA, pointCA, pointDA = [self.wcsA.skyToPixel(coord) for coord in coords]
        self.assertFalse(Box2D(self.bboxA).contains(point0A))
        self.assertTrue(Box2D(self.bboxA).contains(pointAA))
        self.assertFalse(Box2D(self.bboxA).contains(pointBA))
        self.assertTrue(Box2D(self.bboxA).contains(pointCA))
        self.assertTrue(Box2D(self.bboxA).contains(pointDA))

        # The same points, in Epoch B's coordinates
        point0B, pointAB, pointBB, pointCB, pointDB = [self.wcsB.skyToPixel(coord) for coord in coords]
        self.assertFalse(Box2D(self.bboxB).contains(point0B))
        self.assertFalse(Box2D(self.bboxB).contains(pointAB))
        self.assertTrue(Box2D(self.bboxB).contains(pointBB))
        self.assertTrue(Box2D(self.bboxB).contains(pointCB))
        self.assertTrue(Box2D(self.bboxB).contains(pointDB))
        self.assertTrue(self.polygonB.contains(pointBB))
        self.assertTrue(self.polygonB.contains(pointCB))
        self.assertFalse(self.polygonB.contains(pointDB))

        # Test that we can't compute throughputs in region 0 (where there are no inputs)
        self.assertRaises(InvalidParameterError, self.curveCoadd.sampleAt, point0, wavelengths)

        # Test throughputs in region A (only Epoch A contributes)
        throughputA1 = self.curveCoadd.sampleAt(pointA, wavelengths)
        throughputA2 = self.curveA.sampleAt(pointAA, wavelengths)
        self.assertFloatsAlmostEqual(throughputA1, throughputA2)

        # Test throughputs in region B (only Epoch B contributes)
        throughputB1 = self.curveCoadd.sampleAt(pointB, wavelengths)
        throughputB2 = self.curveB.sampleAt(pointBB, wavelengths)
        self.assertFloatsAlmostEqual(throughputB1, throughputB2)

        # Test throughputs in region C (both epochs contribute)
        throughputC1 = self.curveCoadd.sampleAt(pointC, wavelengths)
        throughputC2 = self.curveA.sampleAt(pointCA, wavelengths)
        throughputC3 = self.curveB.sampleAt(pointCB, wavelengths)
        self.assertFloatsAlmostEqual(throughputC1, throughputC2*0.75 + throughputC3*0.25)

        # Test throughputs in region D (only Epoch A contributes)
        throughputD1 = self.curveCoadd.sampleAt(pointD, wavelengths)
        throughputD2 = self.curveA.sampleAt(pointDA, wavelengths)
        self.assertFloatsAlmostEqual(throughputD1, throughputD2)

    def testPersistence(self):
        wavelengths = np.linspace(4000, 7000, 200)
        with lsst.utils.tests.getTempFilePath(".fits") as filename:
            self.curveCoadd.writeFits(filename)
            roundtripped = TransmissionCurve.readFits(filename)
        for i in range(10):
            point = self.makeRandomPoint(self.bboxCoadd)
            try:
                throughput1 = self.curveCoadd.sampleAt(point, wavelengths)
            except InvalidParameterError:
                self.assertRaises(InvalidParameterError, roundtripped.sampleAt, point, wavelengths)
            else:
                throughput2 = roundtripped.sampleAt(point, wavelengths)
                self.assertFloatsAlmostEqual(throughput1, throughput2)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
