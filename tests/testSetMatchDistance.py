#!/usr/bin/env python

# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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
import unittest

import numpy

import lsst.utils.tests as tests
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
from lsst.meas.algorithms import LoadReferenceObjectsTask
from lsst.meas.base import SingleFrameMeasurementTask
from lsst.meas.algorithms import setMatchDistance

class BaseTestCase(unittest.TestCase):
    """A test case for setMatchDistance

    Use involves setting one class attribute:
    * MatchClass: match class, e.g. ReferenceMatch or SourceMatch

    This test is a bit messy because it exercises two templates of MatchClass

    The test creates source and reference object catalogs that intentionally have
    some separation in on-sky coordinates. The reference catalog is set using
    a uniform grid of pixel positions and a simple WCS to compute on-sky coordinates.
    The source catalog is created by using a distored version of the same grid
    of pixel positions, which is converted to on-sky coordinates using the same WCS.
    """
    MatchClass = None

    def setUp(self):
        crval = afwCoord.IcrsCoord(afwGeom.PointD(44., 45.))
        crpix = afwGeom.PointD(0, 0)
        
        arcsecPerPixel = 1/3600.0
        CD11 = arcsecPerPixel
        CD12 = 0
        CD21 = 0
        CD22 = arcsecPerPixel
        
        self.tanWcs = afwImage.makeWcs(crval, crpix, CD11, CD12, CD21, CD22)

        S = 300
        N = 5

        if self.MatchClass == afwTable.ReferenceMatch:
            refSchema = LoadReferenceObjectsTask.makeMinimalSchema(
                filterNameList = ["r"], addFluxSigma=True, addIsPhotometric=True)
            self.refCat = afwTable.SimpleCatalog(refSchema)
        elif self.MatchClass == afwTable.SourceMatch:
            refSchema = afwTable.SourceTable.makeMinimalSchema()
            self.refCat = afwTable.SourceCatalog(refSchema)
        else:
            raise RuntimeError("Unsupported MatchClass=%r" % (self.MatchClass,))
        srcSchema = afwTable.SourceTable.makeMinimalSchema()
        SingleFrameMeasurementTask(schema=srcSchema)
        self.refCoordKey = refSchema["coord"].asKey()
        self.srcCoordKey = srcSchema["coord"].asKey()
        self.srcCentroidKey = afwTable.Point2DKey(srcSchema["slot_Centroid"])
        self.sourceCat = afwTable.SourceCatalog(srcSchema)
        self.origSourceCat = afwTable.SourceCatalog(srcSchema) # undistorted copy
        self.matches = []

        for i in numpy.linspace(0., S, N):
            for j in numpy.linspace(0., S, N):
                src = self.sourceCat.addNew()
                refObj = self.refCat.addNew()

                src.set(self.srcCentroidKey, afwGeom.Point2D(i, j))

                c = self.tanWcs.pixelToSky(afwGeom.Point2D(i, j))
                refObj.setCoord(c)

                self.matches.append(self.MatchClass(refObj, src, 0.0))

    def tearDown(self):
        del self.refCat
        del self.origSourceCat
        del self.sourceCat
        del self.matches
        del self.tanWcs

    def testTrivial(self):
        """Add no distortion"""
        self.doTest("testTrivial", lambda x, y: (x, y))

    def testQuadraticX(self):
        """Add quadratic distortion in x"""
        self.doTest("testQuadraticX", lambda x, y: (x + 1e-4*x**2, y))

    def testRadial(self):
        """Add radial distortion"""
        radialTransform = afwGeom.RadialXYTransform([0, 1.02, 1e-6])
        def radialDistortion(x, y):
            x, y = radialTransform.forwardTransform(afwGeom.Point2D(x, y))
            return (x, y)
        self.doTest("testRadial", radialDistortion)

    def doTest(self, name, func):
        """Apply func(x, y) to each source in self.sourceCat, then set coord, compute and check dist
        """
        for refObj, src, d in self.matches:
            origPos = src.get(self.srcCentroidKey)
            x, y = func(*origPos)
            distortedPos = afwGeom.Point2D(*func(*origPos))
            src.set(self.srcCentroidKey, distortedPos)
            src.set(self.srcCoordKey, self.tanWcs.pixelToSky(distortedPos))

        setMatchDistance(self.matches)
        maxDistErr = afwGeom.Angle(0)
        for refObj, source, distRad in self.matches:
            sourceCoord = source.get(self.srcCoordKey)
            refCoord = refObj.get(self.refCoordKey)
            predDist = sourceCoord.angularSeparation(refCoord)
            distErr = abs(predDist - distRad*afwGeom.radians)
            maxDistErr = max(distErr, maxDistErr)

        self.assertLess(maxDistErr.asArcseconds(), 1e-7)

def makeTestCase(_MatchClass):
    class SetMatchDistanceTestCase(BaseTestCase):
        MatchClass = _MatchClass
    return SetMatchDistanceTestCase

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(makeTestCase(afwTable.ReferenceMatch))
    suites += unittest.makeSuite(makeTestCase(afwTable.SourceMatch))
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
