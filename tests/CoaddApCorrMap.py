#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2014 LSST Corporation.
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
Tests for CoaddApCorr code

Run with:
   python CoaddApCorr.py
"""

import os, sys
from math import *
import numpy
import unittest
import eups
import lsst.utils.tests as utilsTests
import lsst.pex.exceptions as pexExceptions
import lsst.pex.logging as logging

import math
import pdb
import numpy

import lsst.pex.config as pexConfig
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.afw.detection as afwDetection
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.coord as afwCoord
import lsst.pipe.base as pipeBase
import lsst.afw.cameraGeom as cameraGeom
import lsst.meas.algorithms as measAlg

try:
    type(verbose)
except NameError:
    verbose = 0
    display = False


class CoaddApCorrMapTest(unittest.TestCase):

    def test(self):
        """Check that we can create and use a coadd ApCorrMap"""
        coaddBox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.Extent2I(100, 100))
        scale = 5.0e-5 # deg/pix; for CD matrix
        coord = afwCoord.Coord(0.0*afwGeom.degrees, 0.0*afwGeom.degrees)
        center = afwGeom.Point2D(afwGeom.Extent2D(coaddBox.getDimensions())*0.5)
        coaddWcs = afwImage.makeWcs(coord, afwGeom.Point2D(0, 0), scale, 0.0, 0.0, scale)
        schema = afwTable.ExposureTable.makeMinimalSchema()
        weightKey = schema.addField("customweightname", type="D", doc="Coadd weight")
        catalog = afwTable.ExposureCatalog(schema)

        # Non-overlapping
        num = 5
        inputBox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.Extent2I(10, 10))
        pointList = []
        overlapping = []
        for i in range(num):
            value = numpy.array([[1]], dtype=float) # Constant with value = i+1
            apCorrMap = afwImage.ApCorrMap()
            bf = afwMath.ChebyshevBoundedField(inputBox, value*(i + 1))
            apCorrMap.set("only", bf)

            point = afwGeom.Point2D(0, 0) - afwGeom.Extent2D(coaddBox.getDimensions())*(i+0.5)/num
            wcs = afwImage.makeWcs(coord, point, scale, 0.0, 0.0, scale)
            center = afwGeom.Box2D(inputBox).getCenter()
            pointList.append(coaddWcs.skyToPixel(wcs.pixelToSky(center)))
            record = catalog.getTable().makeRecord()
            record.setWcs(wcs)
            record.setBBox(inputBox)
            record.setApCorrMap(apCorrMap)
            record.set(weightKey, i + 1)
            record['id'] = i
            catalog.append(record)

            # An overlapping record
            record = catalog.getTable().makeRecord()
            record.setWcs(wcs)
            record.setBBox(inputBox)
            apCorrMap = afwImage.ApCorrMap()
            bf = afwMath.ChebyshevBoundedField(inputBox, value*(i + 2))
            apCorrMap.set("only", bf)
            record.setApCorrMap(apCorrMap)
            record.set(weightKey, i + 2)
            record['id'] = i + num
            catalog.append(record)

        apCorrMap = measAlg.makeCoaddApCorrMap(catalog, coaddBox, coaddWcs, "customweightname")
        self.assertApCorrMap(apCorrMap, pointList)

        filename = "tests/coaddApCorrMap.fits"
        exposure = afwImage.ExposureF(1, 1)
        exposure.getInfo().setApCorrMap(apCorrMap)
        exposure.writeFits(filename)
        exposure = afwImage.ExposureF(filename)
        self.assertApCorrMap(exposure.getInfo().getApCorrMap(), pointList)
        os.unlink(filename)

    def assertApCorrMap(self, apCorrMap, pointList):
        for i, point in enumerate(pointList):
            weights = [i+1, i+2]
            values = [i+1, i+2]
            expected = sum((w*v for w,v in zip(weights, values)), 0.0) / sum(weights)
            actual = apCorrMap["only"].evaluate(point)
            self.assertEqual(actual, expected)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(CoaddApCorrMapTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the utilsTests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
