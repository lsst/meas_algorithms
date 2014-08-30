#!/usr/bin/env python

# 
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

import re
import os
import glob
import math
import pdb                          # we may want to say pdb.set_trace()
import unittest

import eups
import lsst.pex.exceptions as pexExceptions
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.afw.detection as afwDetection
import lsst.meas.algorithms as algorithms
import lsst.utils.tests as utilsTests

import testLib

try:
    type(verbose)
except NameError:
    verbose = 0

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class CentroidTestCase(unittest.TestCase):
    """A test case for centroiding"""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testMeasureCentroid(self):
        """Test that we can instantiate and play with SillyMeasureCentroid"""

        for imageFactory in (afwImage.MaskedImageF,
                             afwImage.MaskedImageD,
                             ):
            im = imageFactory(afwGeom.ExtentI(100, 100))
            exp = afwImage.makeExposure(im)
            control1 = testLib.SillyCentroidControl()
            control1.name = "silly1"
            control1.priority = 0.0
            control1.param = 0
            control2 = testLib.SillyCentroidControl()
            control2.name = "silly2"
            control2.priority = 1.0
            control2.param = 1
            control3 = testLib.SillyCentroidControl()
            control3.name = "silly3"
            control3.priority = 2.0
            control3.param = 2
            schema = afwTable.SourceTable.makeMinimalSchema()
            schema.setVersion(0)
            builder = algorithms.MeasureSourcesBuilder()
            builder.addAlgorithm(control1)
            builder.addAlgorithm(control2)
            builder.addAlgorithm(control3)
            centroider =  builder.build(schema)
            table = afwTable.SourceTable.make(schema)
            source = table.makeRecord()
            x, y = 10, 20
            centroider.apply(source, exp, afwGeom.Point2D(x, y))
            table.defineCentroid(control1.name)
            self.assertEqual(x, source.getX() - 0)
            self.assertEqual(y, source.getY() - 0)
            table.defineCentroid(control2.name)
            self.assertEqual(x, source.getX() - 1)
            self.assertEqual(y, source.getY() - 1)
            table.defineCentroid(control3.name)
            self.assertEqual(x, source.getX() - 2)
            self.assertEqual(y, source.getY() - 2)

    def testMeasureCentroid(self):
        """Test that we can use our silly centroid through the usual Tasks"""
        algorithms.AlgorithmRegistry.register("centroid.silly", testLib.SillyCentroidControl)

        x, y = 10, 20

        im = afwImage.MaskedImageF(afwGeom.ExtentI(512, 512))
        im.set(0)
        arr = im.getImage().getArray()
        arr[y,x] = 1
        exp = afwImage.makeExposure(im)

        schema = afwTable.SourceTable.makeMinimalSchema()
        schema.setVersion(0)

        detConfig = algorithms.SourceDetectionConfig()
        detConfig.thresholdValue = 0.5
        detConfig.thresholdType = "value"
        measConfig = algorithms.SourceMeasurementConfig()
        measConfig.algorithms.names.add("centroid.silly")
        measConfig.slots.centroid = "centroid.silly"
        measConfig.algorithms["centroid.silly"].param = 5
        measConfig.doReplaceWithNoise = False

        det = algorithms.SourceDetectionTask(schema=schema, config=detConfig)
        meas = algorithms.SourceMeasurementTask(schema, config=measConfig)

        table = afwTable.SourceTable.make(schema)
        sources = det.makeSourceCatalog(table, exp, doSmooth=False, sigma=1.0).sources
        self.assertEqual(len(sources), 1)
        meas.run(exp, sources)
        self.assertEqual(len(sources), 1)
        self.assertEqual(sources[0].getY(), y + 5)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(CentroidTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)
 
if __name__ == "__main__":
    run(True)
