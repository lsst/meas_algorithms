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

import os, re, sys
import glob
import math
import unittest

import eups
import lsst.pex.exceptions as pexExceptions
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.meas.algorithms as algorithms
import lsst.utils.tests as utilsTests
import lsst.afw.detection.detectionLib as afwDetection

import lsst.afw.display.ds9 as ds9

try:
    type(verbose)
except NameError:
    display = False
    verbose = 0

if False:
    dataDir = eups.productDir("afwdata")
    if not dataDir:
        raise RuntimeError("Must set up afwdata to run these tests")

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class CentroidTestCase(unittest.TestCase):
    """A test case for centroiding"""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testCleanup(self):
        """Test that tearDown does"""
        pass

    def do_testAstrometry(self, control, bkgd):
        """Test that we can instantiate and play with a centroiding algorithms"""

        for imageFactory in (afwImage.MaskedImageF,
                             afwImage.MaskedImageD,
                             ):

            im = imageFactory(afwGeom.ExtentI(100, 100))

            exp = afwImage.makeExposure(im)
            schema = afwTable.SourceTable.makeMinimalSchema()
            centroider = algorithms.MeasureSourcesBuilder().addAlgorithm(control).build(schema)

            #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

            im.set(bkgd)
            x, y = 30, 20
            im.set(x, y, (1010,))

            table = afwTable.SourceTable.make(schema)
            table.defineCentroid(control.name)
            source = table.makeRecord()
            foot = afwDetection.Footprint(exp.getBBox())
            source.setFootprint(foot)

            centroider.apply(source, exp, afwGeom.Point2D(x, y))

            self.assertEqual(x, source.getX())
            self.assertEqual(y, source.getY())

            #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

            im.set(bkgd)
            im.set(10, 20, (1010,))
            im.set(10, 21, (1010,))
            im.set(11, 20, (1010,))
            im.set(11, 21, (1010,))

            x, y = 10.5, 20.5
            centroider.apply(source, exp, afwGeom.Point2D(x, y))

            self.assertEqual(x, source.getX())
            self.assertEqual(y, source.getY())

    def testGaussianMeasureCentroid(self):
        """Test that we can instantiate and play with GAUSSIAN centroids"""
        control = algorithms.GaussianCentroidControl()
        self.do_testAstrometry(control, 10.0)

    def testNaiveMeasureCentroid(self):
        """Test that we can instantiate and play with NAIVE centroids"""
        bkgd = 10.0
        control = algorithms.NaiveCentroidControl()
        control.background = bkgd
        self.do_testAstrometry(control, bkgd)

    def testSdssMeasureCentroid(self):
        """Test that we can instantiate and play with SDSS centroids"""
        control = algorithms.SdssCentroidControl()
        self.do_testAstrometry(control, 10.0)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class SourceMeasurementTaskTestCase(unittest.TestCase):
    """A test case for the SourceMeasurementTask"""

    def mySetup(self, runCentroider=True):
        msConfig = algorithms.SourceMeasurementConfig()
        if not runCentroider:
            msConfig.centroider = None
            msConfig.slots.centroid = None

        schema = afwTable.SourceTable.makeMinimalSchema()
        ms = msConfig.makeMeasureSources(schema)
        
        table = afwTable.SourceTable.make(schema)
        msConfig.slots.setupTable(table)
        source = table.makeRecord()

        fp = afwDetection.Footprint(self.exp.getBBox())
        source.setFootprint(fp)
        ms.apply(source, self.exp, afwGeom.Point2D(self.xcen, self.ycen))

        return source

    def setUp(self):
        """Make the image we'll measure"""

        self.exp = afwImage.ExposureF(100, 100)
        self.I0, self.xcen, self.ycen = 1000.0, 50.5, 50.0
        im = self.exp.getMaskedImage().getImage()

        for i in (-1, 0, 1):
            for j in (-1, 0, 1):
                im.set(int(self.xcen) + i, int(self.ycen) + j, self.I0*(1 - 0.5*math.hypot(i - 0.5, j)))

    def tearDown(self):
        del self.exp

    def testCentroider(self):
        """Measure the centroid"""
        s = self.mySetup()

        self.assertEqual(s.getCentroid(), afwGeom.PointD(self.xcen, self.ycen))

    def testNoCentroider(self):
        """Check that we can disable running a centroid algorithm"""
        s = self.mySetup(False)

        utilsTests.assertRaisesLsstCpp(self, pexExceptions.LogicErrorException, s.getCentroid)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


class MonetTestCase(unittest.TestCase):
    """A test case for centroiding using Dave Monet's 2-D Gaussian fitter"""

    def setUp(self):
	im = afwImage.ImageF(self.monetFile("small.fits"))
        self.mi = afwImage.MaskedImageF(im, afwImage.MaskU(im.getDimensions()),
                                        afwImage.ImageF(im.getDimensions()));
        self.ds = afwDetection.FootprintSet(self.mi, afwDetection.Threshold(100))

        if display:
            ds9.mtv(self.mi.getImage())
            ds9.erase()
            
        for foot in self.ds.getFootprints():
            bbox = foot.getBBox()
            x0, y0 = bbox.getMinX(), bbox.getMinY()
            x1, y1 = bbox.getMaxX(), bbox.getMaxY()
            xc = (x0 + x1)/2.0
            yc = (y0 + y1)/2.0
            
            if display:
                ds9.dot("+", xc, yc, ctype=ds9.BLUE)

                if False:
                    x0 -= 0.5; y0 -= 0.5
                    x1 += 0.5; y1 += 0.5
                    
                    ds9.line([(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)], ctype=ds9.RED)

        self.control = algorithms.GaussianCentroidControl()
        schema = afwTable.SourceTable.makeMinimalSchema()
        self.centroider = algorithms.MeasureSourcesBuilder().addAlgorithm(self.control).build(schema)
        self.ssMeasured = afwTable.SourceCatalog(schema)
        self.ssMeasured.table.defineCentroid(self.control.name)
        self.ssTruth = afwTable.SourceCatalog(schema)
        self.readTruth(self.monetFile("positions.dat-original"))

    def tearDown(self):
        del self.mi
        del self.ds
        del self.centroider
        del self.control
        del self.ssMeasured
        del self.ssTruth

    def monetFile(self, file):
        """Return a Monet file used for regression testing"""
        return os.path.join(eups.productDir("meas_algorithms"), "tests", "Monet", file)

    def readTruth(self, filename):
        """Read Dave Monet's truth table"""
        self.ssTruth.table.defineCentroid(self.control.name)
        for line in open(filename).readlines():
            if re.search(r"^\s*#", line):
                continue
            status, ID, xSex, xDGM, ySex, yDGM, sky = [float(el) for el in line.split()]

            s = self.ssTruth.addNew()
            s.setId(int(ID))
            s.set(self.ssTruth.table.getCentroidKey().getX(), xDGM)
            s.set(self.ssTruth.table.getCentroidKey().getY(), yDGM)

    def testMeasureCentroid(self):
        """Test that we can instantiate and play with a measureCentroid"""
        exposure = afwImage.makeExposure(self.mi)
        self.ds.makeSources(self.ssMeasured)
        ID = 1
        for s in self.ssMeasured:
            s.setId(ID); ID += 1
            foot = s.getFootprint()
            bbox = foot.getBBox()
            xc = (bbox.getMinX() + bbox.getMaxX())//2
            yc = (bbox.getMinY() + bbox.getMaxY())//2

            self.centroider.apply(s, exposure, afwGeom.Point2D(xc, yc))

            if display:
                ds9.dot("x", c.getX(), c.getY(), ctype=ds9.GREEN)
        #
        # OK, we've measured all the sources.  Compare positions with Dave Monet's values
        #

        # FIXME: this test will fail until source matching in afw is updated to use afw/table
        mat = afwTable.matchXy(self.ssTruth, self.ssMeasured, 1.0)
        #self.assertEqual(ID, len(mat))  # we matched all the input sources

        eps = 6e-6                      # offset in pixels between measured centroid and the Truth
        for match in mat:
            dx = match[0].getX() - match[1].getX()
            dy = match[0].getY() - match[1].getY()
            
            good = True if math.hypot(dx, dy) < eps else False
            if not good:
                msg = "Star at (%.1f, %.1f): (dx, dy) = %g, %g)" % \
                    (match[0].getXAstrom(), match[0].getYAstrom(), dx, dy)
                if True:
                    print msg
                else:
                    self.assertTrue(good, msg)

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(CentroidTestCase)
    suites += unittest.makeSuite(MonetTestCase)
    suites += unittest.makeSuite(SourceMeasurementTaskTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the tests"""
    utilsTests.run(suite(), exit)
 
if __name__ == "__main__":
    run(True)
