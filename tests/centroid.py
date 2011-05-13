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
import lsst.pex.policy as pexPolicy
import lsst.pex.exceptions as pexExceptions
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
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

    def testInvalidMeasureCentroid(self):
        """Test that we cannot instantiate an unknown measureCentroid"""

        def getInvalid():
            centroider = algorithms.makeMeasureAstrometry(None)
            centroider.addAlgorithm("XXX")

        utilsTests.assertRaisesLsstCpp(self, pexExceptions.NotFoundException, getInvalid)

    def do_testAstrometry(self, algorithmName):
        """Test that we can instantiate and play with a centroiding algorithms"""

        for imageFactory in (afwImage.MaskedImageF,
                             afwImage.MaskedImageI,
                             ):

            im = imageFactory(afwGeom.ExtentI(100, 100))

            centroider = algorithms.makeMeasureAstrometry(afwImage.makeExposure(im))
            centroider.addAlgorithm(algorithmName)

            bkgd = 10
            centroider.configure(pexPolicy.Policy(pexPolicy.PolicyString("NAIVE.background: %f" % bkgd)))

            #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

            im.set(bkgd)
            im.set(10, 20, (1010,))
            x, y = 10, 20
            c = centroider.measure(afwDetection.Peak(x, y)).find(algorithmName)
            self.assertEqual(x, c.getX())
            self.assertEqual(y, c.getY())

            #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

            im.set(bkgd)
            im.set(10, 20, (1010,))
            im.set(10, 21, (1010,))
            im.set(11, 20, (1010,))
            im.set(11, 21, (1010,))

            x, y = 10.5, 20.5
            c = centroider.measure(afwDetection.Peak(x, y)).find(algorithmName)

            self.assertEqual(x, c.getX())
            self.assertEqual(y, c.getY())

    def testGaussianMeasureCentroid(self):
        """Test that we can instantiate and play with GAUSSIAN centroids"""

        self.do_testAstrometry("GAUSSIAN")

    def testNaiveMeasureCentroid(self):
        """Test that we can instantiate and play with NAIVE centroids"""

        self.do_testAstrometry("NAIVE")

    def testSdssMeasureCentroid(self):
        """Test that we can instantiate and play with SDSS centroids"""

        self.do_testAstrometry("SDSS")

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class MonetTestCase(unittest.TestCase):
    """A test case for centroiding using Dave Monet's 2-D Gaussian fitter"""

    def setUp(self):
	im = afwImage.ImageF(self.monetFile("small.fits"))
        self.mi = afwImage.MaskedImageF(im, afwImage.MaskU(im.getDimensions()),
                                        afwImage.ImageF(im.getDimensions()));
        self.ds = afwDetection.makeFootprintSet(self.mi, afwDetection.Threshold(100))

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

        self.readTruth(self.monetFile("positions.dat-original"))
        self.ssMeasured = afwDetection.SourceSet()

    def tearDown(self):
        del self.mi
        del self.ds

    def monetFile(self, file):
        """Return a Monet file used for regression testing"""
        return os.path.join(eups.productDir("meas_algorithms"), "tests", "Monet", file)

    def readTruth(self, filename):
        """Read Dave Monet's truth table"""
        self.ssTruth = afwDetection.SourceSet()
        for line in open(filename).readlines():
            if re.search(r"^\s*#", line):
                continue
            status, ID, xSex, xDGM, ySex, yDGM, sky = [float(el) for el in line.split()]

            s = afwDetection.Source()
            s.setId(int(ID))
            s.setXAstrom(xDGM)
            s.setYAstrom(yDGM)
        
            self.ssTruth.append(s)

    def testMeasureCentroid(self):
        """Test that we can instantiate and play with a measureCentroid"""
 
        algorithmName = "GAUSSIAN"
        centroider = algorithms.makeMeasureAstrometry(afwImage.makeExposure(self.mi))
        centroider.addAlgorithm(algorithmName)

        ID = 1
        for foot in self.ds.getFootprints():
            bbox = foot.getBBox()
            xc = (bbox.getMinX() + bbox.getMaxX())//2
            yc = (bbox.getMinY() + bbox.getMaxY())//2

            c = centroider.measure(afwDetection.Peak(xc, yc)).find(algorithmName)

            s = afwDetection.Source()
            s.setId(ID); ID += 1
            s.setXAstrom(c.getX())
            s.setYAstrom(c.getY())

            self.ssMeasured.append(s)

            if display:
                ds9.dot("x", c.getX(), c.getY(), ctype=ds9.GREEN)
        #
        # OK, we've measured all the sources.  Compare positions with Dave Monet's values
        #
        mat = afwDetection.matchXy(self.ssTruth, self.ssMeasured, 1.0)
        #self.assertEqual(ID, len(mat))  # we matched all the input sources

        eps = 6e-6                      # offset in pixels between measured centroid and the Truth
        for match in mat:
            dx = match[0].getXAstrom() - match[1].getXAstrom()
            dy = match[0].getYAstrom() - match[1].getYAstrom()
            
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
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the tests"""
    utilsTests.run(suite(), exit)
 
if __name__ == "__main__":
    run(True)
