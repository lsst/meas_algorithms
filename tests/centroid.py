#!/usr/bin/env python
import re
import os
import glob
import math
import pdb                          # we may want to say pdb.set_trace()
import unittest

import eups
import lsst.pex.exceptions as pexExceptions
import lsst.afw.image as afwImage
import lsst.meas.algorithms as algorithms; centroid = algorithms
import lsst.utils.tests as utilsTests
import lsst.afw.detection.detectionLib as detect

try:
    type(verbose)
except NameError:
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

    def testCentroidClass(self):
        """Test that we can instantiate and play with Centroid, the class"""

        x, xErr = 10, 1
        y, yErr = 20, 2
        covar = -1

        c = centroid.Centroid(x, y)
        self.assertEqual(x, c.getX())
        self.assertEqual(y, c.getY())

        c = centroid.Centroid(centroid.xyAndError(x, xErr), centroid.xyAndError(y, yErr), covar)
        self.assertEqual((x, xErr), c.getX(1))
        self.assertEqual((y, yErr), c.getY(1))
        self.assertEqual(covar, c.getCovar())

        tmp = 1234
        c.setX(tmp); self.assertEqual(c.getX(), tmp); tmp += 0.5
        c.setX((0, tmp)); self.assertEqual(c.getX(0)[1], tmp); tmp += 0.5
        c.setXErr(tmp); self.assertEqual(c.getXErr(), tmp); tmp += 0.5
        c.setY(tmp); self.assertEqual(c.getY(), tmp); tmp += 0.5
        c.setY((0, tmp)); self.assertEqual(c.getY(0)[1], tmp); tmp += 0.5
        c.setYErr(tmp); self.assertEqual(c.getYErr(), tmp); tmp += 0.5
        c.setCovar(tmp); self.assertEqual(c.getCovar(), tmp); tmp += 0.5

    def testInvalidCentroider(self):
        """Test that we cannot instantiate an unknown Centroider"""

        def getInvalid():
            centroider = centroid.createCentroider("XXX")

        utilsTests.assertRaisesLsstCpp(self, pexExceptions.NotFoundException, getInvalid)

    def do_testCentroider(self, centroiderType):
        """Test that we can instantiate and play with a Centroider"""
        centroider = centroid.createCentroider(centroiderType)

        im = afwImage.ImageF(100, 100)

        #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        bkgd = 10; im.set(bkgd)
        im.set(10, 20, 1010)
        x, y = 10, 20
        c = centroider.apply(im, int(x), int(y), None, bkgd)
        self.assertEqual(x, c.getX())
        self.assertEqual(y, c.getY())

        #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        bkgd = 10; im.set(bkgd)
        im.set(10, 20, 1010)
        im.set(10, 21, 1010)
        im.set(11, 20, 1010)
        im.set(11, 21, 1010)

        x, y = 10.5, 20.5
        c = centroider.apply(im, int(x), int(y), None, bkgd)

        self.assertEqual(x, c.getX())
        self.assertEqual(y, c.getY())

        #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        def centroidEmptySky():
            centroider.apply(im, int(x), int(y))

        im.set(0)
        utilsTests.assertRaisesLsstCpp(self, pexExceptions.RuntimeErrorException, centroidEmptySky)

    def testNaiveCentroider(self):
        """Test that we can instantiate and play with NaiveCentroider"""

        self.do_testCentroider("NAIVE")

    def testSDSSCentroider(self):
        """Test that we can instantiate and play with SDSSCentroider"""

        self.do_testCentroider("SDSS")

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
