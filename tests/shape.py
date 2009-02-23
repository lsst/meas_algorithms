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
import lsst.afw.math as afwMath
import lsst.meas.algorithms as algorithms
import lsst.utils.tests as utilsTests
import lsst.afw.detection.detectionLib as detect

import lsst.afw.display.ds9 as ds9

try:
    type(verbose)
except NameError:
    verbose = 0

if False:
    dataDir = eups.productDir("afwdata")
    if not dataDir:
	raise RuntimeError("Must set up afwdata to run these tests")

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class ShapeTestCase(unittest.TestCase):
    """A test case for centroiding"""

    def setUp(self):
	pass

    def tearDown(self):
        pass

    def testCleanup(self):
        """Test that tearDown does"""
        pass

    def testShapeClass(self):
        """Test that we can instantiate and play with Shape, the class"""

        m0 = 1000
        mxx, mxy, myy = 2.0, 1.0, 2.0

        e1 = (mxx - myy)/(mxx + myy)
        e2 = 2*mxy/(mxx + myy)
        rms = math.sqrt(0.5*(mxx + myy))

        c = algorithms.Shape(m0, mxx, mxy, myy)
        self.assertEqual(m0, c.getM0())
        self.assertEqual(mxx, c.getMxx())
        self.assertEqual(mxy, c.getMxy())
        self.assertEqual(myy, c.getMyy())

        self.assertEqual(e1, c.getE1())
        self.assertEqual(e2, c.getE2())
        self.assertEqual(rms, c.getRms())

    def testInvalidmeasureShape(self):
        """Test that we cannot instantiate an unknown measureShape"""

        def getInvalid():
            shapeFinder = algorithms.createMeasureShape("XXX")

        utilsTests.assertRaisesLsstCpp(self, pexExceptions.NotFoundException, getInvalid)

    def do_testmeasureShape(self, shapeFinderType):
        """Test that we can instantiate and play with a measureShape"""
        shapeFinder = algorithms.createMeasureShape(shapeFinderType)

        im = afwImage.ImageF(100, 100)

        #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        bkgd = 100; im.set(bkgd)
        x, y = 30, 40
        im.set(x, y, 1000 + bkgd)

        #
        # Add a Gaussian to the image
        #
        sigma_xx, sigma_xy, sigma_yy, ksize = math.pow(1.5, 2), 0, math.pow(2.5, 2), 15

        if True:
            k = afwMath.AnalyticKernel(ksize, ksize,
                                       afwMath.GaussianFunction2D(math.sqrt(sigma_xx), math.sqrt(sigma_yy)))
            cim = im.Factory(im.getDimensions())
            afwMath.convolve(cim, im, k, True)
            im = cim
        else:
            for dx in range(-ksize/2, ksize/2 + 1):
                for dy in range(-ksize/2, ksize/2 + 1):
                    I = 1000*math.exp(-0.5*(dx*dx/sigma_xx + dy*dy/sigma_yy))
                    im.set(x + dx, y + dy, bkgd + I)

        msk = afwImage.MaskU(im.getDimensions()); msk.set(0)
        var = afwImage.ImageF(im.getDimensions()); var.set(10)
        im = afwImage.MaskedImageF(im, msk, var)
        del msk; del var
            
        if False:
            ds9.mtv(im)

        s = shapeFinder.apply(im, int(x), int(y), None, bkgd)
        self.assertAlmostEqual(x, s.getCentroid().getX(), 6)
        self.assertAlmostEqual(y, s.getCentroid().getY(), 6)

        print "M_xx:  %.5f %.5f" % (s.getMxx(), sigma_xx)
        self.assertTrue(abs(s.getMxx() - sigma_xx) < 1e-3*(1 + sigma_xx))
        print "M_xy:  %.5f %.5f" % (s.getMxy(), sigma_xy)
        self.assertTrue(abs(s.getMxy() - sigma_xy) < 1e-3*(1 + sigma_xy))
        print "M_yy:  %.5f %.5f" % (s.getMyy(), sigma_yy)
        self.assertTrue(abs(s.getMyy() - sigma_yy) < 1e-3*(1 + sigma_yy))

        self.assertTrue(abs((sigma_xx - sigma_yy)/(sigma_xx + sigma_yy) - s.getE1()) < 1e-3)
        self.assertTrue(abs(2*sigma_xy/(sigma_xx + sigma_yy) - s.getE2()) < 1e-3)
        self.assertTrue(abs(math.sqrt(0.5*(sigma_xx + sigma_yy)) - s.getRms()) < 1e-3)

        #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        bkgd = 10; im.set(bkgd)
        im.set(10, 20, 1010)
        im.set(10, 21, 1010)
        im.set(11, 20, 1010)
        im.set(11, 21, 1010)

        x, y = 10.5, 20.5
        c = shapeFinder.apply(im, int(x), int(y), None, bkgd)

        if False:                       # these are left over from centroid.py
            self.assertEqual(x, c.getE1())
            self.assertEqual(y, c.getE2())

        #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        def centroidEmptySky():
            shapeFinder.apply(im, int(x), int(y))

        im.set(0)
        if False:                       # these are left over from centroid.py
            utilsTests.assertRaisesLsstCpp(self, pexExceptions.RuntimeErrorException, centroidEmptySky)

    def testSDSSmeasureShape(self):
        """Test that we can instantiate and play with SDSSmeasureShape"""

        self.do_testmeasureShape("SDSS")

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(ShapeTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)
 
if __name__ == "__main__":
    run(True)
