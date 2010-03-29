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

    def testCentroidClass(self):
        """Test that we can instantiate and play with Centroid, the class"""

        x, xErr = 10, 1
        y, yErr = 20, 2
        covar = -1

        c = algorithms.Centroid(x, y)
        self.assertEqual(x, c.getX())
        self.assertEqual(y, c.getY())

        c = algorithms.Centroid(algorithms.xyAndError(x, xErr), algorithms.xyAndError(y, yErr), covar)
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

    def testInvalidmeasureCentroid(self):
        """Test that we cannot instantiate an unknown measureCentroid"""

        def getInvalid():
            centroider = algorithms.createMeasureCentroid("XXX")

        utilsTests.assertRaisesLsstCpp(self, pexExceptions.NotFoundException, getInvalid)

    def do_testmeasureCentroid(self, centroiderType):
        """Test that we can instantiate and play with a measureCentroid"""

        for imageFactory in (afwImage.ImageF, afwImage.ImageI):

            im = imageFactory(100, 100)
            if imageFactory == afwImage.ImageF: # only ImageF supports the old no-Image API; this is
                                                # set in the centroid.i swig interface
                centroider = algorithms.createMeasureCentroid(centroiderType, im)
                centroider = algorithms.createMeasureCentroid(centroiderType)
            else:
                centroider = algorithms.createMeasureCentroid(centroiderType, im)

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

    def testNaivemeasureCentroid(self):
        """Test that we can instantiate and play with NaivemeasureCentroid"""

        self.do_testmeasureCentroid("NAIVE")

    def testSDSSmeasureCentroid(self):
        """Test that we can instantiate and play with SDSSmeasureCentroid"""

        self.do_testmeasureCentroid("SDSS")

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
            x0, y0, x1, y1 = bbox.getX0(), bbox.getY0(), bbox.getX1(), bbox.getY1()
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

    def tearDown(self):
        pass

    def testMeasureCentroid(self):
        """Test that we can instantiate and play with a measureCentroid"""
        centroider = algorithms.createMeasureCentroid("GAUSSIAN")

        ID = 1
        for foot in self.ds.getFootprints():
            bbox = foot.getBBox()
            xc = (bbox.getX0() + bbox.getX1())//2
            yc = (bbox.getY0() + bbox.getY1())//2

            c = centroider.apply(self.mi.getImage(), xc, yc, None, 0)

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
