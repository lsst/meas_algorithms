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
                             afwImage.MaskedImageI,
                             ):
            im = imageFactory(100, 100)

            centroider =  algorithms.makeMeasureAstrometry(afwImage.makeExposure(im))
            centroider.addAlgorithm("SILLY")
            
            x, y = 10, 20
            c = centroider.measure(afwDetection.Peak(x, y)).find("SILLY")
            self.assertEqual(x, c.getX() - 1)
            self.assertEqual(y, c.getY() - 1)

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
