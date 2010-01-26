#!/usr/bin/env python
# -*- lsst-python -*-
"""
Tests for ticket 1043 - Photometry fails when no PSF is provided
"""

import lsst.meas.algorithms as measAlgorithms
import lsst.afw.image.imageLib as afwImage

import math
import unittest
import lsst.utils.tests as utilsTests

# math.isnan() available in 2.6, but not 2.5.2
try:
    math.isnan(1.0)
except AttributeError:
    math.isnan = lambda x: x != x

class ticket1043TestCase(unittest.TestCase):

    def setUp(self):
        self.mi         = afwImage.MaskedImageF(100,100)
        self.sincPhotometer = measAlgorithms.createMeasurePhotometry("SINC", 3.0)
        self.naivPhotometer = measAlgorithms.createMeasurePhotometry("NAIVE", 3.0)

    def tearDown(self):
        del self.mi
        del self.sincPhotometer
        del self.naivPhotometer

    def testTicket1043(self):
        self.mi.set(50, 50, (1, 0x0, 1))
        sincPhot = self.sincPhotometer.apply(self.mi, 50.0, 50.0)
        naivPhot = self.naivPhotometer.apply(self.mi, 50.0, 50.0)

        # make sure aperture photometry works
        knownSincApFlux = 1.0111901760101318      
        self.assertEqual(naivPhot.getApFlux(), 1.0)
        self.assertEqual(sincPhot.getApFlux(), knownSincApFlux)

        # make sure psf photometry returns nan
        self.assertTrue(math.isnan(naivPhot.getPsfFlux()))
        self.assertTrue(math.isnan(sincPhot.getPsfFlux()))
        
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(ticket1043TestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)


def run(exit = False):
    """Run the tests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
 
