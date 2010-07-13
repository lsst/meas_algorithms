#!/usr/bin/env python
# -*- lsst-python -*-
"""
Tests for ticket 1043 - Photometry fails when no PSF is provided
"""

import lsst.pex.policy as pexPolicy
import lsst.meas.algorithms as measAlgorithms
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetection

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
        self.mi = afwImage.MaskedImageF(100, 100)
        self.mi.set(0, 0x0, 1)

        self.measurePhotom = measAlgorithms.MeasurePhotometryF(self.mi)

        for alg in ("NAIVE", "PSF", "SINC",):
            self.measurePhotom.addAlgorithm(alg)

        pol = pexPolicy.Policy(pexPolicy.PolicyString(
            """#<?cfg paf policy?>
            NAIVE.radius: 10.0
            SINC.radius: 3.0
            """
            ))
            
        self.measurePhotom.configure(pol)

    def tearDown(self):
        del self.mi
        del self.measurePhotom

    def testTicket1043(self):
        self.mi.set(50, 50, (1, 0x0, 1))
        peak = afwDetection.Peak(50, 50)

        photom = self.measurePhotom.measure(peak)

        # make sure aperture photometry works
        knownSincApFlux = 1.0111901760101318
        
        self.assertEqual(photom.find("NAIVE").getFlux(), 1.0)
        self.assertEqual(photom.find("SINC").getFlux(),  knownSincApFlux)
        print photom.find("PSF").getFlux()

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
 
