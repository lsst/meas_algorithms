#!/usr/bin/env python
"""
Tests for measuring things

Run with:
   python MeasureSources.py
or
   python
   >>> import MeasureSources; MeasureSources.run()
"""

import os, sys, unittest
from math import *
import lsst.utils.tests as tests
import lsst.pex.exceptions as pexExceptions
import lsst.pex.logging as pexLogging
import lsst.pex.policy as pexPolicy
import lsst.afw.detection as afwDetection
import lsst.afw.image as afwImage
import lsst.meas.algorithms as measAlg

try:
    type(verbose)
except NameError:
    verbose = 0
pexLogging.Trace_setVerbosity("afwDetection.Measure", verbose)

try:
    type(display)
except NameError:
    display = False

import lsst.afw.display.ds9 as ds9

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class MeasureSourcesTestCase(unittest.TestCase):
    """A test case for Measure"""

    def setUp(self):
        pass
        
    def tearDown(self):
        pass

    def testNaiveMeasure(self):
        mi = afwImage.MaskedImageF(100, 200)
        mi.set(10)
        #
        # Create our measuring engine
        #
        algorithms = ["NAIVE",]
        mp = measAlg.makeMeasurePhotometry(afwImage.makeExposure(mi))
        for a in algorithms:
            mp.addAlgorithm(a)

        pol = pexPolicy.Policy(pexPolicy.PolicyString(
            """#<?cfg paf policy?>
            NAIVE.radius: 10.0
            """
            ))

        mp.configure(pol)

        p = mp.measure(afwDetection.Peak(30, 50))

        if False:
            n = p.find(algorithms[0])

            print n.getAlgorithm(), n.getFlux()
            sch = p.getSchema()
            print [(x.getName(), x.getType(), n.get(x.getName())) for x in n.getSchema()]
            print [(x.getAlgorithm(), x.getFlux()) for x in p]
            print [(c.getAlgorithm(), [(x.getName(), x.getType(), n.get(x.getName()))
                                       for x in c.getSchema()]) for c in p]

        aName = algorithms[0]
        flux = 3170.0

        def findInvalid():
            return p.find("InvaliD")

        tests.assertRaisesLsstCpp(self, pexExceptions.NotFoundException, findInvalid)

        n = p.find(aName)
        self.assertEqual(n.getAlgorithm(), aName)
        self.assertEqual(n.getFlux(), flux)

        sch = n.getSchema()
        schEl = sch[0]
        self.assertEqual(schEl.getName(), "flux")
        self.assertEqual(n.get(schEl.getName()), flux)

        schEl = sch[1]
        self.assertEqual(schEl.getName(), "fluxErr")

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(MeasureSourcesTestCase)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
