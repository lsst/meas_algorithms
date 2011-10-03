#!/usr/bin/env python
"""
Tests for measuring things

Run with:
   python MeasureSources.py
or
   python
   >>> import MeasureSources; MeasureSources.run()
"""

import math, os, sys, unittest
import lsst.utils.tests as tests
import lsst.pex.exceptions as pexExceptions
import lsst.pex.logging as pexLogging
import lsst.pex.policy as pexPolicy
import lsst.afw.detection as afwDet
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.meas.algorithms as measAlg

try:
    type(verbose)
except NameError:
    verbose = 0
pexLogging.Trace_setVerbosity("afwDet.Measure", verbose)

try:
    type(display)
except NameError:
    display = False

import lsst.afw.display.ds9 as ds9

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class MeasureMultipleTestCase(unittest.TestCase):
    """A test case for Measuring multiple images"""

    def setUp(self):
        self.imageSize = (100, 200)
        self.psfSize = 20
        self.fwhm = 3.0
        self.center = (32.1, 45.6)
        self.footprint = afwDet.Footprint(afwGeom.Point2I(self.center), 2*self.fwhm)
        mi = afwImage.MaskedImageF(afwGeom.Extent2I(self.imageSize))
        mi.set(0)
        self.psf = afwDet.createPsf("DoubleGaussian", self.psfSize, self.psfSize, 
                                          self.fwhm/(2*math.sqrt(2*math.log(2))))
        self.exp = afwImage.makeExposure(mi)
        self.peak = afwDet.Peak(self.center[0], self.center[1])
        
        afwImage.Filter.define(afwImage.FilterProperty("F"))
        afwImage.Filter.define(afwImage.FilterProperty("G"))
        self.filter = afwImage.Filter("F")
        self.exp.setFilter(self.filter)
        self.exp.setPsf(self.psf)

        #self.wcs = afwImage.makeWcs((0,0), self.center, 0.2, 0.0, 0.2, 0.0)
        #self.exp.setWcs(self.wcs)

        psfImage = self.psf.computeImage(afwGeom.Point2D(self.center)).convertF()
        subImage = mi.getImage().Factory(mi.getImage(), 
                                         afwGeom.BoxI(afwGeom.Point2I(int(self.center[0] - self.psfSize/2),
                                                                      int(self.center[1] - self.psfSize/2)),
                                                      afwGeom.Extent2I(self.psfSize)), 
                                         afwImage.LOCAL)
        subImage += psfImage

        self.algName = "PSF"
        self.mp = measAlg.makeMeasurePhotometry(self.exp)
        self.mp.addAlgorithm(self.algName)
###        pol = pexPolicy.Policy(pexPolicy.PolicyString(
###            """#<?cfg paf policy?>
###            PSF.radius: 10.0
###            """
###            ))
###
###        mp.configure(pol)

    def tearDown(self):
        del self.psf
        del self.exp
        del self.filter
        del self.peak
        del self.footprint
        del self.mp


    def testInvalidFind(self):
        s = afwDet.Source(0)
        s.setFootprint(self.footprint)
        p = self.mp.measure(self.exp, self.peak, s)

        def findInvalid():
            return p.find("InvaliD")
        tests.assertRaisesLsstCpp(self, pexExceptions.NotFoundException, findInvalid)


    def testSinglePhotometry(self):
        s = afwDet.Source(0)
        s.setFootprint(self.footprint)
        p = self.mp.measure(self.exp, self.peak, s)

        n = p.find(self.algName)
        self.assertEqual(n.getAlgorithm(), self.algName)
        self.assertAlmostEqual(n.getFlux(), 1.0)

    def testInvalidGroup(self):
        patch = measAlg.makeExposurePatch(self.exp)
        exp = self.exp.Factory(self.exp, True) # Deep copy
        filt = afwImage.Filter("G")
        exp.setFilter(filt)
        group = measAlg.makeExposureGroup(patch)
        def invalidGroup():
            p = measAlg.makeExposurePatch(exp)
            group.append(p)
        tests.assertRaisesLsstCpp(self, pexExceptions.RuntimeErrorException, invalidGroup)

    def testGroupPhotometry(self):
        s = afwDet.Source(0)
        patch = measAlg.makeExposurePatch(self.exp)
        patch.setFootprint(self.footprint)
        patch.setPeak(self.peak)

        group = measAlg.makeExposureGroup(patch)
        group.append(patch)

        p = self.mp.measureGroup(group, s)

        n = p.find(self.algName)
        print n.size()
        self.assertEqual(n.getAlgorithm(), self.algName)
        self.assertAlmostEqual(n.getFlux(), 1.0)

    def testGroupsPhotometry(self):
        s = afwDet.Source(0)
        patch = measAlg.makeExposurePatch(self.exp)
        patch.setFootprint(self.footprint)
        patch.setPeak(self.peak)

        group = measAlg.makeExposureGroup(patch)
        group.append(patch)

        groups = measAlg.ExposureGroupSetF([group, group])

        p = self.mp.measureGroups(groups, s)

        pPsf = p.find(self.algName)
        self.assertEqual(pPsf.size(), groups.size())

        for n in pPsf:
            self.assertEqual(n.getAlgorithm(), self.algName)
            self.assertAlmostEqual(n.getFlux(), 1.0)
  

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(MeasureMultipleTestCase)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
