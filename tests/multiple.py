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
import lsst.afw.coord as afwCoord
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

PSF_PLACES = 5                          # "almost" accuracy for PSF measurement

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class MeasureMultipleTestCase(unittest.TestCase):
    """A test case for Measuring multiple images"""

    def setUp(self):
        self.imageSize = (100, 200)
        self.psfSize = 20
        self.fwhm = 3.0
        self.center = afwGeom.Point2D(32.1, 45.6)
        self.footprint = afwDet.Footprint(afwGeom.Point2I(self.center), 2*self.fwhm)
        mi = afwImage.MaskedImageF(afwGeom.Extent2I(self.imageSize))
        mi.set(0, 0, 0.1)
        self.psf = afwDet.createPsf("DoubleGaussian", self.psfSize, self.psfSize, 
                                          self.fwhm/(2*math.sqrt(2*math.log(2))))
        self.exp = afwImage.makeExposure(mi)
        
        afwImage.Filter.define(afwImage.FilterProperty("F"))
        afwImage.Filter.define(afwImage.FilterProperty("G"))
        self.filter = afwImage.Filter("F")
        self.exp.setFilter(self.filter)
        self.exp.setPsf(self.psf)

        self.wcs = afwImage.makeWcs(afwCoord.Coord(afwGeom.Point2D(180, 0), afwGeom.degrees), 
                                    self.center, 0.2 / 3600, 0.0, 0.0, 0.2 / 3600)
        self.exp.setWcs(self.wcs)

        psfImage = self.psf.computeImage(self.center, False).convertF()
        psfBox = psfImage.getBBox(afwImage.LOCAL)
        psfBox.shift(psfImage.getXY0() - mi.getImage().getXY0());
        subImage = mi.getImage().Factory(mi.getImage(), psfBox, afwImage.LOCAL)
        subImage += psfImage

        if display:
            ds9.mtv(self.exp)
            ds9.dot("+", self.center[0], self.center[1])

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
        del self.wcs
        del self.filter
        del self.footprint
        del self.mp


    def testInvalidFind(self):
        s = afwDet.Source(0)
        s.setFootprint(self.footprint)
        p = self.mp.measure(s, self.exp, self.center)

        def findInvalid():
            return p.find("InvaliD")
        tests.assertRaisesLsstCpp(self, pexExceptions.NotFoundException, findInvalid)


    def testSinglePhotometry(self):
        s = afwDet.Source(0)
        s.setFootprint(self.footprint)
        s.setXAstrom(self.center.getX())
        s.setYAstrom(self.center.getY())
        p = self.mp.measure(s, self.exp, self.center)

        n = p.find(self.algName)
        self.assertEqual(n.getAlgorithm(), self.algName)
        self.assertAlmostEqual(n.getFlux(), 1.0, places=PSF_PLACES)

    def testMultiplePhotometry(self):
        s = afwDet.Source(0)
        s.setFootprint(self.footprint)
        s.setXAstrom(self.center.getX())
        s.setYAstrom(self.center.getY())

        exposures = measAlg.ExposureListF()
        exposures.push_back(self.exp)

        p = self.mp.measure(s, s, self.wcs, exposures)

        n = p.find(self.algName)
        print n.size()
        self.assertEqual(n.getAlgorithm(), self.algName)
        self.assertAlmostEqual(n.getFlux(), 1.0, places=PSF_PLACES)
  

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
