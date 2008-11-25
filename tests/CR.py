#!/usr/bin/env python
"""
Tests for cosmic ray detection

Run with:
   python CR.py
or
   python
   >>> import CR; CR.run()
"""

import pdb                              # we may want to say pdb.set_trace()
import os
import sys
from math import *
import unittest
import eups
import lsst.utils.tests as tests
import lsst.pex.logging as logging
import lsst.pex.policy as policy
import lsst.afw.image.imageLib as afwImage
import lsst.afw.math.mathLib as afwMath
import lsst.afw.detection.detectionLib as afwDetection
import lsst.afw.display.ds9 as ds9
import lsst.detection.detectionLib as detection
import lsst.detection.defects as defects

try:
    type(verbose)
except NameError:
    verbose = 0
    logging.Trace_setVerbosity("detection.CR", verbose)

try:
    type(display)
except NameError:
    display = False

    if display:
        import lsst.afw.display.ds9 as ds9

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class CosmicRayTestCase(unittest.TestCase):
    """A test case for Cosmic Ray detection"""
    def setUp(self):
        self.FWHM = 5                   # pixels
        self.psf = detection.dgPSF(self.FWHM/(2*sqrt(2*log(2))))
            
        self.mi = afwImage.MaskedImageF(os.path.join(eups.productDir("afwdata"), "CFHT", "D4", "cal-53535-i-797722_1"))

        if False:                           # use full image
            self.nCR = 1094                 # number of CRs we should detect
        else:                               # use sub-image
            self.mi = self.mi.Factory(self.mi, afwImage.BBox(afwImage.PointI(824, 140), 256, 256))
            self.nCR = 13

        self.mi.getMask().addMaskPlane("DETECTED")

        self.policy = policy.Policy.createPolicy(os.path.join(eups.productDir("detection"),
                                                              "pipeline", "CosmicRays.paf"))

    def tearDown(self):
        del self.mi
        del self.psf
        del self.policy

    def testDetection(self):
        """Test CR detection"""

        if display:
            frame = 0
            ds9.mtv(self.mi, frame=frame) # raw frame
            if self.mi.getX0() == 0:
                ds9.pan(944, 260)
        #
        # Mask known bad pixels
        #
        badPixels = defects.policyToBadRegionList(os.path.join(os.environ["DETECTION_DIR"], "pipeline/BadPixels.paf"))
        detection.interpolateOverDefects(self.mi, self.psf, badPixels)

        background = afwMath.StatisticsF(self.mi.getImage(), afwMath.MEAN).getValue(afwMath.MEAN)
        crs = detection.findCosmicRays(self.mi, self.psf, background, self.policy)

        if display:
            ds9.mtv(self.mi.getImage(), frame=frame+1)
            if self.mi.getX0() == 0:
                ds9.pan(944, 260)

            ds9.mtv(self.mi, frame=frame+2)
            if self.mi.getX0() == 0:
                ds9.pan(944, 260)

        self.assertEqual(len(crs), self.nCR)
        print "Detected %d CRs" % len(crs)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    if eups.productDir("afwdata"):
        suites += unittest.makeSuite(CosmicRayTestCase)
    else:
        print >> sys.stderr, "afwdata is not setup; skipping CosmicRayTestCase"
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)


def run(exit=False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
