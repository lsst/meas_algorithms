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
import lsst.afw.detection.detectionLib as detection
import lsst.afw.display.ds9 as ds9
import lsst.meas.algorithms as algorithms
import lsst.meas.algorithms.defects as defects

display = True

try:
    type(verbose)
except NameError:
    verbose = 0
logging.Trace_setVerbosity("algorithms.CR", verbose)

try:
    type(display)
except NameError:
    display = False

    if eups.productDir("afwdata"):
        imageFile0 = os.path.join(eups.productDir("afwdata"), "CFHT", "D4", "cal-53535-i-797722_1")
    else:
        imageFile0 = None
    imageFile = imageFile0

if display:
    import lsst.afw.display.ds9 as ds9

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class CosmicRayTestCase(unittest.TestCase):
    """A test case for Cosmic Ray detection"""
    def setUp(self):
        self.FWHM = 5                   # pixels
        self.psf = algorithms.createPSF("DoubleGaussian", 0, 0, self.FWHM/(2*sqrt(2*log(2))))
            
#        self.mi = afwImage.MaskedImageF(os.path.join(eups.productDir("afwdata"), "CFHT", "D4", "cal-53535-i-797722_1"))
#        self.img = afwImage.ImageF("/data/CTIO_Sep09/n1/obj1030_6.fits")
#        self.mask = afwImage.MaskU(self.img.getDimensions())
#        self.var = afwImage.ImageF(self.img.getDimensions())
#        self.mi = afwImage.makeMaskedImage(self.img, self.mask, self.var)
        self.mi = afwImage.MaskedImageF("/data/mophot/work/test74/IP/output/isr/calobj1030/c003-a00")

        self.XY0 = afwImage.PointI(0, 0) # origin of the subimage we use

        if True:                        # use full image, trimmed to data section
            self.XY0 = afwImage.PointI(32, 2)
#            self.mi = self.mi.Factory(self.mi, afwImage.BBox(self.XY0, afwImage.PointI(2079, 4609)))
            self.mi.setXY0(afwImage.PointI(0, 0))
            self.nCR = 1088                 # number of CRs we should detect
        else:                               # use sub-image
            if True:
                self.XY0 = afwImage.PointI(824, 140)
                self.nCR = 10
            else:
                self.XY0 = afwImage.PointI(280, 2750)
                self.nCR = 2
            self.mi = self.mi.Factory(self.mi, afwImage.BBox(self.XY0, 256, 256))
            self.mi.setXY0(afwImage.PointI(0, 0))

        self.mi.getMask().addMaskPlane("DETECTED")

        self.policy = policy.Policy.createPolicy(os.path.join(eups.productDir("meas_algorithms"),
                                                              "policy", "CosmicRays.paf"))
    def tearDown(self):
        del self.mi
        del self.psf
        del self.policy

    def testDetection(self):
        """Test CR detection"""
        #
        # Subtract background
        #
        bctrl = afwMath.BackgroundControl(afwMath.Interpolate.NATURAL_SPLINE);
        bctrl.setNxSample(int(self.mi.getWidth()/256) + 1);
        bctrl.setNySample(int(self.mi.getHeight()/256) + 1);
        bctrl.sctrl.setNumSigmaClip(3.0)  
        bctrl.sctrl.setNumIter(2)

        im = self.mi.getImage()
        try:
            backobj = afwMath.makeBackground(im, bctrl)
        except Exception, e:
            print >> sys.stderr, e,

            bctrl.setInterpStyle(afwMath.Interpolate.CONSTANT)
            backobj = afwMath.makeBackground(im, bctrl)
            
        im -= backobj.getImageF()

        if display:
            frame = 0
            ds9.mtv(self.mi, frame = frame, title="Raw") # raw frame
            if self.mi.getWidth() > 256:
                ds9.pan(944 - self.mi.getX0(), 260 - self.mi.getY0())
        #
        # Mask known bad pixels
        #
        badPixels = defects.policyToBadRegionList(os.path.join(os.environ["MEAS_ALGORITHMS_DIR"],
                                                               "policy", "BadPixels.paf"))
        # did someone lie about the origin of the maskedImage?  If so, adjust bad pixel list
        if self.XY0.getX() != self.mi.getX0() or self.XY0.getY() != self.mi.getY0():
            dx = self.XY0.getX() - self.mi.getX0()
            dy = self.XY0.getY() - self.mi.getY0()
            for bp in badPixels:
                bp.shift(-dx, -dy)

        algorithms.interpolateOverDefects(self.mi, self.psf, badPixels)

        stats = afwMath.makeStatistics(self.mi.getImage(), afwMath.MEANCLIP | afwMath.STDEVCLIP)
        background = stats.getValue(afwMath.MEANCLIP)

        crPolicy = self.policy.getPolicy('CR')
        crs = algorithms.findCosmicRays(self.mi, self.psf, background, crPolicy)

        if display:
            ds9.mtv(self.mi, frame = frame + 1, title="CRs removed")
            if self.mi.getWidth() > 256:
                ds9.pan(944 - self.mi.getX0(), 260 - self.mi.getY0())

        print "Detected %d CRs" % len(crs)
        if display and False:
            for cr in crs:
                bbox = cr.getBBox()
                bbox.shift(-self.mi.getX0(), -self.mi.getY0())
                ds9.line([(bbox.getX0() - 0.5, bbox.getY0() - 0.5),
                          (bbox.getX1() + 0.5, bbox.getY0() - 0.5),
                          (bbox.getX1() + 0.5, bbox.getY1() + 0.5),
                          (bbox.getX0() - 0.5, bbox.getY1() + 0.5),
                          (bbox.getX0() - 0.5, bbox.getY0() - 0.5)], frame = frame + 1)

        if self.nCR is not None:
            self.assertEqual(len(crs), self.nCR)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
#    if eups.productDir("afwdata"):
#        suites += unittest.makeSuite(CosmicRayTestCase)
#    else:
#        print >> sys.stderr, "afwdata is not setup; skipping CosmicRayTestCase"
    suites += unittest.makeSuite(CosmicRayTestCase)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)


def run(exit = False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
