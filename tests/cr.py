#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#

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
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.afw.detection as afwDetection
import lsst.afw.display.ds9 as ds9
import lsst.meas.algorithms as algorithms
import lsst.meas.algorithms.defects as defects

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

import lsst.afw.display.ds9 as ds9

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class CosmicRayTestCase(unittest.TestCase):
    """A test case for Cosmic Ray detection"""
    def setUp(self):
        self.FWHM = 5                   # pixels
        self.psf = afwDetection.createPsf("DoubleGaussian", 29, 29, self.FWHM/(2*sqrt(2*log(2))))
            
        self.mi = afwImage.MaskedImageF(imageFile)
        self.XY0 = afwGeom.PointI(0, 0) # origin of the subimage we use

        if imageFile == imageFile0:
            if True:                        # use full image, trimmed to data section
                self.XY0 = afwGeom.PointI(32, 2)
                self.mi = self.mi.Factory(self.mi, afwGeom.BoxI(self.XY0, afwGeom.PointI(2079, 4609)), afwImage.LOCAL)
                self.mi.setXY0(afwGeom.PointI(0, 0))
                self.nCR = 1076                 # number of CRs we should detect
            else:                               # use sub-image
                if True:
                    self.XY0 = afwGeom.PointI(824, 140)
                    self.nCR = 10
                else:
                    self.XY0 = afwGeom.PointI(280, 2750)
                    self.nCR = 2
                self.mi = self.mi.Factory(self.mi, afwGeom.BoxI(self.XY0, afwGeom.ExtentI(256, 256), afwImage.LOCAL))
                self.mi.setXY0(afwGeom.PointI(0, 0))
        else:
            self.nCR = None

        self.mi.getMask().addMaskPlane("DETECTED")

    def tearDown(self):
        del self.mi
        del self.psf

    def testDetection(self):
        """Test CR detection"""
        #
        # Subtract background
        #
        bctrl = afwMath.BackgroundControl(afwMath.Interpolate.NATURAL_SPLINE);
        bctrl.setNxSample(int(self.mi.getWidth()/256) + 1);
        bctrl.setNySample(int(self.mi.getHeight()/256) + 1);
        bctrl.getStatisticsControl().setNumSigmaClip(3.0)  
        bctrl.getStatisticsControl().setNumIter(2)

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

        crConfig = algorithms.FindCosmicRaysConfig()
        crs = algorithms.findCosmicRays(self.mi, self.psf, background, crConfig)

        if display:
            ds9.mtv(self.mi, frame = frame + 1, title="CRs removed")
            if self.mi.getWidth() > 256:
                ds9.pan(944 - self.mi.getX0(), 260 - self.mi.getY0())

        print "Detected %d CRs" % len(crs)
        if display and False:
            for cr in crs:
                bbox = cr.getBBox()
                bbox.shift(afwGeom.ExtentI(-self.mi.getX0(), -self.mi.getY0()))
                ds9.line([(bbox.getMinX() - 0.5, bbox.getMinY() - 0.5),
                          (bbox.getMaxX() + 0.5, bbox.getMinY() - 0.5),
                          (bbox.getMaxX() + 0.5, bbox.getMaxY() + 0.5),
                          (bbox.getMinX() - 0.5, bbox.getMaxY() + 0.5),
                          (bbox.getMinX() - 0.5, bbox.getMinY() - 0.5)], frame = frame + 1)

        if self.nCR is not None:
            self.assertEqual(len(crs), self.nCR)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    if imageFile0:
        suites += unittest.makeSuite(CosmicRayTestCase)
    else:
        print >> sys.stderr, "afwdata is not setup; skipping CosmicRayTestCase"
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)


def run(exit = False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
