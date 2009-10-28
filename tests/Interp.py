#!/usr/bin/env python
"""
Tests for bad pixel interpolation code

Run with:
   python Interp.py
or
   python
   >>> import Interp; Interp.run()
"""

import pdb                              # we may want to say pdb.set_trace()
import os
from math import *
import unittest
import eups
import lsst.utils.tests as tests
import lsst.pex.logging as logging
import lsst.afw.image as afwImage
import lsst.afw.display.ds9 as ds9
import lsst.meas.algorithms as algorithms
import lsst.meas.algorithms.defects as defects

try:
    type(verbose)
except NameError:
    verbose = 0
    logging.Trace_setVerbosity("algorithms.Interp", verbose)

try:
    type(display)
except NameError:
    display = False

    if display:
        import lsst.afw.display.ds9 as ds9

class interpolationTestCase(unittest.TestCase):
    """A test case for interpolation"""
    def setUp(self):
        self.FWHM = 5
        self.psf = algorithms.createPSF("DoubleGaussian", 0, 0, self.FWHM/(2*sqrt(2*log(2))))
        maskedImageFile = os.path.join(eups.productDir("afwdata"), "CFHT", "D4", "cal-53535-i-797722_1")
            
        self.mi = afwImage.MaskedImageF(maskedImageFile)
        if False:                       # use sub-image?
            self.mi = self.mi.Factory(self.mi, afwImage.BBox(afwImage.PointI(760, 20), 256, 256))
        self.mi.getMask().addMaskPlane("INTERP")

        self.badPixels = defects.policyToBadRegionList(os.path.join(eups.productDir("meas_algorithms"),
                                                                    "pipeline/BadPixels.paf"))

    def tearDown(self):
        del self.mi
        del self.psf
        del self.badPixels

    def testDetection(self):
        """Test Interp algorithms"""

        if display:
            frame = 0
            ds9.mtv(self.mi, frame = frame)

        algorithms.interpolateOverDefects(self.mi, self.psf, self.badPixels)

        if display:
            ds9.mtv(self.mi, frame = frame + 1)

    def test818(self):
        """A test case for #818; the full test is in /lsst/DC3root/ticketFiles/818"""

        badPixels = algorithms.DefectListT()
        defects = [((82, 663), 6, 8),
                   ((83, 659), 9, 6),
                   ((85, 660), 10, 11),
                   ((87, 669), 3, 3),
                   ]

        for xy0, width, height in defects:
            x0, y0 = xy0
            bbox = afwImage.BBox(afwImage.PointI(x0, y0), width, height)
            badPixels.push_back(algorithms.Defect(bbox))

        mi = afwImage.MaskedImageF(517, 800)

        algorithms.interpolateOverDefects(mi, self.psf, badPixels)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    if eups.productDir("afwdata"):
        suites += unittest.makeSuite(interpolationTestCase)
    else:
        print "Skipping interpolation test case as afwdata isn't set up"
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)


def run(exit = False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
