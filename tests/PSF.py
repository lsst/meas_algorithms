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
import lsst.afw.math as afwMath
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

class dgPsfTestCase(unittest.TestCase):
    """A test case for PSFs"""
    def setUp(self):
        self.FWHM = 5
        self.ksize = 25                      # size of desired kernel
        self.psf = algorithms.dgPSF(self.FWHM/(2*sqrt(2*log(2))), 1, 0.1, self.ksize)

    def tearDown(self):
        del self.psf

    def testKernel(self):
        """Test the creation of the PSF's kernel"""

        kim = afwImage.ImageD(self.psf.getKernel().getDimensions())
        self.psf.getKernel().computeImage(kim, False)

        self.assertTrue(kim.getWidth() == self.ksize)
        #
        # Check that the image is as expected
        #
        I0 = kim.get(self.ksize/2, self.ksize/2)
        self.assertAlmostEqual(kim.get(self.ksize/2 + 1, self.ksize/2 + 1), I0*self.psf.getValue(1, 1))
        #
        # Is image normalised?
        #
        stats = afwMath.make_Statistics(kim, afwMath.MEAN)
        self.assertAlmostEqual(self.ksize*self.ksize*stats.getValue(afwMath.MEAN), 1.0)

        if False:
            ds9.mtv(kim)        

    def testKernelConvolution(self):
        """Test convolving with the PSF"""

        im = afwImage.ImageF(100, 100)
        im.set(0)
        im.set(50, 50, 1000)

        cim = im.Factory(im.getDimensions())
        self.psf.convolve(cim, im)

        if False:
            ds9.mtv(cim)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(dgPsfTestCase)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
