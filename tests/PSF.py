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
import lsst.utils.tests as utilsTests
import lsst.pex.exceptions as pexExceptions
import lsst.pex.logging as logging
import lsst.pex.policy as policy
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetection
import lsst.afw.math as afwMath
import lsst.afw.display.ds9 as ds9
import lsst.afw.display.utils as displayUtils
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

class dgPsfTestCase(unittest.TestCase):
    """A test case for dgPSFs"""
    def setUp(self):
        self.FWHM = 5
        self.ksize = 25                      # size of desired kernel
        self.psf = algorithms.createPSF("DGPSF", self.ksize, self.FWHM/(2*sqrt(2*log(2))), 1, 0.1)

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
        stats = afwMath.makeStatistics(kim, afwMath.MEAN)
        self.assertAlmostEqual(self.ksize*self.ksize*stats.getValue(afwMath.MEAN), 1.0)

        if False:
            ds9.mtv(kim)        

    def testKernelConvolution(self):
        """Test convolving with the PSF"""

        for im in (afwImage.ImageF(100, 100), afwImage.MaskedImageF(100, 100)):
            im.set(0)
            im.set(50, 50, 1000)

            cim = im.Factory(im.getDimensions())
            self.psf.convolve(cim, im)

            if False:
                ds9.mtv(cim)
        #
        # Check that a PSF with a zero-sized kernel can't be used to convolve
        #
        def badKernelSize():
            psf = algorithms.createPSF("DGPSF", 0, 1)
            psf.convolve(cim, im)

        utilsTests.assertRaisesLsstCpp(self, pexExceptions.RuntimeErrorException, badKernelSize)

    def testInvalidDgPSF(self):
        """Test parameters of dgPSFs, both valid and not"""
        sigma1, sigma2, b = 1, 0, 0                     # sigma2 may be 0 iff b == 0
        algorithms.createPSF("DGPSF", self.ksize, sigma1, sigma2, b)

        def badSigma1():
            sigma1 = 0
            algorithms.createPSF("DGPSF", self.ksize, sigma1, sigma2, b)

        utilsTests.assertRaisesLsstCpp(self, pexExceptions.DomainErrorException, badSigma1)

        def badSigma2():
            sigma2, b = 0, 1
            algorithms.createPSF("DGPSF", self.ksize, sigma1, sigma2, b)

        utilsTests.assertRaisesLsstCpp(self, pexExceptions.DomainErrorException, badSigma2)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class SpatialModelPsfTestCase(unittest.TestCase):
    """A test case for SpatialModelPsf"""

    def setUp(self):
        width, height = 100, 300
        self.mi = afwImage.MaskedImageF(width, height)
        self.mi.set(100)
        self.mi.getVariance().set(10)
        self.mi.getMask().addMaskPlane("DETECTED")

        self.FWHM = 5
        self.ksize = 15                      # size of desired kernel
        self.psf = algorithms.createPSF("DGPSF", self.ksize, self.FWHM/(2*sqrt(2*log(2))), 1, 0.1)

        for x, y in [(20, 20), (30, 30), (50, 50), (60, 20), (60, 210)]:
            source = afwDetection.Source()

            flux = 10000 - 100*x - 10*y

            self.mi.getImage().set(x, y, flux)

        psf = algorithms.createPSF("DGPSF", self.ksize, self.FWHM/(2*sqrt(2*log(2))), 1, 0.1)
        cim = self.mi.Factory(self.mi.getDimensions())
        psf.convolve(cim, self.mi)
        self.mi = cim

        self.cellSet = afwMath.SpatialCellSet(afwImage.BBox(afwImage.PointI(0, 0), width, height), 100)
        ds = afwDetection.DetectionSetF(self.mi, afwDetection.Threshold(110), "DETECTED")
        objects = ds.getFootprints()
        #
        # Prepare to measure
        #
        moPolicy = policy.Policy.createPolicy(os.path.join(eups.productDir("meas_algorithms"),
                                                           "pipeline", "MeasureSources.paf"))
        
        measureSources = algorithms.makeMeasureSources(afwImage.makeExposure(self.mi), moPolicy, psf)

        sourceList = afwDetection.SourceSet()
        for i in range(len(objects)):
            source = afwDetection.Source()
            sourceList.append(source)

            source.setId(i)
            source.setFlagForDetection(source.getFlagForDetection() | algorithms.Flags.BINNED1);

            measureSources.apply(source, objects[i])

            self.cellSet.insertCandidate(algorithms.makePsfCandidate(source, self.mi))

    def tearDown(self):
        del self.cellSet
        del self.mi

    def testCandidateList(self):
        if False and display:
            ds9.mtv(self.mi)

            for i in range(len(self.cellSet.getCellList())):
                cell = self.cellSet.getCellList()[i]
                x0, y0, x1, y1 = \
                    cell.getBBox().getX0(), cell.getBBox().getY0(), cell.getBBox().getX1(), cell.getBBox().getY1()
                print x0, y0, " ", x1, y1
                x0 -= 0.5; y0 -= 0.5
                x1 += 0.5; y1 += 0.5

                ds9.line([(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)], ctype=ds9.RED)

        self.assertEqual(self.cellSet.getCellList()[0].isUsable(), True)
        self.assertEqual(self.cellSet.getCellList()[1].isUsable(), False)
        self.assertEqual(self.cellSet.getCellList()[2].isUsable(), True)

        stamps = []
        stampInfo = []
        for i in range(len(self.cellSet.getCellList())):
            cell = self.cellSet.getCellList()[i]
            if not cell.isUsable():
                continue

            cand = cell.getCurrentCandidate()
            #
            # Swig doesn't know that we inherited from SpatialCellImageCandidate;  all
            # it knows is that we have a SpatialCellCandidate, and SpatialCellCandidates
            # don't know about getImage;  so cast the pointer to SpatialCellImageCandidate<Image<float> >
            # and all will be well
            #
            cand = afwMath.cast_SpatialCellImageCandidateMF(cand)
            width, height = 15, 17
            cand.setWidth(width); cand.setHeight(height);

            im = cand.getImage()
            stamps.append(im)
            stampInfo.append(i)

            self.assertEqual(im.getWidth(), width)
            self.assertEqual(im.getHeight(), height)
        
        if display:
            mos = displayUtils.Mosaic()
            ds9.mtv(mos.makeMosaic(stamps), frame=1)
            for i in range(len(stampInfo)):
                ds9.dot(stampInfo[i], mos.getBBox(i).getX0(), mos.getBBox(i).getY0(), frame=1, ctype=ds9.RED)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(dgPsfTestCase)
    suites += unittest.makeSuite(SpatialModelPsfTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the utilsTests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
