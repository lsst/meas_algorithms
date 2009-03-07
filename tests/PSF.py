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
        self.psf = algorithms.createPSF("DoubleGaussian", self.ksize, self.ksize, self.FWHM/(2*sqrt(2*log(2))), 1, 0.1)

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
            psf = algorithms.createPSF("DoubleGaussian", 0, 0, 1)
            psf.convolve(cim, im)

        utilsTests.assertRaisesLsstCpp(self, pexExceptions.RuntimeErrorException, badKernelSize)

    def testInvalidDgPSF(self):
        """Test parameters of dgPSFs, both valid and not"""
        sigma1, sigma2, b = 1, 0, 0                     # sigma2 may be 0 iff b == 0
        algorithms.createPSF("DoubleGaussian", self.ksize, self.ksize, sigma1, sigma2, b)

        def badSigma1():
            sigma1 = 0
            algorithms.createPSF("DoubleGaussian", self.ksize, self.ksize, sigma1, sigma2, b)

        utilsTests.assertRaisesLsstCpp(self, pexExceptions.DomainErrorException, badSigma1)

        def badSigma2():
            sigma2, b = 0, 1
            algorithms.createPSF("DoubleGaussian", self.ksize, self.ksize, sigma1, sigma2, b)

        utilsTests.assertRaisesLsstCpp(self, pexExceptions.DomainErrorException, badSigma2)


    def testGetImage(self):
        """Test returning a realisation of the PSF; test the sanity of the SDSS centroider at the same time"""

        xcen = self.psf.getWidth()//2
        ycen = self.psf.getHeight()//2

        centroider = algorithms.createMeasureCentroid("SDSS")

        stamps = []
        trueCenters = []
        centroids = []
        for x, y in ([10, 10], [9.4999, 10.4999], [10.5001, 10.5001]):
            fx, fy = x - int(x), y - int(y)
            if fx >= 0.5:
                fx -= 1.0
            if fy >= 0.5:
                fy -= 1.0

            im = self.psf.getImage(x, y).convertFloat()

            c = centroider.apply(im, xcen, ycen, None, 0.0)

            stamps.append(im.Factory(im, True))
            centroids.append([c.getX(), c.getY()])
            trueCenters.append([xcen + fx, ycen + fy])
            
        if display:
            mos = displayUtils.Mosaic()     # control mosaics
            ds9.mtv(mos.makeMosaic(stamps))

            for i in range(len(trueCenters)):
                bbox = mos.getBBox(i)

                ds9.dot("+",
                        bbox.getX0() + xcen, bbox.getY0() + ycen, ctype=ds9.RED, size=1)
                ds9.dot("+",
                        bbox.getX0() + centroids[i][0], bbox.getY0() + centroids[i][1], ctype=ds9.YELLOW, size=1.5)
                ds9.dot("+",
                        bbox.getX0() + trueCenters[i][0], bbox.getY0() + trueCenters[i][1])

                ds9.dot("%.2f, %.2f" % (trueCenters[i][0], trueCenters[i][1]),
                        bbox.getX0() + xcen, bbox.getY0() + 2)

        for i in range(len(centroids)):
            self.assertAlmostEqual(centroids[i][0], trueCenters[i][0], 4)
            self.assertAlmostEqual(centroids[i][1], trueCenters[i][1], 4)

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
        self.psf = algorithms.createPSF("DoubleGaussian", self.ksize, self.ksize, self.FWHM/(2*sqrt(2*log(2))), 1, 0.1)

        for x, y in [(20, 20), (30, 30), (50, 50), (60, 20), (60, 210)]:
            source = afwDetection.Source()

            flux = 10000 - 100*x - 10*y

            self.mi.getImage().set(x, y, flux)

        psf = algorithms.createPSF("DoubleGaussian", self.ksize, self.ksize, self.FWHM/(2*sqrt(2*log(2))), 1, 0.1)
        cim = self.mi.Factory(self.mi.getDimensions())
        psf.convolve(cim, self.mi)
        self.mi = cim

        self.cellSet = afwMath.SpatialCellSet(afwImage.BBox(afwImage.PointI(0, 0), width, height), 100)
        ds = afwDetection.DetectionSetF(self.mi, afwDetection.Threshold(110), "DETECTED")
        objects = ds.getFootprints()
        #
        # Prepare to measure
        #
        self.moPolicy = policy.Policy.createPolicy(os.path.join(eups.productDir("meas_algorithms"),
                                                                "pipeline", "MeasureSources.paf"))
        
        measureSources = algorithms.makeMeasureSources(afwImage.makeExposure(self.mi), self.moPolicy, psf)

        sourceList = afwDetection.SourceSet()
        for i in range(len(objects)):
            source = afwDetection.Source()
            sourceList.append(source)

            source.setId(i)
            source.setFlagForDetection(source.getFlagForDetection() | algorithms.Flags.BINNED1);

            measureSources.apply(source, objects[i])

            self.cellSet.insertCandidate(algorithms.makePsfCandidate(source, self.mi))

    def testGetPcaKernel(self):
        """Convert our cellSet to a LinearCombinationKernel"""

        nEigenComponents = self.moPolicy.getInt("determinePsf.nEigenComponents")
        spatialOrder  = self.moPolicy.getInt("determinePsf.spatialOrder")
        nStarPerCell = self.moPolicy.getInt("determinePsf.nStarPerCell")
        kernelSize = self.moPolicy.getInt("determinePsf.kernelSize")

        pair = algorithms.createKernelFromPsfCandidates(self.cellSet,
                                                       nEigenComponents, spatialOrder, kernelSize, nStarPerCell)
        kernel, eigenValues = pair[0], pair[1]
                                                              
        psf = algorithms.createPSF("PCA", kernel)
        if False:
            psf.getImage(0,0)               # calculate kernelParameters
            print "PSF parameters:", psf.getKernel().getKernelParameters()

        self.assertTrue(afwMath.cast_AnalyticKernel(psf.getKernel()) is None)
        self.assertTrue(afwMath.cast_LinearCombinationKernel(psf.getKernel()) is not None)
            
        if display:
            eImages = []
            for k in afwMath.cast_LinearCombinationKernel(psf.getKernel()).getKernelList():
                pass

            for k in psf.getKernel().getKernelList():
                im = afwImage.ImageD(k.getDimensions())
                k.computeImage(im, k, True)
                eImages.append(im)

            mos = displayUtils.Mosaic()
            frame = 3
            ds9.mtv(mos.makeMosaic(eImages), frame=frame)
            ds9.dot("Eigen Images", 0, 0, frame=frame)

            #
            # Make a mosaic of PSF candidates
            #
            stamps = []
            stampInfo = []

            for cell in self.cellSet.getCellList():
                for cand in cell:
                    #
                    # Swig doesn't know that we inherited from SpatialCellImageCandidate;  all
                    # it knows is that we have a SpatialCellCandidate, and SpatialCellCandidates
                    # don't know about getImage;  so cast the pointer to PsfCandidate
                    #
                    cand = algorithms.cast_PsfCandidateF(cand)
                    s = cand.getSource()

                    im = cand.getImage()

                    stamps.append(im)
                    stampInfo.append("[%d 0x%x]" % (s.getId(), s.getFlagForDetection()))
        
            ds9.mtv(self.mi, frame=0)
            
            mos = displayUtils.Mosaic()
            frame = 1
            ds9.mtv(mos.makeMosaic(stamps), frame=frame, lowOrderBits=True)
            for i in range(len(stampInfo)):
                ds9.dot(stampInfo[i], mos.getBBox(i).getX0(), mos.getBBox(i).getY0(), frame=frame, ctype=ds9.RED)

            frame = 2
            ds9.mtv(psf.getImage(0, 0), frame=frame)
            ds9.dot("PSF", 0, 0, frame=frame)

    def tearDown(self):
        del self.cellSet
        del self.mi

    def testCandidateList(self):
        if False and display:
            ds9.mtv(self.mi)

            for cell in self.cellSet.getCellList():
                x0, y0, x1, y1 = \
                    cell.getBBox().getX0(), cell.getBBox().getY0(), cell.getBBox().getX1(), cell.getBBox().getY1()
                print x0, y0, " ", x1, y1
                x0 -= 0.5; y0 -= 0.5
                x1 += 0.5; y1 += 0.5

                ds9.line([(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)], ctype=ds9.RED)

        self.assertFalse(self.cellSet.getCellList()[0].empty())
        self.assertTrue(self.cellSet.getCellList()[1].empty())
        self.assertFalse(self.cellSet.getCellList()[2].empty())

        stamps = []
        stampInfo = []
        for cell in self.cellSet.getCellList():
            for cand in cell:
                #
                # Swig doesn't know that we inherited from SpatialCellImageCandidate;  all
                # it knows is that we have a SpatialCellCandidate, and SpatialCellCandidates
                # don't know about getImage;  so cast the pointer to SpatialCellImageCandidate<Image<float> >
                # and all will be well
                #
                cand = afwMath.cast_SpatialCellImageCandidateMF(cell[0])
                width, height = 15, 17
                cand.setWidth(width); cand.setHeight(height);

                im = cand.getImage()
                stamps.append(im)

                self.assertEqual(im.getWidth(), width)
                self.assertEqual(im.getHeight(), height)
        
        if display:
            mos = displayUtils.Mosaic()
            ds9.mtv(mos.makeMosaic(stamps), frame=1)

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
else:
    def dummy_assertRaisesLsstCpp(this, exception, test):
        """Disable assertRaisesLsstCpp as at it fails when run from the python prompt; #656"""
        print "Not running test %s" % test

    utilsTests.assertRaisesLsstCpp = dummy_assertRaisesLsstCpp

    
