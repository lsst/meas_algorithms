#!/usr/bin/env python
"""
Tests for Footprints, DetectionSets, and Measure

Run with:
   python Measure_1.py
or
   python
   >>> import Measure_1; Measure_1.run()
"""

import pdb                              # we may want to say pdb.set_trace()
import os, unittest
from math import *
import eups
import lsst.utils.tests as tests
import lsst.pex.logging as logging
import lsst.pex.policy as policy
import lsst.afw.detection as afwDetection
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.meas.algorithms as algorithms
import lsst.meas.algorithms.defects as defects

try:
    type(verbose)
except NameError:
    verbose = 0
logging.Trace_setVerbosity("afwDetection.Measure", verbose)

try:
    type(display)
except NameError:
    display = False

if display:
    pass
import lsst.afw.display.ds9 as ds9

def toString(*args):
    """toString written in python"""
    if len(args) == 1:
        args = args[0]

    y, x0, x1 = args
    return "%d: %d..%d" % (y, x0, x1)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class MeasureTestCase(unittest.TestCase):
    """A test case for Measure"""
    class Object(object):
        def __init__(self, val, spans):
            self.val = val
            self.spans = spans

        def insert(self, im):
            """Insert self into an image"""
            for sp in self.spans:
                y, x0, x1 = sp
                for x in range(x0, x1+1):
                    im.set(x, y, self.val)

        def __eq__(self, other):
            for osp, sp in zip(other.getSpans(), self.spans):
                if osp.toString() != toString(sp):
                    return False
                
            return True
    
    def setUp(self):
        ms = afwImage.MaskedImageF(14, 10)
        self.mi = afwImage.MaskedImageF(ms, afwImage.BBox(afwImage.PointI(1, 1), 12, 8))
        im = self.mi.getImage()
        #
        # Objects that we should detect.  These are coordinates in the subimage
        #
        self.objects = []
        self.objects += [self.Object(10, [(1, 4, 4), (2, 3, 5), (3, 4, 4)])]
        self.objects += [self.Object(20, [(5, 7, 8), (5, 10, 10), (6, 8, 9)])]
        self.objects += [self.Object(20, [(6, 3, 3)])]

        im.set(0)                       # clear image
        for obj in self.objects:
            obj.insert(im)
        
    def tearDown(self):
        del self.mi

    def testFootprintsMeasure(self):
        """Check that we can measure the objects in a detectionSet"""

        xcentroid = [5.0, 9.4, 4.0]
        ycentroid = [3.0, 6.4, 7.0]
        flux = [50.0, 100.0, 20.0]
        
        ds = afwDetection.DetectionSetF(self.mi, afwDetection.Threshold(10), "DETECTED")

        if display:
            ds9.mtv(self.mi, frame=0)

        objects = ds.getFootprints()
        source = afwDetection.Source()

        for i in range(len(objects)):
            source.setId(i)
            
            algorithms.measureSource(source, self.mi, objects[i], 0.0)

            if display:
                ds9.dot("+", source.getColc() - self.mi.getX0(), source.getRowc() - self.mi.getY0())

            self.assertAlmostEqual(source.getColc(), xcentroid[i], 6)
            self.assertAlmostEqual(source.getRowc(), ycentroid[i], 6)
            self.assertEqual(source.getFlux(), flux[i])

class FindAndMeasureTestCase(unittest.TestCase):
    """A test case detecting and measuring objects"""
    def setUp(self):
        self.mi = afwImage.MaskedImageF(os.path.join(eups.productDir("afwdata"), "CFHT", "D4", "cal-53535-i-797722_1"))

        self.FWHM = 5
        self.psf = algorithms.createPSF("DGPSF", 0, self.FWHM/(2*sqrt(2*log(2))))

        if False:                       # use full image, trimmed to data section
            self.XY0 = afwImage.PointI(32, 2)
            self.mi = self.mi.Factory(self.mi, afwImage.BBox(self.XY0, afwImage.PointI(2079, 4609)))
            self.mi.setXY0(afwImage.PointI(0, 0))
        else:                           # use sub-image
            self.XY0 = afwImage.PointI(824, 140)
            self.mi = self.mi.Factory(self.mi, afwImage.BBox(self.XY0, 256, 256))

        self.mi.getMask().addMaskPlane("DETECTED")

    def tearDown(self):
        del self.mi
        del self.psf

    def testDetection(self):
        """Test object detection"""
        #
        # Fix defects
        #
        #
        # Mask known bad pixels
        #
        badPixels = defects.policyToBadRegionList(os.path.join(eups.productDir("meas_algorithms"),
                                                               "pipeline/BadPixels.paf"))
        # did someone lie about the origin of the maskedImage?  If so, adjust bad pixel list
        if self.XY0.getX() != self.mi.getX0() or self.XY0.getY() != self.mi.getY0():
            dx = self.XY0.getX() - self.mi.getX0()
            dy = self.XY0.getY() - self.mi.getY0()
            for bp in badPixels:
                bp.shift(-dx, -dy)

        algorithms.interpolateOverDefects(self.mi, self.psf, badPixels)
        #
        # Subtract background
        #
        bctrl = afwMath.BackgroundControl(afwMath.NATURAL_SPLINE);
        bctrl.setNxSample(int(self.mi.getWidth()/256) + 1);
        bctrl.setNySample(int(self.mi.getHeight()/256) + 1);
	backobj = afwMath.make_Background(self.mi.getImage(), bctrl)

        img = self.mi.getImage(); img -= backobj.getImageF(); del img
        #
        # Remove CRs
        #
        crPolicy = policy.Policy.createPolicy(os.path.join(eups.productDir("meas_algorithms"),
                                                           "pipeline", "CosmicRays.paf"))
        crs = algorithms.findCosmicRays(self.mi, self.psf, 0, crPolicy)
        #
        # We do a pretty good job of interpolating, so don't propagagate the convolved CR/INTRP bits
        # (we'll keep them for the original CR/INTRP pixels)
        #
        savedMask = self.mi.getMask().Factory(self.mi.getMask(), True)
        saveBits = savedMask.getPlaneBitMask("CR") | \
                   savedMask.getPlaneBitMask("BAD") | \
                   savedMask.getPlaneBitMask("INTRP") # Bits to not convolve
        savedMask &= saveBits

        msk = self.mi.getMask(); msk &= ~saveBits; del msk # Clear the saved bits
        #
        # Smooth image
        #
        FWHM = 5
        psf = algorithms.createPSF("DGPSF", 15, self.FWHM/(2*sqrt(2*log(2))))

        cnvImage = self.mi.Factory(self.mi.getDimensions())
        cnvImage.setXY0(afwImage.PointI(self.mi.getX0(), self.mi.getY0()))
        psf.convolve(cnvImage, self.mi, True, savedMask.getMaskPlane("EDGE"))

        msk = cnvImage.getMask(); msk |= savedMask; del msk # restore the saved bits

        threshold = afwDetection.Threshold(3, afwDetection.Threshold.STDEV)
        #
        # Only search the part of the frame that was PSF-smoothed
        #        
        llc = afwImage.PointI(psf.getKernel().getWidth()/2, psf.getKernel().getHeight()/2)
        urc = afwImage.PointI(cnvImage.getWidth() - 1, cnvImage.getHeight() - 1) - llc;
        middle = cnvImage.Factory(cnvImage, afwImage.BBox(llc, urc))
        ds = afwDetection.DetectionSetF(middle, threshold, "DETECTED")
        del middle
        #
        # Reinstate the saved (e.g. BAD) (and also the DETECTED | EDGE) bits in the unsmoothed image
        #
        savedMask <<= cnvImage.getMask()
        msk = self.mi.getMask(); msk |= savedMask; del msk
        del savedMask

        if display:
            ds9.mtv(self.mi, frame=0)
            ds9.mtv(cnvImage, frame=1)

        objects = ds.getFootprints()
        #
        # Time to actually measure
        #
        moPolicy = policy.Policy.createPolicy(os.path.join(eups.productDir("meas_algorithms"),
                                                           "pipeline", "MeasureObjects.paf"))
        
        for i in range(len(objects)):
            source = afwDetection.Source()
            source.setId(i)
            source.setFlagForDetection(source.getFlagForDetection() | algorithms.Flags.BINNED1);

            algorithms.measureSource(source, self.mi, objects[i], moPolicy, 0.0, psf)

            if source.getFlagForDetection() & algorithms.Flags.EDGE:
                continue

            if display:
                ds9.dot("+", source.getXAstrom() - self.mi.getX0(), source.getYAstrom() - self.mi.getY0())

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    #suites += unittest.makeSuite(MeasureTestCase)
    suites += unittest.makeSuite(FindAndMeasureTestCase)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
