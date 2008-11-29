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
import eups
import lsst.utils.tests as tests
import lsst.pex.logging as logging
import lsst.afw.detection as afwDetection
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.meas.algorithms as algorithms

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
        self.ms = afwImage.MaskedImageF(ms, afwImage.BBox(afwImage.PointI(1, 1), 12, 8))
        im = self.ms.getImage()
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
        del self.ms

    def testFootprintsMeasure(self):
        """Check that we can measure the objects in a detectionSet"""

        xcentroid = [5.0, 9.4, 4.0]
        ycentroid = [3.0, 6.4, 7.0]
        flux = [50.0, 100.0, 20.0]
        
        ds = afwDetection.DetectionSetF(self.ms, afwDetection.Threshold(10), "DETECTED")

        if display:
            ds9.mtv(self.ms, frame=0)

        objects = ds.getFootprints()
        source = afwDetection.Source()

        for i in range(len(objects)):
            source.setId(i)
            
            algorithms.measureSource(source, self.ms, objects[i], 0.0)

            if display:
                ds9.dot("+", source.getColc() - self.ms.getX0(), source.getRowc() - self.ms.getY0())

            self.assertAlmostEqual(source.getColc(), xcentroid[i], 6)
            self.assertAlmostEqual(source.getRowc(), ycentroid[i], 6)
            self.assertEqual(source.getFlux(), flux[i])

class FindAndMeasureTestCase(unittest.TestCase):
    """A test case detecting and measuring objects"""
    def setUp(self):
        self.ms = afwImage.MaskedImageF(os.path.join(eups.productDir("afwdata"), "CFHT", "D4", "cal-53535-i-797722_1"))

        if False:                           # use full image
            pass
        else:                               # use sub-image
            self.ms = self.ms.Factory(self.ms, afwImage.BBox(afwImage.PointI(824, 140), 256, 256))

        self.ms.getMask().addMaskPlane("DETECTED")

    def tearDown(self):
        del self.ms

    def testDetection(self):
        """Test object detection"""

        stats = afwMath.StatisticsF(self.ms.getImage(), afwMath.MEAN)
        img = self.ms.getImage()
        img -= stats.getValue(afwMath.MEAN); del img

        if False and display:
            frame = 0
            ds9.mtv(self.ms, frame=frame) # raw frame
                    
        threshold = afwDetection.Threshold(5, afwDetection.Threshold.STDEV)
        if False:
            #
            # Smooth image
            #
            kFunc =  afwMath.GaussianFunction2D(2.5, 2.5)
            kSize = 6
            k = afwMath.AnalyticKernel(kSize, kSize, kFunc)
            cnvImage = self.ms.Factory(self.ms.getDimensions())
            cnvImage.setXY0(afwImage.PointI(self.ms.getX0(), self.ms.getY0()))
            afwMath.convolve(cnvImage, self.ms, k, True)

            ds = afwDetection.DetectionSetF(cnvImage, threshold, "DETECTED")
            if display:
                ds9.mtv(self.ms, frame=0)
                ds9.mtv(cnvImage, frame=1)
        else:
            ds = afwDetection.DetectionSetF(self.ms, threshold, "DETECTED")
            if display:
                ds9.mtv(self.ms, frame=0)

        objects = ds.getFootprints()
        source = afwDetection.Source()

        for i in range(len(objects)):
            source.setId(i)
            
            algorithms.measureSource(source, self.ms, objects[i], 0.0)

            if display:
                ds9.dot("+", source.getColc() - self.ms.getX0(), source.getRowc() - self.ms.getY0())

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(MeasureTestCase)
    suites += unittest.makeSuite(FindAndMeasureTestCase)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
