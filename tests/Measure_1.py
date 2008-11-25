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
import unittest
import lsst.utils.tests as tests
import lsst.pex.logging as logging
import lsst.afw.image as imageLib
import lsst.afw.detection as afwDetection
import lsst.detection.detectionLib as detection

try:
    type(verbose)
except NameError:
    verbose = 0
logging.Trace_setVerbosity("detection.Measure", verbose)

try:
    type(display)
except NameError:
    display = False

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
        self.ms = imageLib.MaskedImageF(12, 8)
        im = self.ms.getImage()
        #
        # Objects that we should detect
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

        xcentroid = [4.0, 8.4, 3.0]
        ycentroid = [2.0, 5.4, 6.0]
        flux = [50.0, 100.0, 20.0]
        
        ds = afwDetection.DetectionSetF(self.ms, afwDetection.Threshold(10), "DETECTED")

        if display:
            import lsst.afw.display.ds9 as ds9
            ds9.mtv(self.ms, frame=0)

        objects = ds.getFootprints()
        measure = detection.MeasureF(self.ms)
        diaptr = afwDetection.Source()

        for i in range(len(objects)):
            diaptr.setId(i)
            measure.measureSource(diaptr, objects[i], 0.0)

            self.assertAlmostEqual(diaptr.getColc(), xcentroid[i] + 0.5, 6)
            self.assertAlmostEqual(diaptr.getRowc(), ycentroid[i] + 0.5, 6)
            self.assertEqual(diaptr.getFlux(), flux[i])

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(MeasureTestCase)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
