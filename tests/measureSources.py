#!/usr/bin/env python
"""
Tests for measuring things

Run with:
   python MeasureSources.py
or
   python
   >>> import MeasureSources; MeasureSources.run()
"""

import math, os, sys, unittest
import lsst.utils.tests as tests
import lsst.pex.exceptions as pexExceptions
import lsst.pex.logging as pexLogging
import lsst.pex.policy as pexPolicy
import lsst.afw.detection as afwDetection
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.meas.algorithms as measAlg

try:
    type(verbose)
except NameError:
    verbose = 0
pexLogging.Trace_setVerbosity("afwDetection.Measure", verbose)

try:
    type(display)
except NameError:
    display = False

import lsst.afw.display.ds9 as ds9

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class MeasureSourcesTestCase(unittest.TestCase):
    """A test case for Measure"""

    def setUp(self):
        pass
        
    def tearDown(self):
        pass

    def testNaiveMeasure(self):
        mi = afwImage.MaskedImageF(afwGeom.ExtentI(100, 200))
        mi.set(10)
        #
        # Create our measuring engine
        #
        exp = afwImage.makeExposure(mi)
        
        control = measAlg.NaiveFluxControl()
        control.radius = 10.0

        mp = measAlg.MeasureSources()
        mp.addAlgorithm(control)
        table = afwTable.SourceTable.make(mp.getSchema())
        source = table.makeRecord()
        mp.apply(source, exp, afwGeom.Point2D(30, 50))
        flux = 3170.0
        self.assertEqual(source.get(control.name), flux)

    def testApertureMeasure(self):
        mi = afwImage.MaskedImageF(afwGeom.ExtentI(100, 200))
        mi.set(10)
        #
        # Create our measuring engine
        #

        radii =  ( 1.0,   5.0,   10.0)  # radii to use
        fluxes = [50.0, 810.0, 3170.0]  # corresponding correct fluxes
        
        control = measAlg.ApertureFluxControl()
        control.radii = radii
        
        exp = afwImage.makeExposure(mi)
        mp = measAlg.MeasureSources()
        mp.addAlgorithm(control)

        table = afwTable.SourceTable.make(mp.getSchema())
        source = table.makeRecord()

        mp.apply(source, exp, afwGeom.Point2D(30, 50))
        measured = source[control.name]
        for i, f in enumerate(fluxes):
            self.assertEqual(f, measured[i])

    def testEllipticalGaussian(self):
        """Test measuring the properties of an elliptical Gaussian"""

        width, height = 200, 200
        xcen, ycen = 0.5*width, 0.5*height
        #
        # Make the object
        #
        gal = afwImage.ImageF(afwGeom.ExtentI(width, height))
        a, b, theta = float(10), float(5), 20
        flux = 1e4
        I0 = flux/(2*math.pi*a*b)

        c, s = math.cos(math.radians(theta)), math.sin(math.radians(theta))
        for y in range(height):
            for x in range(width):
                dx, dy = x - xcen, y - ycen
                u =  c*dx + s*dy
                v = -s*dx + c*dy
                val = I0*math.exp(-0.5*((u/a)**2 + (v/b)**2))
                if val < 0:
                    val = 0
                gal.set(x, y, val)

        objImg = afwImage.makeExposure(afwImage.makeMaskedImage(gal))
        objImg.getMaskedImage().getVariance().set(1.0)
        del gal

        if display:
            frame = 0
            ds9.mtv(objImg, frame=frame, title="Elliptical")

        self.assertAlmostEqual(1.0, afwMath.makeStatistics(objImg.getMaskedImage().getImage(),
                                                           afwMath.SUM).getValue()/flux)
        #
        # Now measure some annuli
        #
        
        center = afwGeom.Point2D(xcen, ycen)

        for r1, r2 in [(0,      0.45*a),
                       (0.45*a, 1.0*a),
                       ( 1.0*a, 2.0*a),
                       ( 2.0*a, 3.0*a),
                       ( 3.0*a, 5.0*a),
                       ( 3.0*a, 10.0*a),
                       ]:
            control = measAlg.SincFluxControl()
            control.radius1 = r1
            control.radius2 = r2
            control.angle = math.radians(theta)
            control.ellipticity = 1 - b/a
            mp = measAlg.MeasureSources()
            mp.addAlgorithm(control)
            table = afwTable.SourceTable.make(mp.getSchema())
            source = table.makeRecord()

            if display:                 # draw the inner and outer boundaries of the aperture
                Mxx = 1
                Myy = (b/a)**2

                mxx, mxy, myy = c**2*Mxx + s**2*Myy, c*s*(Mxx - Myy), s**2*Mxx + c**2*Myy
                for r in (r1, r2):
                    ds9.dot("@:%g,%g,%g" % (r**2*mxx, r**2*mxy, r**2*myy), xcen, ycen, frame=frame)

            mp.apply(source, objImg, center)

            self.assertAlmostEqual(math.exp(-0.5*(r1/a)**2) - math.exp(-0.5*(r2/a)**2),
                                   source.get(control.name)/flux, 5)

        control = measAlg.GaussianFluxControl()
        mp = measAlg.MeasureSources()
        mp.addAlgorithm(control)
        table = afwTable.SourceTable.make(mp.getSchema())
        source = table.makeRecord()
        mp.apply(source, objImg, center)
        # we haven't provided a PSF, so the built-in aperture correction won't work...but we'll get
        # a result anyway
        self.assertEqual(source.get(control.name + ".flags"), True)
        self.assertEqual(source.get(control.name + ".flags.badapcorr"), True)
        gflux = source.get(control.name)
        err = gflux/flux - 1
        if abs(err) > 1.5e-5:
            self.assertEqual(gflux, flux, ("%g, %g: error is %g" % (gflux, flux, err)))

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(MeasureSourcesTestCase)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
