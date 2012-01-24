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
Tests for Footprints, FootprintSets, and Measure

Run with:
   python Measure_1.py
or
   python
   >>> import Measure_1; Measure_1.run()
"""

import os, sys, unittest
import math; from math import *
import eups
import lsst.utils.tests as tests
import lsst.pex.logging as logging
import lsst.pex.policy as policy
import lsst.afw.detection as afwDetection
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
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

        def insert(self, im, dx = 0, dy = 0):
            """Insert self into an image"""
            for sp in self.spans:
                y, x0, x1 = sp
                for x in range(x0, x1 + 1):
                    im.set(x + dx, y + dy, self.val)

        def __eq__(self, other):
            for osp, sp in zip(other.getSpans(), self.spans):
                if osp.toString() != toString(sp):
                    return False
                
            return True
    
    def setUp(self):
        ms = afwImage.MaskedImageF(afwGeom.ExtentI(31, 27))
        var = ms.getVariance(); var.set(1); del var
        bbox = afwGeom.BoxI(afwGeom.PointI(1,1), afwGeom.ExtentI(24, 20))
        self.mi = afwImage.MaskedImageF(ms, bbox, afwImage.LOCAL)
        self.exposure = afwImage.makeExposure(self.mi)
        im = self.mi.getImage()
        #
        # Objects that we should detect.  These are coordinates in the subimage
        #
        self.objects = []
        self.objects += [self.Object(10, [(1, 4, 4), (2, 3, 5), (3, 4, 4)])]
        self.objects += [self.Object(20, [(5, 7, 8), (5, 10, 10), (6, 8, 9)])]
        self.objects += [self.Object(20, [(8, 3, 3)])]

        im.set(0)                       # clear image
        for obj in self.objects:
            obj.insert(im, 5, 5)
        #
        # Add a few more pixels to make peaks that we can centroid around
        #
        for x, y in [(9, 7), (13, 11)]:
            im.set(x, y, 1 + im.get(x, y))
        
    def tearDown(self):
        del self.mi
        del self.exposure

    def testFootprintsMeasure(self):
        """Check that we can measure the objects in a detectionSet"""

        xcentroid = [10.0, 14.0,        9.0]
        ycentroid = [8.0, 11.5061728,  14.0]
        flux = [51.0, 101.0,         20.0]
        wflux = [51.0, 101.0,        20.0]
        
        ds = afwDetection.FootprintSetF(self.mi, afwDetection.Threshold(10), "DETECTED")

        if display:
            ds9.mtv(self.mi, frame=0)
            ds9.mtv(self.mi.getVariance(), frame=1)

        objects = ds.getFootprints()
        source = afwDetection.Source(0, afwDetection.Footprint())

        msPolicy = policy.Policy.createPolicy(policy.PolicyString(
        """#<?cfg paf policy?>
        astrometry: {
            GAUSSIAN: {
                enabled: false
            }
            NAIVE: {
                enabled: true
            }
            SDSS: {
                enabled: false
            }
        }
        photometry: {
            PSF: {
                enabled: true
            }
            NAIVE: {
                radius: 3.0
            }
            SINC: {
                enabled: false
            }
        }
        shape: {
            SDSS: {
                enabled: true
            }
        }
        source: {
            astrom: "NAIVE"
            psfFlux: "PSF"
            apFlux: "NAIVE"
        }"""))

        sigma = 1e-10; psf = afwDetection.createPsf("DoubleGaussian", 11, 11, sigma) # i.e. a single pixel
        self.exposure.setPsf(psf)

        measureSources = algorithms.makeMeasureSources(self.exposure, msPolicy)

        for i in range(len(objects)):
            source.setId(i)
            source.setFootprint(objects[i])

            measureSources.measure(source, self.exposure)

            xc, yc = source.getXAstrom() - self.mi.getX0(), source.getYAstrom() - self.mi.getY0()
            if display:
                ds9.dot("+", xc, yc)

            self.assertAlmostEqual(source.getXAstrom(), xcentroid[i], 6)
            self.assertAlmostEqual(source.getYAstrom(), ycentroid[i], 6)
            self.assertEqual(source.getApFlux(), flux[i])
            self.assertAlmostEqual(source.getApFluxErr(), math.sqrt(29), 6) # 29 pixels in 3pixel circular ap.
            # We're using a delta-function PSF, so the psfFlux should be the pixel under the centroid,
            # iff the object's centred in the pixel
            if xc == int(xc) and yc == int(yc):
                self.assertAlmostEqual(source.getPsfFlux(),
                                       self.exposure.getMaskedImage().getImage().get(int(xc + 0.5),
                                                                                 int(yc + 0.5)))
                self.assertAlmostEqual(source.getPsfFluxErr(),
                                       self.exposure.getMaskedImage().getVariance().get(int(xc + 0.5),
                                                                                    int(yc + 0.5)))
            
            
class FindAndMeasureTestCase(unittest.TestCase):
    """A test case detecting and measuring objects"""
    def setUp(self):
        self.mi = afwImage.MaskedImageF(os.path.join(eups.productDir("afwdata"),
                                                     "CFHT", "D4", "cal-53535-i-797722_1"))

        self.FWHM = 5
        self.psf = afwDetection.createPsf("DoubleGaussian", 15, 15, self.FWHM/(2*sqrt(2*log(2))))

        if False:                       # use full image, trimmed to data section
            self.XY0 = afwGeom.PointI(32, 2)
            self.mi = self.mi.Factory(self.mi, afwGeom.BoxI(self.XY0, afwGeom.PointI(2079, 4609)), afwImage.LOCAL)
            self.mi.setXY0(afwGeom.PointI(0, 0))
        else:                           # use sub-image
            self.XY0 = afwGeom.PointI(824, 140)
            self.mi = self.mi.Factory(self.mi, afwGeom.BoxI(self.XY0, afwGeom.ExtentI(256, 256)), afwImage.LOCAL)

        self.mi.getMask().addMaskPlane("DETECTED")
        self.exposure = afwImage.makeExposure(self.mi)

    def tearDown(self):
        del self.mi
        del self.psf
        del self.exposure

    def testDetection(self):
        """Test object detection"""
        #
        # Fix defects
        #
        # Mask known bad pixels
        #
        badPixels = defects.policyToBadRegionList(os.path.join(eups.productDir("meas_algorithms"),
                                                               "policy/BadPixels.paf"))
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
        bgGridSize = 64  # was 256 ... but that gives only one region and the spline breaks
        bctrl = afwMath.BackgroundControl(afwMath.Interpolate.NATURAL_SPLINE);
        bctrl.setNxSample(int(self.mi.getWidth()/bgGridSize) + 1);
        bctrl.setNySample(int(self.mi.getHeight()/bgGridSize) + 1);
        backobj = afwMath.makeBackground(self.mi.getImage(), bctrl)
        
        img = self.mi.getImage(); img -= backobj.getImageF(); del img
        #
        # Remove CRs
        #
        crPolicy = policy.Policy.createPolicy(os.path.join(eups.productDir("meas_algorithms"),
                                                           "policy",
                                                           "CrRejectDictionary.paf"))
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
        psf = afwDetection.createPsf("DoubleGaussian", 15, 15, self.FWHM/(2*sqrt(2*log(2))))

        cnvImage = self.mi.Factory(self.mi.getBBox(afwImage.PARENT))
        kernel = psf.getKernel()
        afwMath.convolve(cnvImage, self.mi, kernel, afwMath.ConvolutionControl())

        msk = cnvImage.getMask(); msk |= savedMask; del msk # restore the saved bits

        threshold = afwDetection.Threshold(3, afwDetection.Threshold.STDEV)
        #
        # Only search the part of the frame that was PSF-smoothed
        #        
        llc = afwGeom.PointI(psf.getKernel().getWidth()/2, psf.getKernel().getHeight()/2)
        urc = afwGeom.PointI(cnvImage.getWidth() -llc[0] - 1, cnvImage.getHeight() - llc[1] - 1)
        middle = cnvImage.Factory(cnvImage, afwGeom.BoxI(llc, urc), afwImage.LOCAL)
        ds = afwDetection.FootprintSetF(middle, threshold, "DETECTED")
        del middle
        #
        # Reinstate the saved (e.g. BAD) (and also the DETECTED | EDGE) bits in the unsmoothed image
        #
        savedMask <<= cnvImage.getMask()
        msk = self.mi.getMask(); msk |= savedMask; del msk
        del savedMask

        if display:
            ds9.mtv(self.mi, frame = 0)
            ds9.mtv(cnvImage, frame = 1)

        objects = ds.getFootprints()
        #
        # Time to actually measure
        #
        msPolicy = policy.Policy.createPolicy(policy.DefaultPolicyFile("meas_algorithms",
            "tests/MeasureSources.paf"))
        msPolicy = msPolicy.getPolicy("measureSources")
        measureSources = algorithms.makeMeasureSources(self.exposure, msPolicy)

        sourceList = afwDetection.SourceSet()
        for i in range(len(objects)):
            source = afwDetection.Source(i, objects[i])
            sourceList.append(source)

            source.setFlagForDetection(source.getFlagForDetection() | algorithms.Flags.BINNED1);

            try:
                measureSources.apply(source)
            except Exception, e:
                print "RHL", e

            if source.getFlagForDetection() & algorithms.Flags.EDGE:
                continue

            if display:
                ds9.dot("+", source.getXAstrom() - self.mi.getX0(), source.getYAstrom() - self.mi.getY0())

class GaussianPsfTestCase(unittest.TestCase):
    """A test case detecting and measuring Gaussian PSFs"""
    def setUp(self):
        FWHM = 5
        psf = afwDetection.createPsf("DoubleGaussian", 15, 15, FWHM/(2*sqrt(2*log(2))))
        mi = afwImage.MaskedImageF(afwGeom.ExtentI(100, 100))

        self.xc, self.yc, self.flux = 45, 55, 1000.0
        mi.getImage().set(self.xc, self.yc, self.flux)

        cnvImage = mi.Factory(mi.getDimensions())
        afwMath.convolve(cnvImage, mi, psf.getKernel(), afwMath.ConvolutionControl())

        self.exp = afwImage.makeExposure(cnvImage)
        self.exp.setPsf(psf)

        if display and False:
            ds9.mtv(self.exp)

    def tearDown(self):
        del self.exp

    def testPsfFlux(self):
        """Test that fluxes are measured correctly"""
        #
        # Total flux in image
        #
        flux = afwMath.makeStatistics(self.exp.getMaskedImage(), afwMath.SUM).getValue()
        self.assertAlmostEqual(flux/self.flux, 1.0)

        #
        # Various algorithms
        #
        photoAlgorithms = ["NAIVE", "PSF", "SINC",]
        mp = algorithms.makeMeasurePhotometry(self.exp)
        for a in photoAlgorithms:
            mp.addAlgorithm(a)

        rad = 10.0
        pol = policy.Policy(policy.PolicyString(
            """#<?cfg paf policy?>
            NAIVE.radius: %f
            SINC.radius:  %f
            """ % (rad, rad)
            ))
        mp.configure(pol)        

        source = afwDetection.Source(0, afwDetection.Footprint())

        for a in photoAlgorithms:
            photom = mp.measure(source, self.exp, afwGeom.Point2D(self.xc, self.yc)).find(a)
            self.assertAlmostEqual(photom.getFlux()/self.flux, 1.0, 4, "Measuring with %s: %g v. %g" %
                                   (a, photom.getFlux(), self.flux))

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(MeasureTestCase)
    if eups.productDir("afwdata"):
        suites += unittest.makeSuite(FindAndMeasureTestCase)
    else:
        print >> sys.stderr, "You must set up afwdata to run the CFHT-based tests"
    suites += unittest.makeSuite(GaussianPsfTestCase)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
