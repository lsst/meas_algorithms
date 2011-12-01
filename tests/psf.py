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
Tests for PSF code

Run with:
   python psf.py
or
   python
   >>> import psf; psf.run()
"""

import os, sys
from math import *
import numpy
import unittest
import eups
import lsst.utils.tests as utilsTests
import lsst.pex.exceptions as pexExceptions
import lsst.pex.logging as logging
import lsst.pex.policy as policy
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetection
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.afw.display.ds9 as ds9
import lsst.daf.base as dafBase
import lsst.afw.display.utils as displayUtils
import lsst.meas.algorithms as algorithms
import lsst.meas.algorithms.defects as defects
import lsst.meas.algorithms.utils as maUtils

try:
    type(verbose)
except NameError:
    verbose = 0
    logging.Trace_setVerbosity("algorithms.Interp", verbose)
    display = not False

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def psfVal(ix, iy, x, y, sigma1, sigma2, b):
    return (exp(-0.5*((ix - x)**2 + (iy - y)**2)/sigma1**2) +
            b*exp(-0.5*((ix - x)**2 + (iy - y)**2)/sigma2**2))/(1 + b)

class SpatialModelPsfTestCase(unittest.TestCase):
    """A test case for SpatialModelPsf"""

    @staticmethod
    def measure(objects, exposure):
        """Measure a set of Footprints, returning a sourceList"""
        moPolicy = policy.Policy()

        moPolicy.add("astrometry.SDSS", policy.Policy())
        moPolicy.add("source.astrom",  "SDSS")
        moPolicy.add("source.shape",  "SDSS")

        moPolicy.add("photometry.PSF", policy.Policy())
        moPolicy.add("photometry.NAIVE.radius", 3.0)
        moPolicy.add("source.psfFlux", "PSF")
        moPolicy.add("source.apFlux",  "NAIVE")

        moPolicy.add("shape.SDSS", policy.Policy())

        measureSources = algorithms.makeMeasureSources(exposure, moPolicy)

        sourceList = afwDetection.SourceSet()
        for i, object in enumerate(objects):
            source = afwDetection.Source()
            sourceList.append(source)

            source.setId(i)
            source.setFlagForDetection(source.getFlagForDetection() | algorithms.Flags.BINNED1);
            source.setFootprint(object)

            measureSources.measure(source, exposure)

        return sourceList

    def setUp(self):
        if True:
            width, height = 110, 301
        elif not True:
            width, height = 2*200, 2*300
        else:
            width, height = 50, 5*100
        self.mi = afwImage.MaskedImageF(afwGeom.ExtentI(width, height))
        self.mi.set(0)
        sd = 3                          # standard deviation of image
        self.mi.getVariance().set(sd*sd)
        self.mi.getMask().addMaskPlane("DETECTED")

        self.FWHM = 5
        self.ksize = 35                      # size of desired kernel

        sigma1 = 1.75
        sigma2 = 2*sigma1

        self.exposure = afwImage.makeExposure(self.mi)
        self.exposure.setPsf(afwDetection.createPsf("DoubleGaussian", self.ksize, self.ksize,
                                                    1.5*sigma1, # self.FWHM/(2*sqrt(2*log(2))),
                                                    1, 0.1))

        rand = afwMath.Random()               # make these tests repeatable by setting seed

        im = self.mi.getImage()
        afwMath.randomGaussianImage(im, rand) # N(0, 1)
        im *= sd                              # N(0, sd^2)
        del im

        xarr, yarr = [], []
        if width == 110:
            for x, y in [(20, 20), (60, 20), 
                         (30, 35), (50, 50),
                         (20, 90), (70, 160), (25, 265), (75, 275), (85, 30),
                         (50, 130), (70, 80),
                         (60, 210), (20, 210),
                         ]:
                xarr.append(x)
                yarr.append(y)
        elif width == 2*200:
            for i in range(8*10):
                x = 40*(1 + i%8)
                y = 50*(1 + i//8)
                
                xarr.append(x)
                yarr.append(y)
        else:
            dx, dy = 25, 35
            for i in range((height - self.ksize)//dy):
                x = dx*(1 + 0*abs(i - 2))
                y = dy*(1 + i)

                xarr.append(x)
                yarr.append(y)

        for x, y in zip(xarr, yarr):
            source = afwDetection.Source()

            flux = 10000 - 20*x - 10*(y/float(height))**2
            flux = 10000

            b = 0.2*(1e-2*(y - 0) + 0* 0.5*1e-2*x)

            if True:                    # center offset
                if not True:
                    dx = 0.5
                    dy = 0.5
                else:
                    dx = rand.uniform() - 0.5
                    dy = rand.uniform() - 0.5
            else:
                dx, dy = 0.0, 0.0

            if False:
                psf = afwDetection.createPsf("DoubleGaussian", self.ksize, self.ksize, sigma1, sigma2, b)
                im = psf.computeImage(afwGeom.PointD(0,0), False).convertF()
                im /= im.get(self.ksize//2, self.ksize//2)

                im *= flux
                bbox = afwGeom.BoxI( \
                        afwGeom.PointI(x - self.ksize/2, y - self.ksize/2),
                        afwGeom.ExtentI(self.ksize))
                smi = self.mi.Factory(self.mi, bbox, afwImage.LOCAL)
                im = afwMath.offsetImage(im, dx, dy)

                smi += im
                del psf; del im; del smi
            else:
                totFlux = 0.0
                for iy in range(y - self.ksize//2, y + self.ksize//2 + 1):
                    if iy < 0 or iy >= self.mi.getHeight():
                        continue

                    for ix in range(x - self.ksize//2, x + self.ksize//2 + 1):
                        if ix < 0 or ix >= self.mi.getWidth():
                            continue

                        I = flux*psfVal(ix, iy, x + dx, y + dy, sigma1, sigma2, b)
                        self.mi.getImage().set(ix, iy, self.mi.getImage().get(ix, iy) + rand.poisson(I))
                        self.mi.getVariance().set(ix, iy, self.mi.getVariance().get(ix, iy) + I)

        #
        # Make a kernel with the exactly correct basis functions.  Useful for debugging
        #
        basisKernelList = afwMath.KernelList()
        for sigma in (sigma1, sigma2):
            basisKernel = afwMath.AnalyticKernel(self.ksize, self.ksize,
                                                 afwMath.GaussianFunction2D(sigma, sigma))
            basisImage = afwImage.ImageD(basisKernel.getDimensions())
            basisKernel.computeImage(basisImage, True)
            basisImage /= basisImage.get(self.ksize//2, self.ksize//2)
            if sigma == sigma1:
                basisImage0 = basisImage
            else:
                basisImage -= basisImage0
            basisKernelList.append(afwMath.FixedKernel(basisImage))

        order = 1                                # 1 => up to linear
        spFunc = afwMath.PolynomialFunction2D(order)
        self.exactBasisKernel = afwMath.LinearCombinationKernel(basisKernelList, spFunc)
        self.exactBasisKernel.setSpatialParameters([[1.0] + [0.0]*((order + 1)*(order + 2)//2 - 1)]*2)
        #
        # Create a PSF for measurement
        #
        bbox = afwGeom.BoxI(afwGeom.PointI(0,0), afwGeom.ExtentI(width, height))
        self.cellSet = afwMath.SpatialCellSet(bbox, 100)
        ds = afwDetection.FootprintSetF(self.mi, afwDetection.Threshold(100), "DETECTED")
        self.objects = ds.getFootprints()

        if False and display:
            ds9.mtv(self.mi.getVariance(), title="var")
            import pdb; pdb.set_trace() 
        #
        # Prepare to measure
        #
        self.sourceList = SpatialModelPsfTestCase.measure(self.objects, self.exposure)

        for source in self.sourceList:
            if False and i == 0:
                print "Setting centroids"
                source.setXAstrom(int(source.getXAstrom() + 0.5))
                source.setYAstrom(int(source.getYAstrom() + 0.5))

            try:
                self.cellSet.insertCandidate(algorithms.makePsfCandidate(source, self.mi))
            except Exception, e:
                print e
                continue

    def tearDown(self):
        del self.cellSet
        del self.exposure
        del self.mi
        del self.exactBasisKernel
        del self.objects
        del self.sourceList

    def testPsfDeterminer(self):
        """Test the (PCA) psfDeterminer"""
        secondMomentStarSelectorPolicy = policy.Policy.createPolicy(
            policy.DefaultPolicyFile("meas_algorithms", "policy/secondMomentStarSelectorDictionary.paf"))

        pcaPsfDeterminerPolicy = policy.Policy.createPolicy(
            policy.DefaultPolicyFile("meas_algorithms", "policy/pcaPsfDeterminerDictionary.paf"))
        width, height = self.exposure.getMaskedImage().getDimensions()
        pcaPsfDeterminerPolicy.set("sizeCellX", width)
        pcaPsfDeterminerPolicy.set("sizeCellY", height//3)
        pcaPsfDeterminerPolicy.set("nEigenComponents", 2)
        pcaPsfDeterminerPolicy.set("spatialOrder", 0)
        pcaPsfDeterminerPolicy.set("kernelSizeMin", 31)
        pcaPsfDeterminerPolicy.set("nStarPerCell", 0)
        pcaPsfDeterminerPolicy.set("nStarPerCellSpatialFit", 0) # unlimited
        #
        starSelector = algorithms.makeStarSelector("secondMomentStarSelector", secondMomentStarSelectorPolicy)
        psfCandidateList = starSelector.selectStars(self.exposure, self.sourceList)

        psfDeterminer = algorithms.makePsfDeterminer("pcaPsfDeterminer", pcaPsfDeterminerPolicy)
        
        metadata = dafBase.PropertyList()
        psf, cellSet = psfDeterminer.determinePsf(self.exposure, psfCandidateList, metadata)
        self.exposure.setPsf(psf)
        #
        # Now see well we can subtract the PSF
        #
        self.sourceList = SpatialModelPsfTestCase.measure(self.objects, self.exposure)

        subtracted =  self.mi.Factory(self.mi, True)

        for s in self.sourceList:
            try:
                args = [s.getPsfFlux()] if False else []
                algorithms.subtractPsf(psf, subtracted, s.getXAstrom(), s.getYAstrom(), *args)
            except pexExceptions.LsstCppException, e:
                print e
                import pdb; pdb.set_trace() 
                pass

        if True or display:
            ds9.mtv(subtracted, title="Subtracted", frame=1)
            var = subtracted.getVariance()
            numpy.sqrt(var.getArray(), var.getArray()) # inplace sqrt
            subtracted /= var
            ds9.mtv(subtracted, title="Chi", frame=2)

    def testCandidateList(self):
        self.assertFalse(self.cellSet.getCellList()[0].empty())
        self.assertTrue(self.cellSet.getCellList()[1].empty())
        self.assertFalse(self.cellSet.getCellList()[2].empty())
        self.assertTrue(self.cellSet.getCellList()[3].empty())

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
                width, height = 29, 25
                cand.setWidth(width); cand.setHeight(height);

                im = cand.getImage()
                stamps.append(im)

                self.assertEqual(im.getWidth(), width)
                self.assertEqual(im.getHeight(), height)
        
        if False and display:
            mos = displayUtils.Mosaic()
            mos.makeMosaic(stamps, frame=2)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class psfAttributesTestCase(unittest.TestCase):
    """A test case for psfAttributes"""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testGaussian(self):
        """Check that we can measure a single Gaussian's attributes"""

        sigma0 = 5.0
        aEff0 = 4.0*pi*sigma0**2

        xwid = int(12*sigma0)
        ywid = xwid

        # set the peak of the outer guassian to 0 so this is really a single gaussian.
        psf = afwDetection.createPsf("SingleGaussian", xwid, ywid, sigma0);

        if False and display:
            im = psf.computeImage(afwGeom.PointD(xwid//2, ywid//2))
            ds9.mtv(im, title="N(%g) psf" % sigma0, frame=0)

        psfAttrib = algorithms.PsfAttributes(psf, xwid//2, ywid//2)
        sigma = psfAttrib.computeGaussianWidth(psfAttrib.ADAPTIVE_MOMENT)
        m1    = psfAttrib.computeGaussianWidth(psfAttrib.FIRST_MOMENT)
        m2    = psfAttrib.computeGaussianWidth(psfAttrib.SECOND_MOMENT)
        noise = psfAttrib.computeGaussianWidth(psfAttrib.NOISE_EQUIVALENT)
        bick  = psfAttrib.computeGaussianWidth(psfAttrib.BICKERTON)
        aEff  = psfAttrib.computeEffectiveArea();

        if verbose:
            print "Adaptive            %g v %g" % (sigma0, sigma)
            print "First moment        %g v %g" % (sigma0, m1)
            print "Second moment       %g v %g" % (sigma0, m2)
            print "Noise Equivalent    %g v %g" % (sigma0, sigma)
            print "Bickerton           %g v %g" % (sigma0, bick)
            print "Effective area      %g v %f" % (aEff0, aEff)

        self.assertTrue(abs(sigma0 - sigma) <= 1.0e-2)
        self.assertTrue(abs(sigma0 - m1) <= 3.0e-2)
        self.assertTrue(abs(sigma0 - m2) <= 1.0e-2)
        self.assertTrue(abs(sigma0 - noise) <= 1.0e-2)
        self.assertTrue(abs(sigma0 - bick) <= 1.0e-2)
        self.assertTrue(abs(aEff0 - aEff) <= 1.0e-2)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(SpatialModelPsfTestCase)
    suites += unittest.makeSuite(psfAttributesTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the utilsTests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)











