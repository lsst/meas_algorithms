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
import lsst.afw.display.utils as displayUtils
import lsst.meas.algorithms as algorithms
import lsst.meas.algorithms.defects as defects
import lsst.meas.algorithms.utils as maUtils

try:
    type(verbose)
except NameError:
    verbose = 0
    logging.Trace_setVerbosity("algorithms.Interp", verbose)
    display = False

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def psfVal(ix, iy, x, y, sigma1, sigma2, b):
    return (exp(-0.5*((ix - x)**2 + (iy - y)**2)/sigma1**2) +
            b*exp(-0.5*((ix - x)**2 + (iy - y)**2)/sigma2**2))/(1 + b)

class SpatialModelPsfTestCase(unittest.TestCase):
    """A test case for SpatialModelPsf"""

    def setUp(self):
        if True:
            width, height = 100, 301
        elif not True:
            width, height = 2*200, 2*300
        else:
            width, height = 50, 3*100
        self.mi = afwImage.MaskedImageF(width, height)
        self.mi.set(0)
        sd = 3                          # standard deviation of image
        self.mi.getVariance().set(sd*sd)
        self.mi.getMask().addMaskPlane("DETECTED")

        self.FWHM = 5
        self.ksize = 35                      # size of desired kernel

        self.exposure = afwImage.makeExposure(self.mi)
        self.exposure.setPsf(afwDetection.createPsf("DoubleGaussian", self.ksize, self.ksize,
                                                    self.FWHM/(2*sqrt(2*log(2))), 1, 0.1))

        rand = afwMath.Random()               # make these tests repeatable by setting seed

        im = self.mi.getImage()
        afwMath.randomGaussianImage(im, rand) # N(0, 1)
        im *= sd                              # N(0, sd^2)
        del im

        sigma1 = 1.5
        sigma2 = 2*sigma1

        xarr, yarr = [], []
        if width == 100:
            for x, y in [(20, 20), (60, 20),
                         (30, 35), (50, 50),
                         (50, 130), (70, 80),
                         (60, 210), (20, 210)]:
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
            for i in range(height//dy):
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
                im = psf.computeImage(afwGeom.makePointD(0,0), False).convertF()
                im /= im.get(self.ksize//2, self.ksize//2)

                im *= flux
                smi = self.mi.Factory(self.mi, afwImage.BBox(afwImage.PointI(x - self.ksize/2, y - self.ksize/2),
                                                             self.ksize, self.ksize))

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
                        self.mi.getImage().set(ix, iy, self.mi.getImage().get(ix, iy) +
                                               flux*psfVal(ix, iy, x + dx, y + dy, sigma1, sigma2, b))

                        totFlux += flux*psfVal(ix, iy, x + dx, y + dy, sigma1, sigma2, b)

                print "RHL", x + dx, y + dy, totFlux

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
            basisKernelList.append(afwMath.FixedKernel(basisImage))

        order = 1                                # 1 => up to linear
        spFunc = afwMath.PolynomialFunction2D(order)
        self.exactBasisKernel = afwMath.LinearCombinationKernel(basisKernelList, spFunc)
        self.exactBasisKernel.setSpatialParameters([[1.0] + [0.0]*((order + 1)*(order + 2)//2 - 1)]*2)
        #
        # Create a PSF for measurement
        #
        psf = afwDetection.createPsf("DoubleGaussian", self.ksize, self.ksize,
                                   self.FWHM/(2*sqrt(2*log(2))), 1, 0.1)

        self.cellSet = afwMath.SpatialCellSet(afwImage.BBox(afwImage.PointI(0, 0), width, height), 100)
        ds = afwDetection.FootprintSetF(self.mi, afwDetection.Threshold(100), "DETECTED")
        objects = ds.getFootprints()

        if False and display:
            ds9.mtv(self.mi.getVariance(), title="var"); 0/0
        #
        # Prepare to measure
        #
        moPolicy = policy.Policy()

        moPolicy.add("astrometry.SDSS", policy.Policy())
        moPolicy.add("source.astrom",  "SDSS")
        moPolicy.add("source.shape",  "SDSS")

        moPolicy.add("photometry.PSF", policy.Policy())
        moPolicy.add("photometry.NAIVE.radius", 3.0)
        moPolicy.add("source.psfFlux", "PSF")
        moPolicy.add("source.apFlux",  "NAIVE")

        moPolicy.add("shape.SDSS", policy.Policy())

        self.exposure.setPsf(psf)
        measureSources = algorithms.makeMeasureSources(self.exposure, moPolicy)

        sourceList = afwDetection.SourceSet()
        for i in range(len(objects)):
            source = afwDetection.Source()
            sourceList.append(source)

            source.setId(i)
            source.setFlagForDetection(source.getFlagForDetection() | algorithms.Flags.BINNED1);

            measureSources.apply(source, objects[i])
            if False and i == 0:
                print "Setting centroids"
                source.setXAstrom(int(source.getXAstrom() + 0.5))
                source.setYAstrom(int(source.getYAstrom() + 0.5))

            self.cellSet.insertCandidate(algorithms.makePsfCandidate(source, self.mi))

    def tearDown(self):
        del self.cellSet
        del self.exposure
        del self.mi
        del self.exactBasisKernel

    def testGetPcaKernel(self):
        """Convert our cellSet to a LinearCombinationKernel"""

        nEigenComponents = 2
        spatialOrder  =    1
        kernelSize =       self.ksize
        nStarPerCell =     4
        nStarPerCellSpatialFit = 0
        tolerance =     1e-5
        reducedChi2ForPsfCandidates = 2.0
        nIterForPsf =      5

        width, height = kernelSize, kernelSize
        algorithms.PsfCandidateF.setWidth(width); algorithms.PsfCandidateF.setHeight(height);
        nu = width*height - 1           # number of degrees of freedom/star for chi^2

        reply = ""
        for iter in range(nIterForPsf):
            if display:
                ds9.mtv(self.mi, title="Input image", frame=0)
                #
                # Show the candidates we're using
                #
                for cell in self.cellSet.getCellList():
                    #print "Cell", cell.getBBox()
                    i = 0
                    for cand in cell.begin(False): # don't skip BAD stars
                        cand = algorithms.cast_PsfCandidateF(cand)

                        i += 1
                        source = cand.getSource()

                        xc, yc = source.getXAstrom() - self.mi.getX0(), source.getYAstrom() - self.mi.getY0()

                        if cand.isBad():
                            ds9.dot("+", xc, yc, ctype = ds9.RED)
                        elif i <= nStarPerCell:
                            ds9.dot("+", xc, yc, ctype = ds9.GREEN)
                        else:
                            ds9.dot("+", xc, yc, ctype = ds9.YELLOW)

            pair = algorithms.createKernelFromPsfCandidates(self.cellSet, nEigenComponents, spatialOrder,
                                                            kernelSize, nStarPerCell)

            kernel, eigenValues = pair[0], pair[1]; del pair

            print "lambda", " ".join(["%g" % l for l in eigenValues])

            if not False:                    # fake the input kernel.  Debugging ONLY
                print "Using exact (input) Kernel"
                kernel = self.exactBasisKernel

            pair = algorithms.fitSpatialKernelFromPsfCandidates(kernel, self.cellSet, False,
                                                                nStarPerCellSpatialFit, tolerance)
            status, chi2 = pair[0], pair[1]; del pair
            print "Spatial fit: status = %s,  chi^2 = %.2g" % (status, chi2)

            psf = afwDetection.createPsf("PCA", kernel) # Hurrah!
            #
            # Label PSF candidate stars with bad chi^2 as BAD
            #
            nDiscard = 1
            for cell in self.cellSet.getCellList():
                worstId, worstChi2 = -1, -1
                for cand in cell.begin(True): # only not BAD candidates
                    cand = algorithms.cast_PsfCandidateF(cand)

                    rchi2 = cand.getChi2()/nu

                    if rchi2 < reducedChi2ForPsfCandidates:
                        cand.setStatus(afwMath.SpatialCellCandidate.GOOD)
                        continue

                    if rchi2 > worstChi2:
                        worstId, worstChi2 = cand.getId(), rchi2
                        
                for cand in cell.begin(True): # only not BAD candidates
                    cand = algorithms.cast_PsfCandidateF(cand)
                    if cand.getId() == worstId:
                        cand.setStatus(afwMath.SpatialCellCandidate.BAD)

            self.assertTrue(afwMath.cast_AnalyticKernel(psf.getKernel()) is None)
            self.assertTrue(afwMath.cast_LinearCombinationKernel(psf.getKernel()) is not None)
            #
            # OK, we're done for this iteration.  The rest is fluff
            #
            if not display:
                continue
            
            #print psf.getKernel().toString()

            eImages = []
            for k in afwMath.cast_LinearCombinationKernel(psf.getKernel()).getKernelList():
                im = afwImage.ImageD(k.getDimensions())
                k.computeImage(im, False)
                eImages.append(im)

            mos = displayUtils.Mosaic()
            frame = 3
            if False:
                ds9.mtv(mos.makeMosaic(eImages), title="Eigen Images", frame=frame)
            #
            # Make a mosaic of PSF candidates
            #
            stamps = []; stampInfo = []

            for cell in self.cellSet.getCellList():
                for cand in cell.begin(False):
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

            if False:
                mos = displayUtils.Mosaic()
            else:
                mos = None

            frame = 1
            if mos:
                ds9.mtv(mos.makeMosaic(stamps), title="Stamps", frame=frame, lowOrderBits=True)
                mos.drawLabels(stampInfo, frame=frame)

            if not False:
                maUtils.showPsfCandidates(self.exposure, self.cellSet, psf, frame=frame, normalize=False)
            #
            # Reconstruct the PSF as a function of position
            #
            psfImages = []; labels = []

            nx, ny = 3, 4
            for iy in range(ny):
                for ix in range(nx):
                    x = int((ix + 0.5)*self.mi.getWidth()/nx)
                    y = int((iy + 0.5)*self.mi.getHeight()/ny)

                    im = psf.computeImage(afwGeom.makePointD(x, y))
                    psfImages.append(im.Factory(im, True))
                    labels.append("PSF(%d,%d)" % (int(x), int(y)))

                    if not True:
                        print x, y, "PSF parameters:", psf.getKernel().getKernelParameters()

            frame = 2
            if mos:
                mos.makeMosaic(psfImages, frame = frame, mode = nx)
                mos.drawLabels(labels, frame = frame)

            stamps = []; stampInfo = []

            for cell in self.cellSet.getCellList():
                for cand in cell.begin(False): # include bad candidates
                    cand = algorithms.cast_PsfCandidateF(cand)

                    infoStr = "%d X^2=%.1f" % (cand.getSource().getId(), cand.getChi2()/nu)

                    if cand.isBad():
                        if True:
                            infoStr += "B"
                        else:
                            continue

                    im = cand.getImage()
                    stamps.append(im)
                    stampInfo.append(infoStr)

            if mos:
                try:
                    frame = 5
                    mos.makeMosaic(stamps, frame = frame, title="Psf Candidates")
                    mos.drawLabels(stampInfo, frame = frame)
                except RuntimeError, e:
                    print e

            residuals = self.mi.Factory(self.mi, True)
            for cell in self.cellSet.getCellList():
                for cand in cell.begin(False):
                    #
                    # Swig doesn't know that we inherited from SpatialCellImageCandidate;  all
                    # it knows is that we have a SpatialCellCandidate, and SpatialCellCandidates
                    # don't know about getImage;  so cast the pointer to PsfCandidate
                    #
                    cand = algorithms.cast_PsfCandidateF(cand)
                    s = cand.getSource()

                    algorithms.subtractPsf(psf, residuals, s.getXAstrom(), s.getYAstrom())

            ds9.mtv(residuals, title="Residuals", frame=4)

            if iter < nIterForPsf - 1 and reply != "c":
                while True:
                    try:
                        reply = raw_input("Next iteration? [ync] ")
                    except EOFError:
                        reply = "n"
                        
                    if reply in ("", "c", "n", "y"):
                        break
                    else:
                        print >> sys.stderr, "Unrecognised response: %s" % reply

                if reply == "n":
                    break
            
    def testCandidateList(self):
        if False and display:
            ds9.mtv(self.mi)

            for cell in self.cellSet.getCellList():
                x0, y0 = cell.getBBox().getX0(), cell.getBBox().getY0()
                x1, y1 = cell.getBBox().getX1(), cell.getBBox().getY1()
                
                print x0, y0, " ", x1, y1
                x0 -= 0.5; y0 -= 0.5
                x1 += 0.5; y1 += 0.5

                ds9.line([(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)], ctype = ds9.RED)

        self.assertFalse(self.cellSet.getCellList()[0].empty())
        self.assertFalse(self.cellSet.getCellList()[1].empty())
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
        
        if display:
            mos = displayUtils.Mosaic()
            mos.makeMosaic(stamps, frame = 1)

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
            im = psf.computeImage(afwGeom.makePointD(xwid//2, ywid//2))
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

class RHLTestCase(unittest.TestCase):
    """A test case for SpatialModelPsf"""

    def calcDoubleGaussian(self, im, x, y, amp, sigma1, sigma2 = 1.0, b = 0):
        """Insert a DoubleGaussian into the image centered at (x, y)"""
        import math

        x = x - im.getX0(); y = y - im.getY0()

        for ix in range(im.getWidth()):
            for iy in range(im.getHeight()):
                r2 = math.pow(x - ix, 2) + math.pow(y - iy, 2)
                val = math.exp(-r2/(2.0*pow(sigma1, 2))) + b*math.exp(-r2/(2.0*pow(sigma2, 2)))
                im.set(ix, iy, amp/(1 + b)*val)

    def setUp(self):
        width, height = 300, 250
        self.mi = afwImage.MaskedImageF(width, height)
        self.mi.set(0)
        self.mi.getVariance().set(10)
        self.mi.getMask().addMaskPlane("DETECTED")

        self.ksize = 45                      # size of desired kernel

        for x, y in [(120, 120), (160, 120), ]:
            flux = 10000 # - 0*x - 10*(y - 10)

            sigma = 3
            dx, dy = 0.50, 0.50
            dx, dy = -0.50, -0.50
            #dx, dy = 0, 0

            smi = self.mi.getImage().Factory(self.mi.getImage(),
                                             afwImage.BBox(afwImage.PointI(x - self.ksize/2,
                                                                           y - self.ksize/2),
                                                           self.ksize, self.ksize))
            
            im = afwImage.ImageF(self.ksize, self.ksize)
            self.calcDoubleGaussian(im, self.ksize/2 + dx, self.ksize/2 + dy, 1.0, sigma, 1, 0.1)

            #im /= afwMath.makeStatistics(im, afwMath.MEAN).getValue()*im.getHeight()*im.getWidth()
            im /= afwMath.makeStatistics(im, afwMath.MAX).getValue()
            im *= flux

            smi += im
            del im; del smi

        self.FWHM = 5
        psf = afwDetection.createPsf("DoubleGaussian", self.ksize, self.ksize,
                                   self.FWHM/(2*sqrt(2*log(2))), 1, 0.1)

        self.cellSet = afwMath.SpatialCellSet(afwImage.BBox(afwImage.PointI(0, 0), width, height), 100)
        ds = afwDetection.FootprintSetF(self.mi, afwDetection.Threshold(10), "DETECTED")
        objects = ds.getFootprints()
        #
        # Prepare to measure
        #
        moPolicy = policy.Policy()
        moPolicy.add("centroidAlgorithm", "SDSS")
        moPolicy.add("shapeAlgorithm", "SDSS")
        moPolicy.add("photometryAlgorithm", "NAIVE")
        moPolicy.add("apRadius", 3.0)
 
        measureSources = algorithms.makeMeasureSources(afwImage.makeExposure(self.mi), moPolicy)

        sourceList = afwDetection.SourceSet()
        for i in range(len(objects)):
            source = afwDetection.Source()
            sourceList.append(source)

            source.setId(i)

            measureSources.apply(source, objects[i])
            source.setFlagForDetection(source.getFlagForDetection() | algorithms.Flags.BINNED1);

            source.setXAstrom(1e-2*int(100*source.getXAstrom() + 0.5)) # get exact centroids
            source.setYAstrom(1e-2*int(100*source.getYAstrom() + 0.5))

            if not False:
                print source.getXAstrom(), source.getYAstrom(), source.getPsfFlux(), \
                      maUtils.explainDetectionFlags(source.getFlagForDetection())

            self.cellSet.insertCandidate(algorithms.makePsfCandidate(source, self.mi))
            
        frame = 1
        ds9.mtv(self.mi, frame=frame, title="Double Gaussian")

    def tearDown(self):
        del self.cellSet
        del self.mi

    def testRHL(self):
        """Convert our cellSet to a LinearCombinationKernel"""

        nEigenComponents = 1
        spatialOrder  =    0
        kernelSize =      35
        nStarPerCell =     2

        width, height = 45, 45
        algorithms.PsfCandidateF.setWidth(width); algorithms.PsfCandidateF.setHeight(height);
        #
        # Show candidates
        #
        if not False:
            stamps = []
            for cell in self.cellSet.getCellList():
                for cand in cell:
                    cand = algorithms.cast_PsfCandidateF(cand)
                    s = cand.getSource()

                    im = cand.getImage()

                    stamps.append(im)

            if False:
                mos = displayUtils.Mosaic()
                frame = 2
                im = mos.makeMosaic(stamps)
                
                imim = im.getImage()
                imim *= 10000/afwMath.makeStatistics(imim, afwMath.MAX).getValue()
                del imim
                
                ds9.mtv(im, frame = frame)

        pair = algorithms.createKernelFromPsfCandidates(self.cellSet, nEigenComponents, spatialOrder,
                                                        kernelSize, nStarPerCell)

        kernel, eigenValues = pair[0], pair[1]; del pair
        
        psf = afwDetection.createPsf("PCA", kernel) # Hurrah!

        if display:
            xy = []
            showModel = not True
            if showModel:
                oim = self.mi.Factory(self.mi, True)
            for cell in self.cellSet.getCellList():
                for cand in cell:
                    source = algorithms.cast_PsfCandidateF(cand).getSource()

                    delta = -0.5
                    delta =  0.0
                    algorithms.subtractPsf(psf, self.mi,
                                           source.getXAstrom() + delta, source.getYAstrom() + delta)

                    xy.append([source.getXAstrom(), source.getYAstrom()])

            frame = 4
            ds9.mtv(self.mi, frame = frame)
            for xc, yc in xy:
                ds9.dot("x", xc - self.mi.getX0(), yc - self.mi.getY0(), ctype = ds9.GREEN, frame = frame)

            if showModel:
                self.mi -= oim
                self.mi *= -1
                ds9.mtv(self.mi, frame = frame + 1)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(SpatialModelPsfTestCase)
    suites += unittest.makeSuite(psfAttributesTestCase)
    #suites += unittest.makeSuite(RHLTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the utilsTests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
