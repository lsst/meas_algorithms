#!/usr/bin/env python
#
# LSST Data Management System
#
# Copyright 2008-2016  AURA/LSST.
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
# see <https://www.lsstcorp.org/LegalNotices/>.
#
from __future__ import absolute_import, division, print_function
from builtins import zip
from builtins import range
import math
import numpy as np
import unittest

from lsst.afw.cameraGeom.testUtils import DetectorWrapper
import lsst.afw.detection as afwDetection
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.afw.display.utils as displayUtils
import lsst.daf.base as dafBase
from lsst.log import Log
import lsst.meas.algorithms as measAlg
from lsst.meas.algorithms.pcaPsfDeterminer import numCandidatesToReject
import lsst.meas.base as measBase
import lsst.utils.tests


try:
    type(display)
    import lsst.afw.display.ds9 as ds9
except NameError:
    display = False

# Change the level to Log.DEBUG or Log.TRACE to see debug messages
Log.getLogger("measurement").setLevel(Log.INFO)
Log.getLogger("psfDeterminer").setLevel(Log.TRACE)


def psfVal(ix, iy, x, y, sigma1, sigma2, b):
    """Return the value at (ix, iy) of a double Gaussian
       (N(0, sigma1^2) + b*N(0, sigma2^2))/(1 + b)
    centered at (x, y)
    """
    return (math.exp(-0.5*((ix - x)**2 + (iy - y)**2)/sigma1**2) +
            b*math.exp(-0.5*((ix - x)**2 + (iy - y)**2)/sigma2**2))/(1 + b)


class SpatialModelPsfTestCase(lsst.utils.tests.TestCase):
    """A test case for SpatialModelPsf"""

    def measure(self, footprintSet, exposure):
        """Measure a set of Footprints, returning a SourceCatalog."""
        table = afwTable.SourceCatalog(self.schema)
        footprintSet.makeSources(table)

        # Then run the default SFM task.  Results not checked
        self.measureTask.run(table, exposure)

        if display:
            ds9.mtv(exposure)

        return table

    def setUp(self):

        self.schema = afwTable.SourceTable.makeMinimalSchema()
        config = measBase.SingleFrameMeasurementConfig()
        config.algorithms.names = ["base_PixelFlags",
                                   "base_SdssCentroid",
                                   "base_GaussianFlux",
                                   "base_SdssShape",
                                   "base_CircularApertureFlux",
                                   "base_PsfFlux",
                                   ]
        config.algorithms["base_CircularApertureFlux"].radii = [3.0]
        config.slots.centroid = "base_SdssCentroid"
        config.slots.psfFlux = "base_PsfFlux"
        config.slots.apFlux = "base_CircularApertureFlux_3_0"
        config.slots.modelFlux = None
        config.slots.instFlux = None
        config.slots.calibFlux = None
        config.slots.shape = "base_SdssShape"

        self.measureTask = measBase.SingleFrameMeasurementTask(self.schema, config=config)

        width, height = 110, 301

        self.mi = afwImage.MaskedImageF(afwGeom.ExtentI(width, height))
        self.mi.set(0)
        sd = 3                          # standard deviation of image
        self.mi.getVariance().set(sd*sd)
        self.mi.getMask().addMaskPlane("DETECTED")

        self.FWHM = 5
        self.ksize = 31                      # size of desired kernel

        sigma1 = 1.75
        sigma2 = 2*sigma1

        self.exposure = afwImage.makeExposure(self.mi)
        self.exposure.setPsf(measAlg.DoubleGaussianPsf(self.ksize, self.ksize,
                                                       1.5*sigma1, 1, 0.1))
        self.exposure.setDetector(DetectorWrapper().detector)

        #
        # Make a kernel with the exactly correct basis functions.  Useful for debugging
        #
        basisKernelList = []
        for sigma in (sigma1, sigma2):
            basisKernel = afwMath.AnalyticKernel(self.ksize, self.ksize,
                                                 afwMath.GaussianFunction2D(sigma, sigma))
            basisImage = afwImage.ImageD(basisKernel.getDimensions())
            basisKernel.computeImage(basisImage, True)
            basisImage /= np.sum(basisImage.getArray())

            if sigma == sigma1:
                basisImage0 = basisImage
            else:
                basisImage -= basisImage0

            basisKernelList.append(afwMath.FixedKernel(basisImage))

        order = 1                                # 1 => up to linear
        spFunc = afwMath.PolynomialFunction2D(order)

        exactKernel = afwMath.LinearCombinationKernel(basisKernelList, spFunc)
        exactKernel.setSpatialParameters([[1.0, 0, 0],
                                          [0.0, 0.5*1e-2, 0.2e-2]])
        self.exactPsf = measAlg.PcaPsf(exactKernel)

        rand = afwMath.Random()               # make these tests repeatable by setting seed

        addNoise = True

        if addNoise:
            im = self.mi.getImage()
            afwMath.randomGaussianImage(im, rand)  # N(0, 1)
            im *= sd                              # N(0, sd^2)
            del im

        xarr, yarr = [], []

        for x, y in [(20, 20), (60, 20),
                     (30, 35),
                     (50, 50),
                     (20, 90), (70, 160), (25, 265), (75, 275), (85, 30),
                     (50, 120), (70, 80),
                     (60, 210), (20, 210),
                     ]:
            xarr.append(x)
            yarr.append(y)

        for x, y in zip(xarr, yarr):
            dx = rand.uniform() - 0.5   # random (centered) offsets
            dy = rand.uniform() - 0.5

            k = exactKernel.getSpatialFunction(1)(x, y)  # functional variation of Kernel ...
            b = (k*sigma1**2/((1 - k)*sigma2**2))       # ... converted double Gaussian's "b"

            # flux = 80000 - 20*x - 10*(y/float(height))**2
            flux = 80000*(1 + 0.1*(rand.uniform() - 0.5))
            I0 = flux*(1 + b)/(2*np.pi*(sigma1**2 + b*sigma2**2))
            for iy in range(y - self.ksize//2, y + self.ksize//2 + 1):
                if iy < 0 or iy >= self.mi.getHeight():
                    continue

                for ix in range(x - self.ksize//2, x + self.ksize//2 + 1):
                    if ix < 0 or ix >= self.mi.getWidth():
                        continue

                    I = I0*psfVal(ix, iy, x + dx, y + dy, sigma1, sigma2, b)
                    Isample = rand.poisson(I) if addNoise else I
                    self.mi.getImage().set(ix, iy, self.mi.getImage().get(ix, iy) + Isample)
                    self.mi.getVariance().set(ix, iy, self.mi.getVariance().get(ix, iy) + I)
        #
        bbox = afwGeom.BoxI(afwGeom.PointI(0, 0), afwGeom.ExtentI(width, height))
        self.cellSet = afwMath.SpatialCellSet(bbox, 100)

        self.footprintSet = afwDetection.FootprintSet(self.mi, afwDetection.Threshold(100), "DETECTED")
        self.catalog = self.measure(self.footprintSet, self.exposure)

        for source in self.catalog:
            try:
                cand = measAlg.makePsfCandidate(source, self.exposure)
                self.cellSet.insertCandidate(cand)

            except Exception as e:
                print(e)
                continue

    def tearDown(self):
        del self.cellSet
        del self.exposure
        del self.mi
        del self.exactPsf
        del self.footprintSet
        del self.catalog
        del self.schema
        del self.measureTask

    def setupDeterminer(self, exposure=None, nEigenComponents=2, starSelectorAlg="secondMoment"):
        """Setup the starSelector and psfDeterminer."""
        if exposure is None:
            exposure = self.exposure

        starSelectorClass = measAlg.starSelectorRegistry[starSelectorAlg]
        starSelectorConfig = starSelectorClass.ConfigClass()

        if starSelectorAlg == "secondMoment":
            starSelectorConfig.clumpNSigma = 5.0
            starSelectorConfig.histSize = 14
            starSelectorConfig.badFlags = ["base_PixelFlags_flag_edge",
                                           "base_PixelFlags_flag_interpolatedCenter",
                                           "base_PixelFlags_flag_saturatedCenter",
                                           "base_PixelFlags_flag_crCenter",
                                           ]
        elif starSelectorAlg == "objectSize":
            starSelectorConfig.sourceFluxField = "base_GaussianFlux_flux"
            starSelectorConfig.badFlags = ["base_PixelFlags_flag_edge",
                                           "base_PixelFlags_flag_interpolatedCenter",
                                           "base_PixelFlags_flag_saturatedCenter",
                                           "base_PixelFlags_flag_crCenter",
                                           ]
            starSelectorConfig.widthStdAllowed = 0.5

        starSelector = starSelectorClass(config=starSelectorConfig, schema=self.schema)

        psfDeterminerTask = measAlg.psfDeterminerRegistry["pca"]
        psfDeterminerConfig = psfDeterminerTask.ConfigClass()
        width, height = exposure.getMaskedImage().getDimensions()
        psfDeterminerConfig.sizeCellX = width
        psfDeterminerConfig.sizeCellY = height//3
        psfDeterminerConfig.nEigenComponents = nEigenComponents
        psfDeterminerConfig.spatialOrder = 1
        psfDeterminerConfig.kernelSizeMin = 31
        psfDeterminerConfig.nStarPerCell = 0
        psfDeterminerConfig.nStarPerCellSpatialFit = 0  # unlimited
        psfDeterminer = psfDeterminerTask(psfDeterminerConfig)

        return starSelector, psfDeterminer

    def subtractStars(self, exposure, catalog, chi_lim=-1):
        """Subtract the exposure's PSF from all the sources in catalog."""
        mi, psf = exposure.getMaskedImage(), exposure.getPsf()

        subtracted = mi.Factory(mi, True)

        for s in catalog:
            xc, yc = s.getX(), s.getY()
            bbox = subtracted.getBBox()
            if bbox.contains(afwGeom.PointI(int(xc), int(yc))):
                try:
                    measAlg.subtractPsf(psf, subtracted, xc, yc)
                except:
                    pass

        chi = subtracted.Factory(subtracted, True)
        var = subtracted.getVariance()
        np.sqrt(var.getArray(), var.getArray())  # inplace sqrt
        chi /= var

        if display:
            ds9.mtv(subtracted, title="Subtracted", frame=1)
            ds9.mtv(chi, title="Chi", frame=2)
            ds9.mtv(psf.computeImage(afwGeom.Point2D(xc, yc)), title="Psf", frame=3)
            ds9.mtv(mi, frame=4, title="orig")
            kern = psf.getKernel()
            kimg = afwImage.ImageD(kern.getWidth(), kern.getHeight())
            kern.computeImage(kimg, True, xc, yc)
            ds9.mtv(kimg, title="kernel", frame=5)

        chi_min, chi_max = np.min(chi.getImage().getArray()), np.max(chi.getImage().getArray())
        if False:
            print(chi_min, chi_max)

        if chi_lim > 0:
            self.assertGreater(chi_min, -chi_lim)
            self.assertLess(chi_max, chi_lim)

    def testPsfDeterminer(self):
        """Test the (PCA) psfDeterminer."""
        for starSelectorAlg in ["secondMoment",
                                "objectSize",
                                ]:
            print("Using %s star selector" % (starSelectorAlg))

            starSelector, psfDeterminer = self.setupDeterminer(starSelectorAlg=starSelectorAlg)
            metadata = dafBase.PropertyList()
            psfCandidateList = starSelector.run(self.exposure, self.catalog).psfCandidates
            psf, cellSet = psfDeterminer.determinePsf(self.exposure, psfCandidateList, metadata)
            self.exposure.setPsf(psf)

            chi_lim = 5.0
            self.subtractStars(self.exposure, self.catalog, chi_lim)

    def testPsfDeterminerSubimageObjectSizeStarSelector(self):
        """Test the (PCA) psfDeterminer on subImages."""
        w, h = self.exposure.getDimensions()
        x0, y0 = int(0.35*w), int(0.45*h)
        bbox = afwGeom.BoxI(afwGeom.PointI(x0, y0), afwGeom.ExtentI(w-x0, h-y0))
        subExp = self.exposure.Factory(self.exposure, bbox, afwImage.LOCAL)

        starSelector, psfDeterminer = self.setupDeterminer(subExp, starSelectorAlg="objectSize")
        metadata = dafBase.PropertyList()
        #
        # Only keep the sources that lie within the subregion (avoiding lots of log messages)
        #

        def trimCatalogToImage(exp, catalog):
            trimmedCatalog = afwTable.SourceCatalog(catalog.table.clone())
            for s in catalog:
                if exp.getBBox().contains(afwGeom.PointI(s.getCentroid())):
                    trimmedCatalog.append(trimmedCatalog.table.copyRecord(s))

            return trimmedCatalog

        psfCandidateList = starSelector.run(subExp, trimCatalogToImage(subExp, self.catalog)).psfCandidates
        psf, cellSet = psfDeterminer.determinePsf(subExp, psfCandidateList, metadata)
        subExp.setPsf(psf)

        # Test how well we can subtract the PSF model.  N.b. using self.exposure is an extrapolation
        for exp, chi_lim in [(subExp, 4.5),
                             (self.exposure.Factory(self.exposure,
                                                    afwGeom.BoxI(afwGeom.PointI(0, 100),
                                                                 (afwGeom.PointI(w-1, h-1))),
                                                    afwImage.LOCAL), 7.5),
                             (self.exposure, 19),
                             ]:
            cat = trimCatalogToImage(exp, self.catalog)
            exp.setPsf(psf)
            self.subtractStars(exp, cat, chi_lim)

    def testPsfDeterminerSubimageSecondMomentStarSelector(self):
        """Test the (PCA) psfDeterminer on subImages."""
        w, h = self.exposure.getDimensions()
        x0, y0 = int(0.35*w), int(0.45*h)
        bbox = afwGeom.BoxI(afwGeom.PointI(x0, y0), afwGeom.ExtentI(w-x0, h-y0))
        subExp = self.exposure.Factory(self.exposure, bbox, afwImage.LOCAL)

        starSelector, psfDeterminer = self.setupDeterminer(subExp, starSelectorAlg="secondMoment")
        metadata = dafBase.PropertyList()
        #
        # Only keep the sources that lie within the subregion (avoiding lots of log messages)
        #

        def trimCatalogToImage(exp, catalog):
            trimmedCatalog = afwTable.SourceCatalog(catalog.table.clone())
            for s in catalog:
                if exp.getBBox().contains(afwGeom.PointI(s.getCentroid())):
                    trimmedCatalog.append(trimmedCatalog.table.copyRecord(s))

            return trimmedCatalog

        psfCandidateList = starSelector.run(subExp, trimCatalogToImage(subExp, self.catalog)).psfCandidates
        psf, cellSet = psfDeterminer.determinePsf(subExp, psfCandidateList, metadata)
        subExp.setPsf(psf)

        # Test how well we can subtract the PSF model.  N.b. using self.exposure is an extrapolation
        for exp, chi_lim in [(subExp, 4.5),
                             (self.exposure.Factory(self.exposure,
                                                    afwGeom.BoxI(afwGeom.PointI(0, 100),
                                                                 (afwGeom.PointI(w-1, h-1))),
                                                    afwImage.LOCAL), 7.5),
                             (self.exposure, 19),
                             ]:
            cat = trimCatalogToImage(exp, self.catalog)
            exp.setPsf(psf)
            self.subtractStars(exp, cat, chi_lim)

    def testPsfDeterminerNEigenObjectSizeStarSelector(self):
        """Test the (PCA) psfDeterminer when you ask for more components than acceptable stars."""
        starSelector, psfDeterminer = self.setupDeterminer(nEigenComponents=3, starSelectorAlg="objectSize")
        metadata = dafBase.PropertyList()
        psfCandidateList = starSelector.run(self.exposure, self.catalog).psfCandidates
        psfCandidateList, nEigen = psfCandidateList[0:4], 2  # only enough stars for 2 eigen-components
        psf, cellSet = psfDeterminer.determinePsf(self.exposure, psfCandidateList, metadata)

        self.assertEqual(psf.getKernel().getNKernelParameters(), nEigen)

    def testPsfDeterminerNEigenSecondMomentStarSelector(self):
        """Test the (PCA) psfDeterminer when you ask for more components than acceptable stars."""
        starSelector, psfDeterminer = self.setupDeterminer(nEigenComponents=3, starSelectorAlg="secondMoment")
        metadata = dafBase.PropertyList()
        psfCandidateList = starSelector.run(self.exposure, self.catalog).psfCandidates
        psfCandidateList, nEigen = psfCandidateList[0:4], 2  # only enough stars for 2 eigen-components
        psf, cellSet = psfDeterminer.determinePsf(self.exposure, psfCandidateList, metadata)

        self.assertEqual(psf.getKernel().getNKernelParameters(), nEigen)

    def testCandidateList(self):
        self.assertFalse(self.cellSet.getCellList()[0].empty())
        self.assertTrue(self.cellSet.getCellList()[1].empty())
        self.assertFalse(self.cellSet.getCellList()[2].empty())
        self.assertTrue(self.cellSet.getCellList()[3].empty())

        stamps = []
        for cell in self.cellSet.getCellList():
            for cand in cell:
                cand = cell[0]
                width, height = 29, 25
                cand.setWidth(width)
                cand.setHeight(height)

                im = cand.getMaskedImage()
                stamps.append(im)

                self.assertEqual(im.getWidth(), width)
                self.assertEqual(im.getHeight(), height)

        if False and display:
            mos = displayUtils.Mosaic()
            mos.makeMosaic(stamps, frame=2)

    def testRejectBlends(self):
        """Test the PcaPsfDeterminerTask blend removal."""
        """
        We give it a single blended source, asking it to remove blends,
        and check that it barfs in the expected way.
        """

        psfDeterminerClass = measAlg.psfDeterminerRegistry["pca"]
        config = psfDeterminerClass.ConfigClass()
        config.doRejectBlends = True
        psfDeterminer = psfDeterminerClass(config=config)

        schema = afwTable.SourceTable.makeMinimalSchema()
        # Use The single frame measurement task to populate the schema with standard keys
        measBase.SingleFrameMeasurementTask(schema)
        catalog = afwTable.SourceCatalog(schema)
        source = catalog.addNew()

        # Make the source blended, with necessary information to calculate pca
        spanShift = afwGeom.Point2I(54, 123)
        spans = afwGeom.SpanSet.fromShape(6, offset=spanShift)
        foot = afwDetection.Footprint(spans, self.exposure.getBBox())
        foot.addPeak(45, 123, 6)
        foot.addPeak(47, 126, 5)
        source.setFootprint(foot)
        centerKey = afwTable.Point2DKey(source.schema['slot_Centroid'])
        shapeKey = afwTable.QuadrupoleKey(schema['slot_Shape'])
        source.set(centerKey, afwGeom.Point2D(46, 124))
        source.set(shapeKey, afwGeom.ellipses.Quadrupole(1.1, 2.2, 1))

        candidates = [measAlg.makePsfCandidate(source, self.exposure)]
        metadata = dafBase.PropertyList()

        with self.assertRaises(RuntimeError) as cm:
            psfDeterminer.determinePsf(self.exposure, candidates, metadata)
        self.assertEqual(str(cm.exception), "All PSF candidates removed as blends")


class PsfCandidateTestCase(lsst.utils.tests.TestCase):
    def testNumToReject(self):
        """Reject the correct number of PSF candidates on each iteration"""
        # Numerical values correspond to the problem case identified in
        # DM-8030.

        numBadCandidates = 5
        totalIter = 3

        for numIter, value in [(0, 1), (1, 3), (2, 5)]:
            self.assertEqual(numCandidatesToReject(numBadCandidates, numIter,
                                                   totalIter), value)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
