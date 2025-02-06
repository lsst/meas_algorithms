# This file is part of meas_algorithms.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import math
import unittest
import time

import numpy as np

import lsst.daf.base as dafBase
import lsst.geom
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.meas.algorithms as measAlg
import lsst.meas.base as measBase

import lsst.afw.cameraGeom as cameraGeom
from lsst.afw.cameraGeom.testUtils import DetectorWrapper

import lsst.utils.tests

try:
    display
except NameError:
    display = False
else:
    import lsst.afw.display as afwDisplay
    afwDisplay.setDefaultMaskTransparency(75)


def plantSources(x0, y0, nx, ny, sky, nObj, wid, detector, rng, useRandom=False):

    pixToTanPix = detector.getTransform(cameraGeom.PIXELS, cameraGeom.TAN_PIXELS)

    img0 = afwImage.ImageF(lsst.geom.ExtentI(nx, ny))
    img = afwImage.ImageF(lsst.geom.ExtentI(nx, ny))

    ixx0, iyy0, ixy0 = wid*wid, wid*wid, 0.0

    edgeBuffer = 40.0*wid

    flux = 1.0e4
    nkx, nky = int(10*wid) + 1, int(10*wid) + 1
    xhwid, yhwid = nkx//2, nky//2

    nRow = int(math.sqrt(nObj))
    xstep = (nx - 1 - 0.0*edgeBuffer)//(nRow+1)
    ystep = (ny - 1 - 0.0*edgeBuffer)//(nRow+1)

    if useRandom:
        nObj = nRow*nRow

    goodAdded0 = []
    goodAdded = []

    for i in range(nObj):

        # get our position
        if useRandom:
            xcen0, ycen0 = rng.uniform(nx), rng.uniform(ny)
        else:
            xcen0, ycen0 = xstep*((i % nRow) + 1), ystep*(int(i/nRow) + 1)
        ixcen0, iycen0 = int(xcen0), int(ycen0)

        # distort position and shape
        pTan = lsst.geom.Point2D(xcen0, ycen0)
        p = pixToTanPix.applyInverse(pTan)
        linTransform = afwGeom.linearizeTransform(pixToTanPix, p).inverted().getLinear()
        m = afwGeom.Quadrupole(ixx0, iyy0, ixy0)
        m.transform(linTransform)

        xcen, ycen = xcen0, ycen0  # p.getX(), p.getY()
        if (xcen < 1.0*edgeBuffer or (nx - xcen) < 1.0*edgeBuffer
                or ycen < 1.0*edgeBuffer or (ny - ycen) < 1.0*edgeBuffer):
            continue
        ixcen, iycen = int(xcen), int(ycen)
        ixx, iyy, ixy = m.getIxx(), m.getIyy(), m.getIxy()

        # plant the object
        tmp = 0.25*(ixx-iyy)**2 + ixy**2
        a2 = 0.5*(ixx+iyy) + np.sqrt(tmp)
        b2 = 0.5*(ixx+iyy) - np.sqrt(tmp)

        theta = 0.5*np.arctan2(2.0*ixy, ixx-iyy)
        a = np.sqrt(a2)
        b = np.sqrt(b2)

        c, s = math.cos(theta), math.sin(theta)
        good0, good = True, True
        for y in range(nky):
            iy = iycen + y - yhwid
            iy0 = iycen0 + y - yhwid

            for x in range(nkx):
                ix = ixcen + x - xhwid
                ix0 = ixcen0 + x - xhwid

                if ix >= 0 and ix < nx and iy >= 0 and iy < ny:
                    dx, dy = ix - xcen, iy - ycen
                    u = c*dx + s*dy
                    v = -s*dx + c*dy
                    I0 = flux/(2*math.pi*a*b)
                    val = I0*math.exp(-0.5*((u/a)**2 + (v/b)**2))
                    if val < 0:
                        val = 0
                    prevVal = img[ix, iy, afwImage.LOCAL]
                    img[ix, iy, afwImage.LOCAL] = val+prevVal
                else:
                    good = False

                if ix0 >= 0 and ix0 < nx and iy0 >= 0 and iy0 < ny:
                    dx, dy = ix - xcen, iy - ycen
                    I0 = flux/(2*math.pi*wid*wid)
                    val = I0*math.exp(-0.5*((dx/wid)**2 + (dy/wid)**2))
                    if val < 0:
                        val = 0
                    prevVal = img0[ix0, iy0, afwImage.LOCAL]
                    img0[ix0, iy0, afwImage.LOCAL] = val+prevVal
                else:
                    good0 = False

        if good0:
            goodAdded0.append([xcen, ycen])
        if good:
            goodAdded.append([xcen, ycen])

    # add sky and noise
    img += sky
    img0 += sky
    noise = afwImage.ImageF(lsst.geom.ExtentI(nx, ny))
    noise0 = afwImage.ImageF(lsst.geom.ExtentI(nx, ny))
    for i in range(nx):
        for j in range(ny):
            noise[i, j, afwImage.LOCAL] = rng.poisson(img[i, j, afwImage.LOCAL])
            noise0[i, j, afwImage.LOCAL] = rng.poisson(img0[i, j, afwImage.LOCAL])

    edgeWidth = int(0.5*edgeBuffer)
    mask = afwImage.Mask(lsst.geom.ExtentI(nx, ny))
    left = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.ExtentI(edgeWidth, ny))
    right = lsst.geom.Box2I(lsst.geom.Point2I(nx - edgeWidth, 0), lsst.geom.ExtentI(edgeWidth, ny))
    top = lsst.geom.Box2I(lsst.geom.Point2I(0, ny - edgeWidth), lsst.geom.ExtentI(nx, edgeWidth))
    bottom = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.ExtentI(nx, edgeWidth))

    for pos in [left, right, top, bottom]:
        msk = afwImage.Mask(mask, pos, deep=False)
        msk.set(msk.getPlaneBitMask('EDGE'))

    expos = afwImage.makeExposure(afwImage.makeMaskedImage(noise, mask, afwImage.ImageF(noise, True)))
    expos0 = afwImage.makeExposure(afwImage.makeMaskedImage(noise0, mask, afwImage.ImageF(noise0, True)))

    im = expos.getMaskedImage().getImage()
    im0 = expos0.getMaskedImage().getImage()
    im -= sky
    im0 -= sky

    return expos, goodAdded, expos0, goodAdded0


class PsfSelectionTestCase(lsst.utils.tests.TestCase):
    """Test the aperture correction."""

    def setUp(self):
        self.rng = np.random.Generator(np.random.MT19937(500))
        self.x0, self.y0 = 0, 0
        self.nx, self.ny = 512, 512  # 2048, 4096
        self.sky = 100.0
        self.nObj = 100

        # make a detector with distortion
        self.detector = DetectorWrapper(
            bbox=lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(self.nx, self.ny)),
            orientation=cameraGeom.Orientation(lsst.geom.Point2D(255.0, 255.0)),
            radialDistortion=0.925,
        ).detector

        # make a detector with no distortion
        self.flatDetector = DetectorWrapper(
            bbox=lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(self.nx, self.ny)),
            orientation=cameraGeom.Orientation(lsst.geom.Point2D(255.0, 255.0)),
            radialDistortion=0.0,
        ).detector

        # detection policies
        detConfig = measAlg.SourceDetectionConfig()
        # Cannot use default background approximation order (6) for such a small image.
        detConfig.background.approxOrderX = 4
        # This test depends on footprints grown with the old Manhattan metric.
        detConfig.isotropicGrow = False

        # measurement policies
        measConfig = measBase.SingleFrameMeasurementConfig()
        measConfig.algorithms.names = [
            "base_SdssCentroid",
            "base_SdssShape",
            "base_GaussianFlux",
            "base_PsfFlux",
        ]
        measConfig.slots.centroid = "base_SdssCentroid"
        measConfig.slots.shape = "base_SdssShape"
        measConfig.slots.psfFlux = "base_PsfFlux"
        measConfig.plugins["base_SdssCentroid"].doFootprintCheck = False
        measConfig.slots.apFlux = None
        measConfig.slots.modelFlux = None
        measConfig.slots.gaussianFlux = None
        measConfig.slots.calibFlux = None

        self.schema = afwTable.SourceTable.makeMinimalSchema()
        self.detTask = measAlg.SourceDetectionTask(config=detConfig, schema=self.schema)
        self.measTask = measBase.SingleFrameMeasurementTask(config=measConfig, schema=self.schema)

        # psf star selector
        starSelectorClass = measAlg.sourceSelectorRegistry["objectSize"]
        starSelectorConfig = starSelectorClass.ConfigClass()
        starSelectorConfig.fluxMin = 5000.0
        starSelectorConfig.badFlags = []
        self.starSelector = starSelectorClass(config=starSelectorConfig)

        self.makePsfCandidates = measAlg.MakePsfCandidatesTask()

        # psf determiner
        psfDeterminerFactory = measAlg.psfDeterminerRegistry["pca"]
        psfDeterminerConfig = psfDeterminerFactory.ConfigClass()
        width, height = self.nx, self.ny
        nEigenComponents = 3
        psfDeterminerConfig.sizeCellX = width//3
        psfDeterminerConfig.sizeCellY = height//3
        psfDeterminerConfig.nEigenComponents = nEigenComponents
        psfDeterminerConfig.spatialOrder = 1
        psfDeterminerConfig.nStarPerCell = 0
        psfDeterminerConfig.nStarPerCellSpatialFit = 0  # unlimited
        self.psfDeterminer = psfDeterminerFactory(psfDeterminerConfig)

    def tearDown(self):
        del self.detTask
        del self.measTask
        del self.schema
        del self.detector
        del self.flatDetector
        del self.starSelector
        del self.psfDeterminer

    def detectAndMeasure(self, exposure):
        """Quick and dirty detection (note: we already subtracted background)
        """
        table = afwTable.SourceTable.make(self.schema)
        # detect
        sources = self.detTask.run(table, exposure).sources
        # ... and measure
        self.measTask.run(sources, exposure)
        return sources

    def testPsfCandidate(self):

        detector = self.detector

        # make an exposure
        print("Planting")
        psfSigma = 1.5
        exposDist, nGoodDist, expos0, nGood0 = plantSources(self.x0, self.y0,
                                                            self.nx, self.ny,
                                                            self.sky, self.nObj, psfSigma, detector, self.rng)

        # set the psf
        kwid = 21
        psf = measAlg.SingleGaussianPsf(kwid, kwid, psfSigma)
        exposDist.setPsf(psf)
        exposDist.setDetector(detector)

        # detect
        print("detection")
        sourceList = self.detectAndMeasure(exposDist)

        # select psf stars
        print("PSF selection")
        stars = self.starSelector.run(sourceList, exposure=exposDist)
        psfCandidateList = self.makePsfCandidates.run(stars.sourceCat, exposDist).psfCandidates

        # determine the PSF
        print("PSF determination")
        metadata = dafBase.PropertyList()
        t0 = time.time()
        psf, cellSet = self.psfDeterminer.determinePsf(exposDist, psfCandidateList, metadata)
        print("... determination time: ", time.time() - t0)
        print("PSF kernel width: ", psf.getKernel().getWidth())

        #######################################################################
        # try to subtract off the stars and check the residuals

        imgOrig = exposDist.getMaskedImage().getImage().getArray()
        maxFlux = imgOrig.max()

        ############
        # first try it with no distortion in the psf
        exposDist.setDetector(self.flatDetector)

        print("uncorrected subtraction")
        subImg = afwImage.MaskedImageF(exposDist.getMaskedImage(), True)
        for s in sourceList:
            x, y = s.getX(), s.getY()
            measAlg.subtractPsf(psf, subImg, x, y)

        if display:
            afwDisplay.Display(frame=1).mtv(exposDist, title=self._testMethodName + ": full")
            afwDisplay.Display(frame=0).mtv(subImg, title=self._testMethodName + ": subtracted")

        img = subImg.getImage().getArray()
        norm = img/math.sqrt(maxFlux)

        smin0, smax0, srms0 = norm.min(), norm.max(), norm.std()

        print("min:", smin0, "max: ", smax0, "rms: ", srms0)

        if False:
            # This section has been disabled as distortion was removed from PsfCandidate and Psf;
            # it will be reintroduced in the future with a different API, at which point this
            # test code should be re-enabled.

            ##############
            # try it with the correct distortion in the psf
            exposDist.setDetector(self.detector)

            print("corrected subtraction")
            subImg = afwImage.MaskedImageF(exposDist.getMaskedImage(), True)
            for s in sourceList:
                x, y = s.getX(), s.getY()
                measAlg.subtractPsf(psf, subImg, x, y)

            if display:
                afwDisplay.Display(frame=2).mtv(exposDist, title=self._testMethodName + ": full")
                afwDisplay.Display(frame=3).mtv(subImg, title=self._testMethodName + ": subtracted")

            img = subImg.getImage().getArray()
            norm = img/math.sqrt(maxFlux)

            smin, smax, srms = norm.min(), norm.max(), norm.std()

            # with proper distortion, residuals should be < 4sigma (even for 512x512 pixels)
            print("min:", smin, "max: ", smax, "rms: ", srms)

            # the distrib of residuals should be tighter
            self.assertLess(smin0, smin)
            self.assertGreater(smax0, smax)
            self.assertGreater(srms0, srms)

    def testDistortedImage(self):

        detector = self.detector

        psfSigma = 1.5
        stars = plantSources(
            self.x0, self.y0, self.nx, self.ny, self.sky, self.nObj, psfSigma, detector, self.rng
        )
        expos, starXy = stars[0], stars[1]

        # add some faint round galaxies ... only slightly bigger than the psf
        gxy = plantSources(
            self.x0, self.y0, self.nx, self.ny, self.sky, 10, 1.07*psfSigma, detector, self.rng
        )
        mi = expos.getMaskedImage()
        mi += gxy[0].getMaskedImage()
        gxyXy = gxy[1]

        kwid = 15  # int(10*psfSigma) + 1
        psf = measAlg.SingleGaussianPsf(kwid, kwid, psfSigma)
        expos.setPsf(psf)

        expos.setDetector(detector)

        ########################
        # try without distorter
        expos.setDetector(self.flatDetector)
        print("Testing PSF selection *without* distortion")
        sourceList = self.detectAndMeasure(expos)
        stars = self.starSelector.run(sourceList, exposure=expos)
        psfCandidateList = self.makePsfCandidates.run(stars.sourceCat, expos).psfCandidates

        ########################
        # try with distorter
        expos.setDetector(self.detector)
        print("Testing PSF selection *with* distortion")
        sourceList = self.detectAndMeasure(expos)
        stars = self.starSelector.run(sourceList, exposure=expos)
        psfCandidateListCorrected = self.makePsfCandidates.run(stars.sourceCat, expos).psfCandidates

        def countObjects(candList):
            nStar, nGxy = 0, 0
            for c in candList:
                s = c.getSource()
                x, y = s.getX(), s.getY()
                for xs, ys in starXy:
                    if abs(x-xs) < 2.0 and abs(y-ys) < 2.0:
                        nStar += 1
                for xg, yg in gxyXy:
                    if abs(x-xg) < 2.0 and abs(y-yg) < 2.0:
                        nGxy += 1
            return nStar, nGxy

        nstar, ngxy = countObjects(psfCandidateList)
        nstarC, ngxyC = countObjects(psfCandidateListCorrected)

        print("uncorrected nStar, nGxy: ", nstar, "/", len(starXy), "   ", ngxy, '/', len(gxyXy))
        print("dist-corrected nStar, nGxy: ", nstarC, '/', len(starXy), "   ", ngxyC, '/', len(gxyXy))

        ########################
        # display
        if display:
            iDisp = 1
            disp = afwDisplay.Display(frame=iDisp)
            disp.mtv(expos, title=self._testMethodName + ": image")
            size = 40
            for c in psfCandidateList:
                s = c.getSource()
                ixx, iyy, ixy = size*s.getIxx(), size*s.getIyy(), size*s.getIxy()
                disp.dot("@:%g,%g,%g" % (ixx, ixy, iyy), s.getX(), s.getY(), ctype=afwDisplay.RED)
            size *= 2.0
            for c in psfCandidateListCorrected:
                s = c.getSource()
                ixx, iyy, ixy = size*s.getIxx(), size*s.getIyy(), size*s.getIxy()
                disp.dot("@:%g,%g,%g" % (ixx, ixy, iyy), s.getX(), s.getY(), ctype=afwDisplay.GREEN)

        # we shouldn't expect to get all available stars without distortion correcting
        self.assertLess(nstar, len(starXy))

        # here we should get all of them, occassionally 1 or 2 might get missed
        self.assertGreaterEqual(nstarC, 0.95*len(starXy))

        # no contamination by small gxys
        self.assertEqual(ngxyC, 0)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
