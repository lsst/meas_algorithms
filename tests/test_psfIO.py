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

import os
import math
import unittest
import tempfile

import numpy as np

import lsst.utils.tests
import lsst.geom
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetection
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
from lsst.log import Log
import lsst.meas.base as measBase
import lsst.meas.algorithms as algorithms

try:
    type(display)
except NameError:
    display = False
else:
    import lsst.afw.display as afwDisplay
    afwDisplay.setDefaultMaskTransparency(75)

# Change the level to Log.DEBUG or Log.TRACE to see debug messages
Log.getLogger("measurement").setLevel(Log.INFO)


def roundTripPsf(key, psf):
    with tempfile.NamedTemporaryFile() as f:
        psf.writeFits(f.name)
        psf2 = type(psf).readFits(f.name)

    return psf2


class SpatialModelPsfTestCase(lsst.utils.tests.TestCase):
    """A test case for SpatialModelPsf"""

    def setUp(self):
        width, height = 100, 300
        self.mi = afwImage.MaskedImageF(lsst.geom.ExtentI(width, height))
        self.mi.set(0)
        self.mi.getVariance().set(10)
        self.mi.getMask().addMaskPlane("DETECTED")

        self.FWHM = 5
        self.ksize = 25                      # size of desired kernel

        self.exposure = afwImage.makeExposure(self.mi)

        psf = roundTripPsf(2, algorithms.DoubleGaussianPsf(self.ksize, self.ksize,
                                                           self.FWHM/(2*math.sqrt(2*math.log(2))), 1, 0.1))
        self.exposure.setPsf(psf)

        for x, y in [(20, 20),
                     # (30, 35), (50, 50),
                     (60, 20), (60, 210), (20, 210)]:

            flux = 10000 - 0*x - 10*y

            sigma = 3 + 0.01*(y - self.mi.getHeight()/2)
            psf = roundTripPsf(3, algorithms.DoubleGaussianPsf(self.ksize, self.ksize, sigma, 1, 0.1))
            im = psf.computeImage().convertF()
            im *= flux
            x0y0 = lsst.geom.PointI(x - self.ksize//2, y - self.ksize//2)
            smi = self.mi.getImage().Factory(self.mi.getImage(),
                                             lsst.geom.BoxI(x0y0, lsst.geom.ExtentI(self.ksize)),
                                             afwImage.LOCAL)

            if False:                   # Test subtraction with non-centered psfs
                im = afwMath.offsetImage(im, 0.5, 0.5)

            smi += im
            del psf
            del im
            del smi

        roundTripPsf(4, algorithms.DoubleGaussianPsf(self.ksize, self.ksize,
                                                     self.FWHM/(2*math.sqrt(2*math.log(2))), 1, 0.1))

        self.cellSet = afwMath.SpatialCellSet(lsst.geom.BoxI(lsst.geom.PointI(0, 0),
                                                             lsst.geom.ExtentI(width, height)), 100)
        ds = afwDetection.FootprintSet(self.mi, afwDetection.Threshold(10), "DETECTED")
        #
        # Prepare to measure
        #
        schema = afwTable.SourceTable.makeMinimalSchema()
        sfm_config = measBase.SingleFrameMeasurementConfig()
        sfm_config.plugins = ["base_SdssCentroid", "base_CircularApertureFlux", "base_PsfFlux",
                              "base_SdssShape", "base_GaussianFlux", "base_PixelFlags"]
        sfm_config.slots.centroid = "base_SdssCentroid"
        sfm_config.slots.shape = "base_SdssShape"
        sfm_config.slots.psfFlux = "base_PsfFlux"
        sfm_config.slots.gaussianFlux = None
        sfm_config.slots.apFlux = "base_CircularApertureFlux_3_0"
        sfm_config.slots.modelFlux = "base_GaussianFlux"
        sfm_config.slots.calibFlux = None
        sfm_config.plugins["base_SdssShape"].maxShift = 10.0
        sfm_config.plugins["base_CircularApertureFlux"].radii = [3.0]
        task = measBase.SingleFrameMeasurementTask(schema, config=sfm_config)
        measCat = afwTable.SourceCatalog(schema)
        # detect the sources and run with the measurement task
        ds.makeSources(measCat)
        task.run(measCat, self.exposure)
        for source in measCat:
            self.cellSet.insertCandidate(algorithms.makePsfCandidate(source, self.exposure))

    def tearDown(self):
        del self.exposure
        del self.cellSet
        del self.mi
        del self.ksize
        del self.FWHM

    def testGetPcaKernel(self):
        """Convert our cellSet to a LinearCombinationKernel"""

        nEigenComponents = 2
        spatialOrder = 1
        kernelSize = 21
        nStarPerCell = 2
        nStarPerCellSpatialFit = 2
        tolerance = 1e-5

        if display:
            disp = afwDisplay.Display(frame=0)
            disp.mtv(self.mi, title=self._testMethodName + ": image")
            #
            # Show the candidates we're using
            #
            for cell in self.cellSet.getCellList():
                i = 0
                for cand in cell:
                    i += 1
                    source = cand.getSource()
                    xc, yc = source.getX() - self.mi.getX0(), source.getY() - self.mi.getY0()
                    if i <= nStarPerCell:
                        disp.dot("o", xc, yc, ctype=afwDisplay.GREEN)
                    else:
                        disp.dot("o", xc, yc, ctype=afwDisplay.YELLOW)

        pair = algorithms.createKernelFromPsfCandidates(self.cellSet, self.exposure.getDimensions(),
                                                        self.exposure.getXY0(), nEigenComponents,
                                                        spatialOrder, kernelSize, nStarPerCell)

        kernel, eigenValues = pair[0], pair[1]
        del pair

        print("lambda", " ".join(["%g" % l for l in eigenValues]))

        pair = algorithms.fitSpatialKernelFromPsfCandidates(kernel, self.cellSet, nStarPerCellSpatialFit,
                                                            tolerance)
        status, chi2 = pair[0], pair[1]
        del pair
        print("Spatial fit: %s chi^2 = %.2g" % (status, chi2))

        psf = roundTripPsf(5, algorithms.PcaPsf(kernel))  # Hurrah!

        self.assertIsNotNone(psf.getKernel())

        self.checkTablePersistence(psf)

        if display:
            # print psf.getKernel().toString()

            eImages = []
            for k in psf.getKernel().getKernelList():
                im = afwImage.ImageD(k.getDimensions())
                k.computeImage(im, False)
                eImages.append(im)

            mos = afwDisplay.utils.Mosaic()
            disp = afwDisplay.Display(frame=3)
            disp.mtv(mos.makeMosaic(eImages), title=self._testMethodName + ": mosaic")
            disp.dot("Eigen Images", 0, 0)
            #
            # Make a mosaic of PSF candidates
            #
            stamps = []
            stampInfo = []

            for cell in self.cellSet.getCellList():
                for cand in cell:
                    s = cand.getSource()
                    im = cand.getMaskedImage()

                    stamps.append(im)
                    stampInfo.append("[%d 0x%x]" % (s.getId(), s["base_PixelFlags_flag"]))

            mos = afwDisplay.utils.Mosaic()
            disp = afwDisplay.Display(frame=1)
            disp.mtv(mos.makeMosaic(stamps), title=self._testMethodName + ": PSF candidates")
            for i in range(len(stampInfo)):
                disp.dot(stampInfo[i], mos.getBBox(i).getMinX(), mos.getBBox(i).getMinY(),
                         ctype=afwDisplay.RED)

            psfImages = []
            labels = []
            if False:
                nx, ny = 3, 4
                for iy in range(ny):
                    for ix in range(nx):
                        x = int((ix + 0.5)*self.mi.getWidth()/nx)
                        y = int((iy + 0.5)*self.mi.getHeight()/ny)

                        im = psf.getImage(x, y)
                        psfImages.append(im.Factory(im, True))
                        labels.append("PSF(%d,%d)" % (int(x), int(y)))

                        if True:
                            print((x, y, "PSF parameters:", psf.getKernel().getKernelParameters()))
            else:
                nx, ny = 2, 2
                for x, y in [(20, 20), (60, 20),
                             (60, 210), (20, 210)]:

                    im = psf.computeImage(lsst.geom.PointD(x, y))
                    psfImages.append(im.Factory(im, True))
                    labels.append("PSF(%d,%d)" % (int(x), int(y)))

                    if True:
                        print(x, y, "PSF parameters:", psf.getKernel().getKernelParameters())
            mos = afwDisplay.utils.Mosaic()
            disp = afwDisplay.Display(frame=2)
            mos.makeMosaic(psfImages, display=disp, mode=nx)
            mos.drawLabels(labels, display=disp)

        if display:
            disp = afwDisplay.Display(frame=0)
            disp.mtv(self.mi, title=self._testMethodName + ": image")

            psfImages = []
            labels = []
            if False:
                nx, ny = 3, 4
                for iy in range(ny):
                    for ix in range(nx):
                        x = int((ix + 0.5)*self.mi.getWidth()/nx)
                        y = int((iy + 0.5)*self.mi.getHeight()/ny)

                        algorithms.subtractPsf(psf, self.mi, x, y)
            else:
                nx, ny = 2, 2
                for x, y in [(20, 20), (60, 20),
                             (60, 210), (20, 210)]:

                    if False:               # Test subtraction with non-centered psfs
                        x += 0.5
                        y -= 0.5

                    # algorithms.subtractPsf(psf, self.mi, x, y)

            afwDisplay.Display(frame=1).mtv(self.mi, title=self._testMethodName + ": image")

    def testCandidateList(self):
        if False and display:
            disp = afwDisplay.Display(frame=0)
            disp.mtv(self.mi, title=self._testMethodName + ": image")

            for cell in self.cellSet.getCellList():
                x0, y0, x1, y1 = (
                    cell.getBBox().getX0(), cell.getBBox().getY0(),
                    cell.getBBox().getX1(), cell.getBBox().getY1())
                print((x0, y0, " ", x1, y1))
                x0 -= 0.5
                y0 -= 0.5
                x1 += 0.5
                y1 += 0.5

                disp.line([(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)], ctype=afwDisplay.RED)

        self.assertFalse(self.cellSet.getCellList()[0].empty())
        self.assertTrue(self.cellSet.getCellList()[1].empty())
        self.assertFalse(self.cellSet.getCellList()[2].empty())

        stamps = []
        for cell in self.cellSet.getCellList():
            for cand in cell:
                cand = cell[0]
                width, height = 15, 17
                cand.setWidth(width)
                cand.setHeight(height)

                im = cand.getMaskedImage()
                stamps.append(im)

                self.assertEqual(im.getWidth(), width)
                self.assertEqual(im.getHeight(), height)

        if display:
            mos = afwDisplay.utils.Mosaic()
            afwDisplay.Display(frame=1).mtv(mos.makeMosaic(stamps), title=self._testMethodName + ": image")

    def checkTablePersistence(self, psf1):
        """Called by testGetPcaKernel to test table-based persistence; it's a pain to
        build a PcaPsf, so we don't want to repeat it all for each test case.

        We just verify here that we get a LinearCombinationKernel; all the details of
        testing that we get the *right* one are tested more thoroughly in afw.
        """
        print("Testing PcaPsf!")
        filename = "PcaPsf.fits"
        psf1.writeFits(filename)
        psf2 = algorithms.PcaPsf.readFits(filename)
        self.assertIsNotNone(psf2)
        self.assertIsNotNone(psf2.getKernel())
        os.remove(filename)


class SingleGaussianPsfTestCase(unittest.TestCase):

    def testTablePersistence(self):
        filename = "SingleGaussianPsf.fits"
        psf1 = algorithms.SingleGaussianPsf(5, 7, 4.2)
        psf1.writeFits(filename)
        psf2 = algorithms.SingleGaussianPsf.readFits(filename)
        self.assertEqual(psf1.getSigma(), psf2.getSigma())
        os.remove(filename)


class DoubleGaussianPsfTestCase(unittest.TestCase):
    """A test case for DoubleGaussianPsf"""

    def comparePsfs(self, psf1, psf2):
        self.assertTrue(isinstance(psf1, algorithms.DoubleGaussianPsf))
        self.assertTrue(isinstance(psf2, algorithms.DoubleGaussianPsf))
        self.assertEqual(psf1.getKernel().getWidth(), psf2.getKernel().getWidth())
        self.assertEqual(psf1.getKernel().getHeight(), psf2.getKernel().getHeight())
        self.assertEqual(psf1.getSigma1(), psf2.getSigma1())
        self.assertEqual(psf1.getSigma2(), psf2.getSigma2())
        self.assertEqual(psf1.getB(), psf2.getB())

    def setUp(self):
        self.ksize = 25                      # size of desired kernel
        FWHM = 5
        self.sigma1 = FWHM/(2*np.sqrt(2*np.log(2)))
        self.sigma2 = 2*self.sigma1
        self.b = 0.1

    def tearDown(self):
        del self.ksize
        del self.sigma1
        del self.sigma2
        del self.b

    def testBoostPersistence(self):
        psf1 = algorithms.DoubleGaussianPsf(self.ksize, self.ksize, self.sigma1, self.sigma2, self.b)
        psf2 = roundTripPsf(1, psf1)
        psf3 = roundTripPsf(1, psf1)
        self.comparePsfs(psf1, psf2)
        self.comparePsfs(psf1, psf3)

    def testFitsPersistence(self):
        psf1 = algorithms.DoubleGaussianPsf(self.ksize, self.ksize, self.sigma1, self.sigma2, self.b)
        filename = "tests/data/psf1-1.fits"
        psf1.writeFits(filename)
        psf2 = algorithms.DoubleGaussianPsf.readFits(filename)
        self.comparePsfs(psf1, psf2)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
