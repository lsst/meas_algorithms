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
import unittest
import numpy as np

import lsst.geom
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
from lsst.meas.algorithms import SourceDetectionTask, SubtractBackgroundTask
from lsst.meas.algorithms.testUtils import plantSources
import lsst.utils.tests

# To plot in ds9, `setup display_ds9` first, and open a ds9 window.
# import lsstDebug
# def DebugInfo(name):
#     debug = lsstDebug.getInfo(name)
#     if name == "lsst.meas.algorithms.detection":
#         debug.display = 2
#     return debug
# lsstDebug.Info = DebugInfo


class SourceDetectionTaskTestCase(lsst.utils.tests.TestCase):

    def _create_exposure(self, bigBbox=False):
        """Return a simulated exposure (and relevant parameters) with Gaussian
        stars.
        """
        if bigBbox:
            bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(1024, 1024))
        else:
            bbox = lsst.geom.Box2I(lsst.geom.Point2I(256, 100), lsst.geom.Extent2I(128, 127))
        minCounts = 5000
        maxCounts = 50000
        starSigma = 1.5
        numX = 5
        numY = 5
        coordList = self.makeCoordList(
            bbox=bbox,
            numX=numX,
            numY=numY,
            minCounts=minCounts,
            maxCounts=maxCounts,
            sigma=starSigma,
        )
        kwid = 11
        sky = 2000
        addPoissonNoise = True
        exposure = plantSources(bbox=bbox, kwid=kwid, sky=sky, coordList=coordList,
                                addPoissonNoise=addPoissonNoise)
        schema = afwTable.SourceTable.makeMinimalSchema()

        return exposure, schema, numX, numY, starSigma

    def _check_detectFootprints(self, exposure, numX, numY, starSigma, task, config, doSmooth=False):
        """Run detectFootprints and check that the output is reasonable,
        for either value of doSmooth.
        """
        taskSigma = 2.2
        res = task.detectFootprints(exposure, doSmooth=doSmooth, sigma=taskSigma)
        self.assertEqual(res.numPos, numX * numY)
        self.assertEqual(res.numNeg, 0)
        self.assertEqual(task.metadata.getScalar("sigma"), taskSigma)
        self.assertEqual(task.metadata.getScalar("doSmooth"), doSmooth)
        self.assertEqual(task.metadata.getScalar("nGrow"), int(taskSigma * config.nSigmaToGrow + 0.5))

        res = task.detectFootprints(exposure, doSmooth=doSmooth, sigma=None)
        taskSigma = task.metadata.getScalar("sigma")
        self.assertLess(abs(taskSigma - starSigma), 0.1)
        self.assertEqual(res.numPos, numX * numY)
        self.assertEqual(res.numNeg, 0)
        return res

    def test_stdev(self):
        """Test that sources are detected on a simulated image with
        thresholdType='stdev'.
        """
        exposure, schema, numX, numY, starSigma = self._create_exposure()

        config = SourceDetectionTask.ConfigClass()
        # don't modify the image after detection.
        config.reEstimateBackground = False
        config.thresholdType = "stdev"
        task = SourceDetectionTask(config=config, schema=schema)

        self._check_detectFootprints(exposure, numX, numY, starSigma, task, config, doSmooth=True)
        self._check_detectFootprints(exposure, numX, numY, starSigma, task, config, doSmooth=False)

    def test_significance_stdev(self):
        """Check the non-smoothed, non-background updated peak significance
        values with thresholdType="stddev".
        """
        exposure, schema, numX, numY, starSigma = self._create_exposure()

        config = SourceDetectionTask.ConfigClass()
        # don't modify the image after detection.
        config.reEstimateBackground = False
        config.doTempLocalBackground = False
        config.thresholdType = "stdev"
        task = SourceDetectionTask(config=config, schema=schema)

        result = self._check_detectFootprints(exposure, numX, numY, starSigma, task, config, doSmooth=False)

        bad = exposure.mask.getPlaneBitMask(config.statsMask)
        sctrl = afwMath.StatisticsControl()
        sctrl.setAndMask(bad)
        stats = afwMath.makeStatistics(exposure.maskedImage, afwMath.STDEVCLIP, sctrl)
        stddev = stats.getValue(afwMath.STDEVCLIP)
        for footprint in result.positive.getFootprints():
            for peak in footprint.peaks:
                point = lsst.geom.Point2I(peak.getIx(), peak.getIy())
                value = exposure.image[point]
                with self.subTest(str(point)):
                    self.assertFloatsAlmostEqual(peak["significance"],
                                                 value/stddev,
                                                 rtol=1e-7,
                                                 msg=str(point))

    def test_pixel_stdev(self):
        """Test that sources are detected on a simulated image with
        thresholdType='pixel_stdev', and that they have the right significance.
        """
        exposure, schema, numX, numY, starSigma = self._create_exposure()

        config = SourceDetectionTask.ConfigClass()
        config.thresholdType = "pixel_stdev"
        config.reEstimateBackground = False
        # TempLocalBackground changes the peak value of the faintest peak,
        # so disable it for this test so that we can calculate an expected
        # answer without having to try to deal with backgrounds.
        config.doTempLocalBackground = False
        task = SourceDetectionTask(config=config, schema=schema)
        # Don't smooth, so that we can directly calculate the s/n from the exposure.
        result = task.detectFootprints(exposure, doSmooth=False)
        self.assertEqual(result.numPos, numX * numY)
        self.assertEqual(result.numNeg, 0)
        # Significance values for `pixel_stdev` should match image/sqrt(variance).
        for footprint in result.positive.getFootprints():
            for peak in footprint.peaks:
                point = lsst.geom.Point2I(peak.getIx(), peak.getIy())
                value = exposure.image[point]
                stddev = np.sqrt(exposure.variance[point])
                with self.subTest(str(point)):
                    self.assertFloatsAlmostEqual(peak["significance"],
                                                 value/stddev,
                                                 rtol=1e-7,
                                                 msg=str(point))

    def makeCoordList(self, bbox, numX, numY, minCounts, maxCounts, sigma):
        """Make a coordList for plantSources."""
        """
        Coords are uniformly spaced in a rectangular grid, with linearly increasing counts
        """
        dX = bbox.getWidth() / float(numX)
        dY = bbox.getHeight() / float(numY)
        minX = bbox.getMinX() + (dX / 2.0)
        minY = bbox.getMinY() + (dY / 2.0)
        dCounts = (maxCounts - minCounts) / (numX * numY - 1)

        coordList = []
        counts = minCounts
        for i in range(numX):
            x = minX + (dX * i)
            for j in range(numY):
                y = minY + (dY * j)
                coordList.append([x, y, counts, sigma])
                counts += dCounts
        return coordList

    def testTempBackgrounds(self):
        """Test that the temporary backgrounds we remove are properly restored"""
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(12345, 67890), lsst.geom.Extent2I(128, 127))
        original = afwImage.ExposureF(bbox)
        rng = np.random.RandomState(123)
        original.image.array[:] = rng.normal(size=original.image.array.shape)
        original.mask.set(0)
        original.variance.set(1.0)

        def checkExposure(original, doTempLocalBackground, doTempWideBackground):
            """Check that the original exposure is unmodified."""
            config = SourceDetectionTask.ConfigClass()
            config.reEstimateBackground = False
            config.thresholdType = "pixel_stdev"
            config.doTempLocalBackground = doTempLocalBackground
            config.doTempWideBackground = doTempWideBackground
            schema = afwTable.SourceTable.makeMinimalSchema()
            task = SourceDetectionTask(config=config, schema=schema)

            exposure = original.clone()
            task.detectFootprints(exposure, sigma=3.21)

            self.assertFloatsEqual(exposure.image.array, original.image.array)
            # Mask is permitted to vary: DETECTED bit gets set
            self.assertFloatsEqual(exposure.variance.array, original.variance.array)

        checkExposure(original, False, False)
        checkExposure(original, True, False)
        checkExposure(original, False, True)
        checkExposure(original, True, True)

    def test_removeBadPixels(self):
        """Test that if we set a NO_DATA region on top of a source,
        One fewer sources is detected than when we don't do that."""
        badMaskPlanes = ["NO_DATA", ]
        exposure_regular, schema, numX, numY, starSigma = self._create_exposure()
        exposure_removed = exposure_regular.clone()

        # Define a region on one test exposure we will set to NO_DATA, blocking out one source
        region = lsst.geom.Box2I(exposure_removed.getXY0(), lsst.geom.Extent2I(25, 35))
        exposure_removed[region].mask.array |= exposure_removed.mask.getPlaneBitMask(badMaskPlanes)

        config = SourceDetectionTask.ConfigClass()
        config.reEstimateBackground = False
        config.thresholdType = "stdev"
        config.excludeMaskPlanes = badMaskPlanes
        task = SourceDetectionTask(config=config, schema=schema)

        # The regular test exposure finds 25 sources
        self._check_detectFootprints(exposure_regular, numX, numY, starSigma, task, config, doSmooth=True)

        # Modify numx and numy, which check_detectFootprints multiplies, to yield 24 instead of 25
        self.assertEqual(numX, numY)
        self._check_detectFootprints(exposure_removed, numX-1, numY+1, starSigma, task, config, doSmooth=True)

    def _create_exposure_bkg(self):
        """Return a simulated exposure for reestimating background."""

        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(1024, 1024))
        kwid = 11
        sky = 2000

        coordList = [[10, 10, 0, 2.0]]

        exposure = plantSources(
            bbox=bbox,
            kwid=kwid,
            sky=sky,
            coordList=coordList,
            addPoissonNoise=True,
        )

        return exposure

    def test_reEstimateBackground(self):
        exp, schema, _, _, _ = self._create_exposure(bigBbox=True)

        # Add in some sky.
        sky = 2000.0
        exp.image.array[:, :] += sky

        # Subtract that sky using the background task.
        bkgTask = SubtractBackgroundTask()
        result = bkgTask.run(exp)
        background = result.background

        config = SourceDetectionTask.ConfigClass()
        config.reEstimateBackground = True
        task = SourceDetectionTask(config=config, schema=schema)

        result = task.detectFootprints(exposure=exp, doSmooth=True, sigma=2.0, background=background)

        # Check that an additional background has been added to the list.
        self.assertEqual(len(background), 2)

        backgroundImage = background.getImage()
        self.assertFloatsAlmostEqual(np.mean(backgroundImage.array), sky, atol=0.1)

    def test_reEstimateBackgroundWithFlatRatio(self):
        exp, schema, _, _, _ = self._create_exposure(bigBbox=True)

        # Add in some sky.
        sky = 2000.0
        exp.image.array[:, :] += sky

        # Subtract that sky using the background task.
        bkgConfig = SubtractBackgroundTask.ConfigClass()
        bkgConfig.doApplyFlatBackgroundRatio = True
        bkgTask = SubtractBackgroundTask(config=bkgConfig)

        ratioImage = exp.image.clone()
        ratioImage.array[:, :] = 0.5

        result = bkgTask.run(exp, backgroundToPhotometricRatio=ratioImage)
        background = result.background

        config = SourceDetectionTask.ConfigClass()
        config.reEstimateBackground = True
        config.doApplyFlatBackgroundRatio = True
        task = SourceDetectionTask(config=config, schema=schema)

        result = task.detectFootprints(
            exposure=exp,
            doSmooth=True,
            sigma=2.0,
            background=background,
            backgroundToPhotometricRatio=ratioImage,
        )

        # Check that an additional background has been added to the list.
        self.assertEqual(len(background), 2)

        backgroundImage = background.getImage()
        self.assertFloatsAlmostEqual(np.mean(backgroundImage.array), sky*2, atol=0.2)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
