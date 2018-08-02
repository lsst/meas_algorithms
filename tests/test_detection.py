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
from lsst.meas.algorithms import SourceDetectionTask
from lsst.meas.algorithms.testUtils import plantSources
import lsst.utils.tests

display = False


class DetectionTestCase(lsst.utils.tests.TestCase):
    """Test the aperture correction."""

    def testBasics(self):
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(256, 100), lsst.geom.Extent2I(128, 127), invert=False)
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
        config = SourceDetectionTask.ConfigClass()
        config.reEstimateBackground = False
        task = SourceDetectionTask(config=config, schema=schema)
        for doSmooth in (False, True):
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
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(12345, 67890), lsst.geom.Extent2I(128, 127), invert=False)
        original = afwImage.ExposureF(bbox)
        rng = np.random.RandomState(123)
        original.image.array[:] = rng.normal(size=original.image.array.shape)
        original.mask.set(0)
        original.variance.set(1.0)

        def checkExposure(original, doTempLocalBackground, doTempWideBackground):
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


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
