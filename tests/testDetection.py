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
from builtins import range
import unittest

import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
from lsst.meas.algorithms import SourceDetectionTask
from lsst.meas.algorithms.testUtils import plantSources
import lsst.utils.tests

display = False


class DetectionTestCase(lsst.utils.tests.TestCase):
    """Test the aperture correction."""

    def testBasics(self):
        bbox = afwGeom.Box2I(afwGeom.Point2I(256, 100), afwGeom.Extent2I(128, 127))
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
            self.assertEqual(task.metadata.get("sigma"), taskSigma)
            self.assertEqual(task.metadata.get("doSmooth"), doSmooth)
            self.assertEqual(task.metadata.get("nGrow"), int(taskSigma * config.nSigmaToGrow + 0.5))

            res = task.detectFootprints(exposure, doSmooth=doSmooth, sigma=None)
            taskSigma = task.metadata.get("sigma")
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


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
