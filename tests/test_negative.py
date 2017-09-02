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
import lsst.daf.base as dafBase
from lsst.meas.algorithms import SourceDetectionTask
from lsst.meas.base import SingleFrameMeasurementTask as SourceMeasurementTask
from lsst.meas.algorithms.testUtils import plantSources
import lsst.utils.tests

try:
    type(display)
    import lsst.afw.display.ds9 as ds9
except NameError:
    display = False


class NegativeMeasurementTestCase(lsst.utils.tests.TestCase):
    """A test case for negative objects."""

    def testBasics(self):
        bbox = afwGeom.Box2I(afwGeom.Point2I(256, 100), afwGeom.Extent2I(128, 127))
        minCounts = 2000
        maxCounts = 20000
        starSigma = 1.5
        numX = 4
        numY = 4
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

        if display:
            ds9.mtv(exposure)

        schema = afwTable.SourceTable.makeMinimalSchema()
        config = SourceDetectionTask.ConfigClass()
        config.reEstimateBackground = False
        config.thresholdPolarity = 'both'
        detection = SourceDetectionTask(config=config, schema=schema)
        algMetadata = dafBase.PropertyList()
        measurement = SourceMeasurementTask(schema=schema, algMetadata=algMetadata)

        table = afwTable.SourceTable.make(schema)
        detections = detection.makeSourceCatalog(table, exposure)
        sources = detections.sources
        fpSets = detections.fpSets

        self.assertEqual(len(sources), numX * numY)
        self.assertEqual(fpSets.numPos, numX * numY / 2)
        self.assertEqual(fpSets.numNeg, numX * numY / 2)

        measurement.run(sources, exposure)

        nGoodCent = 0
        nGoodShape = 0
        for s in sources:
            cent = s.getCentroid()
            shape = s.getShape()

            if cent[0] == cent[0] and cent[1] == cent[1]:
                nGoodCent += 1

            if (shape.getIxx() == shape.getIxx() and
                shape.getIyy() == shape.getIyy() and
                    shape.getIxy() == shape.getIxy()):
                nGoodShape += 1

            if display:
                xy = cent[0] - exposure.getX0(), cent[1] - exposure.getY0()
                ds9.dot('+', *xy)
                ds9.dot(shape, *xy, ctype=ds9.RED)

        self.assertEqual(nGoodCent, numX * numY)
        self.assertEqual(nGoodShape, numX * numY)

    def makeCoordList(self, bbox, numX, numY, minCounts, maxCounts, sigma):
        """Make a coordList for makeExposure."""
        dX = bbox.getWidth() / float(numX)
        dY = bbox.getHeight() / float(numY)
        minX = bbox.getMinX() + (dX / 2.0)
        minY = bbox.getMinY() + (dY / 2.0)
        dCounts = (maxCounts - minCounts) / (numX * numY / 2 - 1)

        coordList = []
        counts = minCounts
        for i in range(numX):
            x = minX + (dX * i)
            for j in range(numY):
                y = minY + (dY * j)
                if j % 2 == 0:
                    coordList.append([x, y, counts, sigma])
                else:
                    coordList.append([x, y, -counts, sigma])
                    counts += dCounts
        return coordList


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
