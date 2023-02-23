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

import unittest

import lsst.geom
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.daf.base as dafBase
from lsst.meas.algorithms import SourceDetectionTask
from lsst.meas.base import SingleFrameMeasurementTask as SourceMeasurementTask
from lsst.meas.algorithms.testUtils import plantSources
import lsst.utils.tests

try:
    display
except NameError:
    display = False
else:
    import lsst.afw.display as afwDisplay
    afwDisplay.setDefaultMaskTransparency(75)


class NegativeMeasurementTestCase(lsst.utils.tests.TestCase):
    """Testing detection and measurement on negative objects.
    """

    def _create_exposure(self):
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(256, 100), lsst.geom.Extent2I(128, 127))
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
        return exposure, numX, numY

    def test_detection_stdev(self):
        """Test detection and measurement on an exposure with negative sources
        for thresholdType="stdev".
        """
        exposure, numX, numY = self._create_exposure()

        if display:
            disp = afwDisplay.Display(frame=1)
            disp.mtv(exposure, title=self._testMethodName + ": image with -ve sources")

        schema = afwTable.SourceTable.makeMinimalSchema()
        config = SourceDetectionTask.ConfigClass()
        config.reEstimateBackground = False
        config.thresholdPolarity = 'both'
        detection = SourceDetectionTask(config=config, schema=schema)
        algMetadata = dafBase.PropertyList()
        measurement = SourceMeasurementTask(schema=schema, algMetadata=algMetadata)

        table = afwTable.SourceTable.make(schema)
        detections = detection.run(table, exposure)
        sources = detections.sources

        self.assertEqual(len(sources), numX*numY)
        self.assertEqual(detections.numPos, numX*numY/2)
        self.assertEqual(detections.numNeg, numX*numY/2)

        measurement.run(sources, exposure)

        nGoodCent = 0
        nGoodShape = 0
        for s in sources:
            cent = s.getCentroid()
            shape = s.getShape()

            if cent[0] == cent[0] and cent[1] == cent[1]:
                nGoodCent += 1

            if (shape.getIxx() == shape.getIxx()
                    and shape.getIyy() == shape.getIyy()
                    and shape.getIxy() == shape.getIxy()):
                nGoodShape += 1

            if display:
                xy = cent[0], cent[1]
                disp.dot('+', *xy)
                disp.dot(shape, *xy, ctype=afwDisplay.RED)

        self.assertEqual(nGoodCent, numX*numY)
        self.assertEqual(nGoodShape, numX*numY)

    def test_significance(self):
        """Test that negative peaks have the right significance for
        thresholdType='stdev' for the non-convolved, non-local-background case.
        """
        exposure, numX, numY = self._create_exposure()

        schema = afwTable.SourceTable.makeMinimalSchema()
        config = SourceDetectionTask.ConfigClass()
        config.thresholdPolarity = 'both'
        # don't modify the image after detection.
        config.reEstimateBackground = False
        config.doTempLocalBackground = False
        detection = SourceDetectionTask(config=config, schema=schema)
        result = detection.detectFootprints(exposure, doSmooth=False)

        bad = exposure.mask.getPlaneBitMask(config.statsMask)
        sctrl = afwMath.StatisticsControl()
        sctrl.setAndMask(bad)
        stats = afwMath.makeStatistics(exposure.maskedImage, afwMath.STDEVCLIP, sctrl)
        stddev = stats.getValue(afwMath.STDEVCLIP)
        # Don't bother checking positive footprints: those are tested more
        # thoroughly in test_detection.py.
        for footprint in result.negative.getFootprints():
            for peak in footprint.peaks:
                point = lsst.geom.Point2I(peak.getIx(), peak.getIy())
                value = exposure.image[point]
                with self.subTest(str(point)):
                    self.assertFloatsAlmostEqual(peak["significance"],
                                                 -value/stddev,  # S/N for negative peak
                                                 rtol=1e-7,
                                                 msg=str(point))

    def makeCoordList(self, bbox, numX, numY, minCounts, maxCounts, sigma):
        """Make a coordList for makeExposure.
        """
        dX = bbox.getWidth()/float(numX)
        dY = bbox.getHeight()/float(numY)
        minX = bbox.getMinX() + (dX/2.0)
        minY = bbox.getMinY() + (dY/2.0)
        dCounts = (maxCounts - minCounts)/(numX*numY/2 - 1)

        coordList = []
        counts = minCounts
        for i in range(numX):
            x = minX + (dX*i)
            for j in range(numY):
                y = minY + (dY*j)
                if j%2 == 0:
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
