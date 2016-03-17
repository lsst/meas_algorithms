import unittest

import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.utils.tests as utilsTests
from lsst.meas.algorithms import SourceDetectionTask
from lsst.meas.base import SingleFrameMeasurementTask as SourceMeasurementTask
from lsst.meas.algorithms.testUtils import plantSources
import lsst.daf.base as dafBase

import lsst.afw.display.ds9 as ds9

try:
    type(display)
except NameError:
    display = False

class NegativeMeasurementTestCase(unittest.TestCase):
    """A test case for negative objects"""
    def testBasics(self):
        bbox = afwGeom.Box2I(afwGeom.Point2I(256, 100), afwGeom.Extent2I(128, 127))
        minCounts = 2000
        maxCounts = 20000
        starSigma = 1.5
        numX = 4
        numY = 4
        coordList = self.makeCoordList(
            bbox = bbox,
            numX = numX,
            numY = numY,
            minCounts = minCounts,
            maxCounts = maxCounts,
            sigma = starSigma,
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
        """Make a coordList for makeExposure"""
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

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(NegativeMeasurementTestCase)

    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the tests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
