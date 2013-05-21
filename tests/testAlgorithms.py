#!/usr/bin/env python

#
# LSST Data Management System
# Copyright 2008-2013 LSST Corporation.
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

import numpy
import unittest

import lsst.pex.exceptions as pexEx
import lsst.pex.logging as pexLog
import lsst.afw.coord as afwCoord
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.meas.algorithms as measAlg
import lsst.utils.tests as utilsTests
import lsst.afw.detection.detectionLib as afwDetection
import testLib

import lsst.afw.display.ds9 as ds9

try:
    type(display)
except NameError:
    display = False

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class AlgorithmsTestCase(unittest.TestCase):
    """A test case for centroiding"""

    @utilsTests.debugger(pexEx.LsstCppException, Exception, AssertionError)
    def testAlgorithms(self):
        """Test that we can instantiate and use algorithms"""

        config = measAlg.SourceMeasurementConfig()
        config.algorithms.names = measAlg.AlgorithmRegistry.all.keys()
        config.algorithms.names.discard(config.centroider.name)
        config.doReplaceWithNoise = False

        config.algorithms.names.discard("flux.peakLikelihood")

        if False:
            log = pexLog.getDefaultLog()
            log.setThreshold(log.DEBUG)

        schema = afwTable.SourceTable.makeMinimalSchema()
        task = measAlg.SourceMeasurementTask(schema, config=config)
        catalog = afwTable.SourceCatalog(schema)
        source = catalog.addNew()
        source.set("id", 12345)

        size = 128
        xStar, yStar = 65.432, 76.543
        width = 3.21
        x0, y0 = 12345, 54321
        x, y = numpy.indices((size, size))
        im = afwImage.MaskedImageF(afwGeom.ExtentI(size, size))
        im.setXY0(afwGeom.Point2I(x0, y0))
        im.getVariance().set(1.0)
        arr = im.getImage().getArray()
        arr[y,x] = numpy.exp(-0.5*((x - xStar)**2 + (y - yStar)**2)/width**2)
        psf = testLib.makeTestPsf(im)
        exp = afwImage.makeExposure(im)
        exp.setPsf(psf)
        exp.setXY0(afwGeom.Point2I(x0, y0))
        scale = 1.0e-5
        wcs = afwImage.makeWcs(afwCoord.Coord(0.0*afwGeom.degrees, 0.0*afwGeom.degrees),
                               afwGeom.Point2D(0.0, 0.0), scale, 0.0, 0.0, scale)
        exp.setWcs(wcs)

        peak = afwDetection.Peak(int(xStar + x0), int(yStar + y0), 1.0)
        bbox = im.getBBox()
        bbox.shift(afwGeom.Extent2I(x0, y0))
        foot = afwDetection.Footprint(peak.getI(), width, bbox)
        afwDetection.setMaskFromFootprint(exp.getMaskedImage().getMask(), foot, 1)
        foot.getPeaks().push_back(peak)
        source.setFootprint(foot)

        if display:
            ds9.mtv(exp, frame=1)

        task.run(exp, catalog)

        for alg in config.algorithms:
            flagName = alg + ".flags"
            if False:
                print (alg, source.get(flagName) if flagName in schema else None,
                       source.get(alg) if alg in schema else None)
            elif flagName in schema:
                self.assertFalse(source.get(alg + ".flags"))

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(AlgorithmsTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the tests"""
    utilsTests.run(suite(), exit)
 
if __name__ == "__main__":
    run(True)
