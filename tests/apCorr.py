#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2014 LSST Corporation.
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

"""
Tests for aperture correction code
"""

import numpy
import unittest
import lsst.utils.tests as utilsTests

import lsst.afw.detection as afwDetection
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.display.ds9 as ds9

from lsst.meas.algorithms import SourceMeasurementConfig, MeasureApCorrTask

try:
    type(verbose)
except NameError:
    display = False


class ApCorrTestCase(utilsTests.TestCase):
    def test(self):
        """Check that we can measure an aperture correction"""
        width, height = 1234, 2345
        edge = 25 # buffer around edge
        numStars = 10 # in each dimension
        size = 15 # Size of PSF
        sigma = 3.0 # PSF sigma

        measureConfig = SourceMeasurementConfig()
        measureConfig.algorithms["flux.sinc"].radius = size
        apCorrConfig = MeasureApCorrTask.ConfigClass()
        apCorrConfig.inputFilterFlag = "use.this"

        schema = afwTable.SourceTable.makeMinimalSchema()
        measurer = measureConfig.makeMeasureSources(schema)
        schema.addField("use.this", type="Flag", doc="use this for aperture correction")
        apCorr = MeasureApCorrTask(config=apCorrConfig, schema=schema)

        # Generate an image with a bunch of sources
        exposure = afwImage.ExposureF(width, height)
        psf = afwDetection.GaussianPsf(size, size, sigma)
        exposure.setPsf(psf)
        image = exposure.getMaskedImage().getImage()
        mask = exposure.getMaskedImage().getMask()
        mask.addMaskPlane("DETECTED")
        detected = mask.getPlaneBitMask("DETECTED")
        xList = numpy.linspace(edge, width - edge, numStars)
        yList = numpy.linspace(edge, height - edge, numStars)
        for x in xList:
            for y in yList:
                psfImage = psf.computeImage(afwGeom.Point2D(x, y)).convertF()
                psfImage *= 10000
                bbox = psfImage.getBBox(afwImage.PARENT)
                subImage = image.Factory(image, bbox)
                subImage += psfImage

        # Detect and measure those sources
        catalog = afwTable.SourceCatalog(schema)
        measureConfig.slots.setupTable(catalog.table)
        feet = afwDetection.FootprintSet(exposure.getMaskedImage(), afwDetection.Threshold(0.1), "DETECTED")
        feet.makeSources(catalog)
        self.assertGreater(len(catalog), 0)
        for source in catalog:
            measurer.applyWithPeak(source, exposure)
            source.set("use.this", True)

        if display:
            ds9.mtv(exposure, frame=1)

        # Cheat in order to trigger clipping in aperture correction measurement
        good = numpy.ones(len(catalog), dtype=bool)
        good[len(catalog)//3] = False
        catalog[apCorrConfig.reference][numpy.where(good == 0)] *= 2

        apCorrMap = apCorr.run(exposure.getBBox(), catalog)

        # Validate answers
        self.assertGreater(len(apCorrMap), 0) # Or we're not really testing anything
        for name in apCorrMap:
            expected = (catalog[apCorrConfig.reference][good]/catalog[name][good]).mean()
            corr = apCorrMap[name]
            for x in xList:
                for y in yList:
                    value = corr.evaluate(x, y)
                    if name.endswith(".err"):
                        self.assertLess(value, 5.0e-4)
                    else:
                        self.assertClose(corr.evaluate(x, y), expected, atol=5.0e-4)

            if not name.endswith(".err"):
                self.assertTrue(numpy.all(catalog["apcorr." + name + ".used"] == good)) # Used only good srcs

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(ApCorrTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the utilsTests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
