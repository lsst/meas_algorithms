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

import lsst.afw.geom
import lsst.afw.image
import lsst.afw.detection
import lsst.afw.table
import lsst.meas.algorithms
import lsst.utils.tests

class FluxTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        self.exposure = lsst.afw.image.ExposureF(201, 201)
        # for convenience, we'll put the source at the origin
        self.exposure.setXY0(lsst.afw.geom.Point2I(-100,-100))
        self.exposure.getMaskedImage().getVariance()[:] = 1.0
        self.psf = lsst.meas.algorithms.SingleGaussianPsf(71, 71, 8.0)
        self.exposure.setPsf(self.psf)
        self.flux = 50.0
        psfImage = self.psf.computeImage()
        box = psfImage.getBBox(lsst.afw.image.PARENT)
        image = self.exposure.getMaskedImage().getImage()
        subImage = image.Factory(image, box, lsst.afw.image.PARENT, False)
        subImage.scaledPlus(self.flux, psfImage.convertF())
        self.footprint = lsst.afw.detection.Footprint(box)
        self.footprint.addPeak(0, 0, float("NaN"))
        self.config = lsst.meas.algorithms.SourceMeasurementConfig()
        self.config.algorithms.names = ["flux.psf", "flux.gaussian", "flux.sinc"]
        self.config.algorithms.names.add("shape.sdss")
        self.config.centroider.name = None
        self.config.slots.centroid = None
        self.config.slots.shape = None
        self.config.slots.calibFlux = None
        self.config.doReplaceWithNoise = False

    def tearDown(self):
        del self.psf
        del self.exposure
        del self.footprint

    def measure(self, radius=None):
        if radius is not None:
            self.config.algorithms["flux.sinc"].radius = radius
        schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        task = lsst.meas.algorithms.SourceMeasurementTask(config=self.config, schema=schema)
        catalog = lsst.afw.table.SourceCatalog(schema)
        source = catalog.addNew()
        source.setFootprint(self.footprint)
        task.run(self.exposure, catalog)
        return source

    def testGaussian(self):
        """Test that we can measure a Gaussian flux"""

        self.config.algorithms["flux.gaussian"].fixed = False
        source = self.measure()
        self.assertClose(self.flux, source.get("flux.gaussian"), rtol=1E-4)

        self.config.algorithms["flux.gaussian"].fixed = True
        source = self.measure()
        self.assertClose(self.flux, source.get("flux.gaussian"), rtol=1E-4)
        self.assertTrue(numpy.isfinite(source.get("flux.gaussian.err")))

def suite():
    """Returns a suite containing all the test cases in this module."""
    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(FluxTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    lsst.utils.tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
