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

class FluxTestCase(unittest.TestCase):

    def assertClose(self, a, b, rtol=1E-5, atol=1E-8):
        self.assert_(numpy.allclose(a, b, rtol=rtol, atol=atol), "\n%s\n!=\n%s" % (a, b))

    def assertNotClose(self, a, b, rtol=1E-5, atol=1E-8):
        self.assertFalse(numpy.allclose(a, b, rtol=rtol, atol=atol), "\n%s\n==\n%s" % (a, b))

    def setUp(self):
        self.exposure = lsst.afw.image.ExposureF(201, 201)
        # for convenience, we'll put the source at the origin
        self.exposure.setXY0(lsst.afw.geom.Point2I(-100,-100))
        self.exposure.getMaskedImage().getVariance()[:] = 1.0
        self.psf = lsst.meas.algorithms.DoubleGaussianPsf(71, 71, 8.0, 15.0, 1.0)
        self.exposure.setPsf(self.psf)
        self.flux = 50.0
        psfImage = self.psf.computeImage()
        box = psfImage.getBBox(lsst.afw.image.PARENT)
        image = self.exposure.getMaskedImage().getImage()
        subImage = image.Factory(image, box, lsst.afw.image.PARENT, False)
        subImage.scaledPlus(self.flux, psfImage.convertF())
        self.footprint = lsst.afw.detection.Footprint(box)
        self.footprint.getPeaks().append(lsst.afw.detection.Peak(0,0))
        self.config = lsst.meas.algorithms.SourceMeasurementConfig()
        self.config.algorithms.names = ["flux.psf", "flux.gaussian", "flux.sinc", "correctfluxes"]
        self.config.algorithms.names.add("shape.sdss")
        self.config.centroider.name = None
        self.config.slots.centroid = None
        self.config.slots.shape = None
        self.config.doReplaceWithNoise = False

    def tearDown(self):
        del self.psf
        del self.exposure
        del self.footprint

    def measure(self, radius=None):
        if radius is not None:
            self.config.algorithms["flux.sinc"].radius = radius
            self.config.algorithms["correctfluxes"].apCorrRadius = radius
        schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        task = lsst.meas.algorithms.SourceMeasurementTask(config=self.config, schema=schema)
        catalog = lsst.afw.table.SourceCatalog(schema)
        source = catalog.addNew()
        source.setFootprint(self.footprint)
        task.run(self.exposure, catalog)
        # flux.psf.psffactor should be 1.0 because it's just the dot product of the PSF with itself
        self.assertClose(source.get("flux.psf.psffactor"), 1.0)
        # flux.gaussian.psffactor should be the dot product of a Gaussian with a double-Gaussian PSF.
        self.assertNotClose(source.get("flux.gaussian.psffactor"), 1.0, rtol=1E-2, atol=1E-2)
        return source

    def testGaussian(self):
        """Test that we can measure a Gaussian flux"""
        
        self.config.algorithms["correctfluxes"].doApCorr = False
        
        self.config.algorithms["flux.gaussian"].fixed = False
        source = self.measure()
        self.assertClose(self.flux, source.get("flux.gaussian"))

        self.config.algorithms["flux.gaussian"].fixed = True
        source = self.measure()
        self.assertClose(self.flux, source.get("flux.gaussian"))
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
