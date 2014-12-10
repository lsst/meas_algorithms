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

class ApertureFluxTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        self.exposure = lsst.afw.image.ExposureF(201, 201)
        # for convenience, we'll put the source at the origin
        self.exposure.setXY0(lsst.afw.geom.Point2I(-100,-100))
        self.psf = lsst.meas.algorithms.DoubleGaussianPsf(71, 71, 8.0, 15.0, 1.0)
        self.exposure.setPsf(self.psf)
        self.flux = 50.0
        self.radii = [3.0, 6.0, 9.0, 12.0, 15.0, 98.0, 200.0]
        psfImage = self.psf.computeImage()
        box = psfImage.getBBox(lsst.afw.image.PARENT)
        image = self.exposure.getMaskedImage().getImage()
        subImage = image.Factory(image, box, lsst.afw.image.PARENT, False)
        subImage.scaledPlus(self.flux, psfImage.convertF())
        self.footprint = lsst.afw.detection.Footprint(box)
        self.footprint.getPeaks().append(lsst.afw.detection.Peak(0,0))
        self.config = lsst.meas.algorithms.SourceMeasurementConfig()
        self.config.algorithms.names = ["flux.sinc", "flux.naive", "flux.aperture"]
        self.config.algorithms["flux.sinc"].radius = 3.0
        self.config.algorithms["flux.naive"].radius = 15.0
        self.config.algorithms["flux.aperture"].maxSincRadius = 10.0
        self.config.algorithms["flux.aperture"].radii = self.radii
        self.config.centroider.name = None
        self.config.slots.centroid = None
        self.config.slots.shape = None
        self.config.slots.apFlux = None
        self.config.slots.modelFlux = None
        self.config.slots.psfFlux = None
        self.config.slots.instFlux = None
        self.config.slots.calibFlux = None
        self.config.doReplaceWithNoise = False

    def tearDown(self):
        del self.psf
        del self.exposure
        del self.footprint

    def measure(self):
        schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        task = lsst.meas.algorithms.SourceMeasurementTask(config=self.config, schema=schema)
        catalog = lsst.afw.table.SourceCatalog(schema)
        source = catalog.addNew()
        source.setFootprint(self.footprint)
        task.run(self.exposure, catalog)
        return source

    def testApertureFlux(self):
        """Check that ApertureFlux produces values consistent with SincFlux and NaiveFlux
        for small and large radii (respectively), and that large radii include all the flux
        in the image.
        """
        source = self.measure()
        fluxArrayKey = source.schema.find("flux.aperture").key
        # we failed because the last radius was too big, but the others are ok
        self.assertTrue(source.get("flux.aperture.flags"))
        n = source.get("flux.aperture.nProfile")
        self.assertTrue(numpy.isnan(source.get(fluxArrayKey[n])))
        self.assertEqual(n, len(self.radii) - 1)
        self.assertEqual(source.get(fluxArrayKey[0]), source.get("flux.sinc"))
        self.assertEqual(source.get(fluxArrayKey[4]), source.get("flux.naive"))
        self.assertClose(source.get(fluxArrayKey[5]), self.flux, rtol=1E-7)
        for i1, i2 in zip(range(0,n-1), range(1, n)):
            self.assertLess(source.get(fluxArrayKey[i1]), source.get(fluxArrayKey[i2]))

def suite():
    """Returns a suite containing all the test cases in this module."""
    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(ApertureFluxTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    lsst.utils.tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
