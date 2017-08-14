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
import numpy as np

import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
from lsst.afw.math import ChebyshevBoundedField
import lsst.pex.config
import lsst.meas.algorithms.measureApCorr as measureApCorr
from lsst.meas.base.apCorrRegistry import addApCorrName
import lsst.meas.base.tests
import lsst.utils.tests


def apCorrDefaultMap(value=None, bbox=None):
    default_coefficients = np.ones((1, 1), dtype=float)
    default_coefficients /= value
    default_apCorrMap = ChebyshevBoundedField(bbox, default_coefficients)
    default_fill = afwImage.ImageF(bbox)
    default_apCorrMap.fillImage(default_fill)
    return(default_fill)


class MeasureApCorrTestCase(lsst.meas.base.tests.AlgorithmTestCase, lsst.utils.tests.TestCase):

    def makeCatalog(self, apCorrScale=1.0, numSources=5):
        sourceCat = afwTable.SourceCatalog(self.schema)
        fluxName = self.name + "_flux"
        flagName = self.name + "_flag"
        fluxSigmaName = self.name + "_fluxSigma"
        apFluxName = self.apname + "_flux"
        apFlagName = self.apname + "_flag"
        apFluxSigmaName = self.apname + "_fluxSigma"
        fluxKey = self.schema.find(fluxName).key
        flagKey = self.schema.find(flagName).key
        fluxSigmaKey = self.schema.find(fluxSigmaName).key
        apFluxKey = self.schema.find(apFluxName).key
        apFlagKey = self.schema.find(apFlagName).key
        apFluxSigmaKey = self.schema.find(apFluxSigmaName).key
        centroidKey = afwTable.Point2DKey(self.schema["slot_Centroid"])
        inputFilterFlagKey = self.schema.find(self.meas_apCorr_task.config.starSelector.active.field).key
        x = np.random.rand(numSources)*self.exposure.getWidth() + self.exposure.getX0()
        y = np.random.rand(numSources)*self.exposure.getHeight() + self.exposure.getY0()
        for _i in range(numSources):
            source_test_flux = 5.1
            source_test_centroid = afwGeom.Point2D(x[_i], y[_i])
            source = sourceCat.addNew()
            source.set(fluxKey, source_test_flux)
            source.set(apFluxKey, source_test_flux * apCorrScale)
            source.set(centroidKey, source_test_centroid)
            source.set(fluxSigmaKey, 0.)
            source.set(apFluxSigmaKey, 0.)
            source.set(flagKey, False)
            source.set(apFlagKey, False)
            source.set(inputFilterFlagKey, True)
        return(sourceCat)

    def setUp(self):
        schema = afwTable.SourceTable.makeMinimalSchema()
        name = "test"
        apname = "testAp"
        calib_flag_name = "cal_source_use"
        addApCorrName(apname)
        schema.addField(name + "_flux", type=float)
        schema.addField(name + "_fluxSigma", type=float)
        schema.addField(name + "_flag", type="Flag")
        schema.addField(apname + "_flux", type=float)
        schema.addField(apname + "_fluxSigma", type=float)
        schema.addField(apname + "_flag", type="Flag")
        schema.addField(calib_flag_name, type="Flag")
        schema.addField(name + "_Centroid_x", type=float)
        schema.addField(name + "_Centroid_y", type=float)
        schema.getAliasMap().set('slot_Centroid', name + '_Centroid')
        config = measureApCorr.MeasureApCorrTask.ConfigClass()
        config.refFluxName = name
        config.starSelector.active.field = calib_flag_name
        self.meas_apCorr_task = measureApCorr.MeasureApCorrTask(schema=schema, config=config)
        self.name = name
        self.apname = apname
        self.schema = schema
        self.exposure = lsst.afw.image.ExposureF(10, 10)

    def apCorrDefaultMap(value=None, bbox=None):
        default_coefficients = np.ones((1, 1), dtype=float)
        default_coefficients /= value
        default_apCorrMap = ChebyshevBoundedField(bbox, default_coefficients)
        default_fill = afwImage.ImageF(bbox)
        default_apCorrMap.fillImage(default_fill)
        return(default_fill)

    def tearDown(self):
        del self.schema
        del self.meas_apCorr_task
        del self.exposure

    def testAddFields(self):
        """Instantiating the task should add one field to the schema."""
        self.assertIn("apcorr_" + self.name + "_used", self.schema.getNames())

    def testReturnApCorrMap(self):
        """The measureApCorr task should return a structure with a single key 'apCorrMap'."""
        struct = self.meas_apCorr_task.run(catalog=self.makeCatalog(), exposure=self.exposure)
        self.assertEqual(list(struct.getDict().keys()), ['apCorrMap'])

    def testApCorrMapKeys(self):
        """An apCorrMap structure should have two keys, based on the name supplied to addApCorrName()."""
        apfluxName = self.apname + "_flux"
        apfluxSigmaName = self.apname + "_fluxSigma"
        struct = self.meas_apCorr_task.run(catalog=self.makeCatalog(), exposure=self.exposure)
        key_names = [apfluxName, apfluxSigmaName]
        self.assertEqual(set(struct.apCorrMap.keys()), set(key_names))

    def testTooFewSources(self):
        """ If there are too few sources, check that an exception is raised."""
        catalog = afwTable.SourceCatalog(self.schema)
        with self.assertRaises(RuntimeError):
            self.meas_apCorr_task.run(catalog=catalog, exposure=self.exposure)
        # With the measurement algorithm declared as something that might fail, should not get an exception
        self.meas_apCorr_task.config.allowFailure.append(self.apname)
        self.meas_apCorr_task.run(catalog=catalog, exposure=self.exposure)

    def testSourceNotUsed(self):
        """ Check that a source outside the bounding box is flagged as not used (False)."""
        fluxName = self.name + "_flux"
        apCorrFlagKey = self.schema.find("apcorr_" + self.name + "_used").key
        sourceCat = self.makeCatalog()
        source = sourceCat.addNew()
        source_test_flux = 5.1
        source_test_centroid = afwGeom.Point2D(15, 7.1)
        fluxKey = self.schema.find(fluxName).key
        centroidKey = afwTable.Point2DKey(self.schema["slot_Centroid"])
        source.set(fluxKey, source_test_flux)
        source.set(centroidKey, source_test_centroid)
        self.meas_apCorr_task.run(catalog=sourceCat, exposure=self.exposure)
        self.assertFalse(sourceCat[apCorrFlagKey][-1])

    def testSourceUsed(self):
        """Check that valid sources inside the bounding box that are used have their flags set to True."""
        inputFilterFlagKey = self.schema.find(self.meas_apCorr_task.config.starSelector.active.field).key
        sourceCat = self.makeCatalog()
        self.meas_apCorr_task.run(catalog=sourceCat, exposure=self.exposure)
        self.assertTrue(sourceCat[inputFilterFlagKey].all())

    def testApertureMeasOnes(self):
        """ Check that sources with aperture fluxes exactly the same as their catalog fluxes
            returns an aperture correction map of 1s"""
        apFluxName = self.apname + "_flux"
        sourceCat = self.makeCatalog()
        struct = self.meas_apCorr_task.run(catalog=sourceCat, exposure=self.exposure)
        default_fill = apCorrDefaultMap(value=1.0, bbox=self.exposure.getBBox())
        test_fill = afwImage.ImageF(self.exposure.getBBox())
        struct.apCorrMap[apFluxName].fillImage(test_fill)
        np.testing.assert_allclose(test_fill.getArray(), default_fill.getArray())

    def testApertureMeasTens(self):
        """Check that aperture correction scales source fluxes in the correct direction."""
        apCorr_factor = 10.
        sourceCat = self.makeCatalog(apCorrScale=apCorr_factor)
        apFluxName = self.apname + "_flux"
        struct = self.meas_apCorr_task.run(catalog=sourceCat, exposure=self.exposure)
        default_fill = apCorrDefaultMap(value=apCorr_factor, bbox=self.exposure.getBBox())
        test_fill = afwImage.ImageF(self.exposure.getBBox())
        struct.apCorrMap[apFluxName].fillImage(test_fill)
        np.testing.assert_allclose(test_fill.getArray(), default_fill.getArray())


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
