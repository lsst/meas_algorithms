#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2016 LSST Corporation.
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

import unittest
import numpy
import lsst.utils.tests
import lsst.meas.base.tests
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
import lsst.pex.config
import lsst.meas.algorithms.measureApCorr as measureApCorr
from lsst.meas.base.apCorrRegistry import addApCorrName
from lsst.afw.math import ChebyshevBoundedField


def apCorrTestSourceCatalog(schema=None, name=None, apname=None, config=None, apCorrScale=1.0,
                            numSources=None, dimension=None):
    sourceCat = afwTable.SourceCatalog(schema)
    fluxName = name + "_flux"
    flagName = name + "_flag"
    fluxSigmaName = name + "_fluxSigma"
    apFluxName = apname + "_flux"
    apFlagName = apname + "_flag"
    apFluxSigmaName = apname + "_fluxSigma"
    fluxKey = schema.find(fluxName).key
    flagKey = schema.find(flagName).key
    fluxSigmaKey = schema.find(fluxSigmaName).key
    apFluxKey = schema.find(apFluxName).key
    apFlagKey = schema.find(apFlagName).key
    apFluxSigmaKey = schema.find(apFluxSigmaName).key
    centroidKey = afwTable.Point2DKey(schema["slot_Centroid"])
    inputFilterFlagKey = schema.find(config.inputFilterFlag).key
    x = numpy.random.rand(numSources) * dimension
    y = numpy.random.rand(numSources) * dimension
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


def apCorrDefaultMap(value=None, bbox=None):
    default_coefficients = numpy.ones((1, 1), dtype=float)
    default_coefficients /= value
    default_apCorrMap = ChebyshevBoundedField(bbox, default_coefficients)
    default_fill = afwImage.ImageF(bbox)
    default_apCorrMap.fillImage(default_fill)
    return(default_fill)


class MeasureApCorrTestCase(lsst.meas.base.tests.AlgorithmTestCase):

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
        task = measureApCorr.MeasureApCorrTask
        config = task.ConfigClass()
        config.refFluxName = name
        config.inputFilterFlag = calib_flag_name
        self.meas_apCorr_task = task(schema=schema, config=config)
        self.name = name
        self.apname = apname
        self.schema = schema

    def tearDown(self):
        del self.schema
        del self.meas_apCorr_task

    def testAddFields(self):
        """Instantiating the task should add one field to the schema"""
        self.assertTrue("apcorr_" + self.name + "_used" in self.schema.getNames())

    def testReturnApCorrMap(self):
        """The measureApCorr task should return a structure with a single key "apCorrMap"""
        struct = self.meas_apCorr_task.run(catalog=[], bbox=afwGeom.Box2I())
        self.assertEqual(struct.getDict().keys(), ['apCorrMap'])

    def testApCorrMapKeys(self):
        """An apCorrMap structure should have two keys, based on the name supplied to addApCorrName()"""
        apfluxName = self.apname + "_flux"
        apfluxSigmaName = self.apname + "_fluxSigma"
        struct = self.meas_apCorr_task.run(catalog=[], bbox=afwGeom.Box2I())
        key_names = [apfluxName, apfluxSigmaName]
        self.assertEqual(set(struct.apCorrMap.keys()), set(key_names))

    def testTooFewSources(self):
        """ If there are too few sources, check that a default of a field of 1s is returned."""
        apFluxName = self.apname + "_flux"
        catalog = afwTable.SourceCatalog(self.schema)
        # For this test, pick a minimal bounding box that will return an array
        bbox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.ExtentI(2, 2))
        struct = self.meas_apCorr_task.run(catalog=catalog, bbox=bbox)
        default_coefficients = numpy.ones((1, 1), dtype=float)
        default_apCorrMap = ChebyshevBoundedField(bbox, default_coefficients)
        default_fill = afwImage.ImageF(bbox)
        default_apCorrMap.fillImage(default_fill)
        test_fill = afwImage.ImageF(bbox)
        struct.apCorrMap[apFluxName].fillImage(test_fill)
        numpy.testing.assert_allclose(test_fill.getArray(), default_fill.getArray())

    def testSourceNotUsed(self):
        """ Check that a source outside the bounding box is flagged as not used (False)"""
        fluxName = self.name + "_flux"
        apCorrFlagKey = self.schema.find("apcorr_" + self.name + "_used").key
        catalog = afwTable.SourceCatalog(self.schema)
        source = catalog.addNew()

        source_test_flux = 5.1
        source_test_centroid = afwGeom.Point2D(5, 7.1)
        fluxKey = self.schema.find(fluxName).key
        centroidKey = afwTable.Point2DKey(self.schema["slot_Centroid"])
        source.set(fluxKey, source_test_flux)
        source.set(centroidKey, source_test_centroid)
        # For this test, pick a minimal bounding box that will return an array
        bbox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.ExtentI(2, 2))
        self.meas_apCorr_task.run(catalog=catalog, bbox=bbox)
        self.assertFalse(catalog[apCorrFlagKey])

    def testSourceUsed(self):
        """ Check that valid sources inside the bounding box that are used have their flags set to True"""
        bbox_size = 10
        sourceCat = apCorrTestSourceCatalog(schema=self.schema, name=self.name, apname=self.apname,
                                            config=self.meas_apCorr_task.config, apCorrScale=1.0,
                                            numSources=5, dimension=bbox_size)
        bbox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.ExtentI(bbox_size, bbox_size))
        inputFilterFlagKey = self.schema.find(self.meas_apCorr_task.config.inputFilterFlag).key
        self.meas_apCorr_task.run(catalog=sourceCat, bbox=bbox)
        self.assertTrue(sourceCat[inputFilterFlagKey].all())

    def testApertureMeasOnes(self):
        """ Check that sources with aperture fluxes exactly the same as their catalog fluxes
            returns an aperture correction map of 1s"""
        bbox_size = 10
        sourceCat = apCorrTestSourceCatalog(schema=self.schema, name=self.name, apname=self.apname,
                                            config=self.meas_apCorr_task.config, apCorrScale=1.0,
                                            numSources=5, dimension=bbox_size)
        apFluxName = self.apname + "_flux"

        bbox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.ExtentI(bbox_size, bbox_size))
        struct = self.meas_apCorr_task.run(catalog=sourceCat, bbox=bbox)
        default_fill = apCorrDefaultMap(value=1.0, bbox=bbox)
        test_fill = afwImage.ImageF(bbox)
        struct.apCorrMap[apFluxName].fillImage(test_fill)
        numpy.testing.assert_allclose(test_fill.getArray(), default_fill.getArray())

    def testApertureMeasTens(self):
        """ Check that aperture correction scales source fluxes in the correct direction"""
        apCorr_factor = 10.
        bbox_size = 10
        sourceCat = apCorrTestSourceCatalog(schema=self.schema, name=self.name, apname=self.apname,
                                            config=self.meas_apCorr_task.config, apCorrScale=apCorr_factor,
                                            numSources=5, dimension=bbox_size)

        apFluxName = self.apname + "_flux"
        bbox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.ExtentI(bbox_size, bbox_size))
        struct = self.meas_apCorr_task.run(catalog=sourceCat, bbox=bbox)
        default_fill = apCorrDefaultMap(value=apCorr_factor, bbox=bbox)
        test_fill = afwImage.ImageF(bbox)
        struct.apCorrMap[apFluxName].fillImage(test_fill)
        numpy.testing.assert_allclose(test_fill.getArray(), default_fill.getArray())


def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(MeasureApCorrTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)


def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
