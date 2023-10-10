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
import unittest
import numpy as np

import lsst.geom
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
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
    return default_fill


class MeasureApCorrTestCase(lsst.meas.base.tests.AlgorithmTestCase, lsst.utils.tests.TestCase):

    def makeCatalog(self, apCorrScale=1.0, numSources=5):
        sourceCat = afwTable.SourceCatalog(self.schema)
        inputFilterFlagKey = self.schema.find(self.meas_apCorr_task.config.sourceSelector.active.field).key

        centroidKey = afwTable.Point2DKey(self.schema["slot_Centroid"])
        x = np.random.rand(numSources)*self.exposure.getWidth() + self.exposure.getX0()
        y = np.random.rand(numSources)*self.exposure.getHeight() + self.exposure.getY0()
        for _i in range(numSources):
            source_test_instFlux = 5.1
            source_test_centroid = lsst.geom.Point2D(x[_i], y[_i])
            source = sourceCat.addNew()
            source.set(centroidKey, source_test_centroid)
            source.set(inputFilterFlagKey, True)

        for name in self.names:
            fluxName = name + "_instFlux"
            flagName = name + "_flag"
            fluxErrName = name + "_instFluxErr"
            apFluxName = name + self.apNameStr + "_instFlux"
            apFlagName = name + self.apNameStr + "_flag"
            apFluxErrName = name + self.apNameStr + "_instFluxErr"
            fluxKey = self.schema.find(fluxName).key
            flagKey = self.schema.find(flagName).key
            fluxErrKey = self.schema.find(fluxErrName).key
            apFluxKey = self.schema.find(apFluxName).key
            apFlagKey = self.schema.find(apFlagName).key
            apFluxErrKey = self.schema.find(apFluxErrName).key
            for source in sourceCat:
                source.set(fluxKey, source_test_instFlux)
                source.set(apFluxKey, source_test_instFlux * apCorrScale)
                source.set(fluxErrKey, 0.)
                source.set(apFluxErrKey, 0.)
                source.set(flagKey, False)
                source.set(apFlagKey, False)
        return sourceCat

    def setUp(self):
        schema = afwTable.SourceTable.makeMinimalSchema()
        apNameStr = "Ap"
        calib_flag_name = "cal_source_use"
        # Add fields in anti-sorted order to try to impose a need for sorting
        # in the addition of the apCorr fields (may happen by fluke, but this
        # is the best we can do to test this here.
        names = ["test2", "test1"]
        for name in names:
            apName = name + apNameStr
            addApCorrName(apName)
            schema.addField(name + "_instFlux", type=float)
            schema.addField(name + "_instFluxErr", type=float)
            schema.addField(name + "_flag", type="Flag")
            schema.addField(apName + "_instFlux", type=float)
            schema.addField(apName + "_instFluxErr", type=float)
            schema.addField(apName + "_flag", type="Flag")
        schema.addField(names[0] + "_Centroid_x", type=float)
        schema.addField(names[0] + "_Centroid_y", type=float)
        schema.getAliasMap().set("slot_Centroid", names[0] + "_Centroid")
        schema.addField(calib_flag_name, type="Flag")
        config = measureApCorr.MeasureApCorrTask.ConfigClass()
        config.refFluxName = names[0]
        config.sourceSelector.active.field = calib_flag_name
        self.meas_apCorr_task = measureApCorr.MeasureApCorrTask(schema=schema, config=config)
        self.names = names
        self.apNameStr = apNameStr
        self.schema = schema
        self.exposure = lsst.afw.image.ExposureF(10, 10)

    def tearDown(self):
        del self.schema
        del self.meas_apCorr_task
        del self.exposure

    def testAddFields(self):
        """Instantiating the task should add one field to the schema."""
        for name in self.names:
            self.assertIn("apcorr_" + name + self.apNameStr + "_used", self.schema.getNames())
        sortedNames = sorted(self.names)
        key0 = self.schema.find("apcorr_" + sortedNames[0] + self.apNameStr + "_used").key
        key1 = self.schema.find("apcorr_" + sortedNames[1] + self.apNameStr + "_used").key
        # Check that the apCorr fields were added in a sorted order (not
        # foolproof as this could have happened by fluke, but it's the best
        # we can do to test this here (having added the two fields in an anti-
        # sorted order).
        self.assertLess(key0.getOffset() + key0.getBit(), key1.getOffset() + key1.getBit())

    def testReturnApCorrMap(self):
        """The measureApCorr task should return a structure with a single key 'apCorrMap'."""
        struct = self.meas_apCorr_task.run(catalog=self.makeCatalog(), exposure=self.exposure)
        self.assertEqual(list(struct.getDict().keys()), ['apCorrMap'])

    def testApCorrMapKeys(self):
        """An apCorrMap structure should have two keys per name supplied to addApCorrName()."""
        key_names = []
        for name in self.names:
            apFluxName = name + self.apNameStr + "_instFlux"
            apFluxErrName = name + self.apNameStr + "_instFluxErr"
            struct = self.meas_apCorr_task.run(catalog=self.makeCatalog(), exposure=self.exposure)
            key_names.append(apFluxName)
            key_names.append(apFluxErrName)
        self.assertEqual(set(struct.apCorrMap.keys()), set(key_names))

    def testTooFewSources(self):
        """ If there are too few sources, check that an exception is raised."""
        catalog = afwTable.SourceCatalog(self.schema)
        with self.assertRaises(RuntimeError):
            self.meas_apCorr_task.run(catalog=catalog, exposure=self.exposure)
        # With the measurement algorithm declared as something that might fail, should not get an exception
        for name in self.names:
            self.meas_apCorr_task.config.allowFailure.append(name + self.apNameStr)
        self.meas_apCorr_task.run(catalog=catalog, exposure=self.exposure)

    def testSourceNotUsed(self):
        """ Check that a source outside the bounding box is flagged as not used (False)."""
        fluxName = self.names[0] + "_instFlux"
        apCorrFlagKey = self.schema.find("apcorr_" + self.names[0] + "_used").key
        sourceCat = self.makeCatalog()
        source = sourceCat.addNew()
        source_test_instFlux = 5.1
        source_test_centroid = lsst.geom.Point2D(15, 7.1)
        fluxKey = self.schema.find(fluxName).key
        centroidKey = afwTable.Point2DKey(self.schema["slot_Centroid"])
        source.set(fluxKey, source_test_instFlux)
        source.set(centroidKey, source_test_centroid)
        self.meas_apCorr_task.run(catalog=sourceCat, exposure=self.exposure)
        self.assertFalse(sourceCat[apCorrFlagKey][-1])

    def testSourceUsed(self):
        """Check that valid sources inside the bounding box that are used have their flags set to True."""
        inputFilterFlagKey = self.schema.find(self.meas_apCorr_task.config.sourceSelector.active.field).key
        sourceCat = self.makeCatalog()
        self.meas_apCorr_task.run(catalog=sourceCat, exposure=self.exposure)
        self.assertTrue(sourceCat[inputFilterFlagKey].all())

    def testApertureMeasOnes(self):
        """ Check that sources with aperture fluxes exactly the same as their catalog fluxes
            returns an aperture correction map of 1s"""
        apFluxName = self.names[0] + self.apNameStr + "_instFlux"
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
        apFluxName = self.names[0] + self.apNameStr + "_instFlux"
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
