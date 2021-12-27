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

import itertools
import unittest

import lsst.afw.table as afwTable
import lsst.log
from lsst.meas.algorithms import ReferenceObjectLoaderBase, getRefFluxField, getRefFluxKeys
from lsst.meas.algorithms.loadReferenceObjects import hasNanojanskyFluxUnits, convertToNanojansky
import lsst.pex.config
import lsst.utils.tests


class TrivialLoader(ReferenceObjectLoaderBase):
    """Minimal subclass of ReferenceObjectLoaderBase to allow instantiation
    """

    def loadSkyCircle(self, ctrCoord, radius, filterName):
        pass


class TestReferenceObjectLoaderBase(lsst.utils.tests.TestCase):
    """Test case for ReferenceObjectLoaderBase abstract base class.

    Only methods with concrete implementations are tested (hence not loadSkyCircle)
    """

    def testFilterMapVsAnyFilterMapsToThis(self):
        config = TrivialLoader.ConfigClass()
        # check that a filterMap-only config passes validation
        config.filterMap = {"b": "a"}
        try:
            config.validate()
        except lsst.pex.config.FieldValidationError:
            self.fail("`filterMap`-only LoadReferenceObjectsConfig should not fail validation.")

        # anyFilterMapsToThis and filterMap are mutually exclusive
        config.anyFilterMapsToThis = "c"
        with self.assertRaises(lsst.pex.config.FieldValidationError):
            config.validate()

        # check that a anyFilterMapsToThis-only config passes validation
        config.filterMap = {}
        try:
            config.validate()
        except lsst.pex.config.FieldValidationError:
            self.fail("`anyFilterMapsToThis`-only LoadReferenceObjectsConfig should not fail validation.")

    def testMakeMinimalSchema(self):
        """Make a schema and check it."""
        for filterNameList in (["r"], ["foo", "_bar"]):
            for (addIsPhotometric, addIsResolved, addIsVariable,
                 coordErrDim, addProperMotion, properMotionErrDim,
                 addParallax) in itertools.product(
                    (False, True), (False, True), (False, True),
                    (-1, 0, 1, 2, 3, 4), (False, True), (-1, 0, 1, 2, 3, 4),
                    (False, True)):
                argDict = dict(
                    filterNameList=filterNameList,
                    addIsPhotometric=addIsPhotometric,
                    addIsResolved=addIsResolved,
                    addIsVariable=addIsVariable,
                    coordErrDim=coordErrDim,
                    addProperMotion=addProperMotion,
                    properMotionErrDim=properMotionErrDim,
                    addParallax=addParallax,
                )
                if coordErrDim not in (0, 2, 3) or \
                        (addProperMotion and properMotionErrDim not in (0, 2, 3)):
                    with self.assertRaises(ValueError):
                        ReferenceObjectLoaderBase.makeMinimalSchema(**argDict)
                else:
                    refSchema = ReferenceObjectLoaderBase.makeMinimalSchema(**argDict)
                    self.assertTrue("coord_ra" in refSchema)
                    self.assertTrue("coord_dec" in refSchema)
                    for filterName in filterNameList:
                        fluxField = filterName + "_flux"
                        self.assertIn(fluxField, refSchema)
                        self.assertNotIn("x" + fluxField, refSchema)
                        fluxErrField = fluxField + "Err"
                        self.assertIn(fluxErrField, refSchema)
                        self.assertEqual(getRefFluxField(refSchema, filterName), filterName + "_flux")
                    self.assertEqual("resolved" in refSchema, addIsResolved)
                    self.assertEqual("variable" in refSchema, addIsVariable)
                    self.assertEqual("photometric" in refSchema, addIsPhotometric)
                    self.assertEqual("photometric" in refSchema, addIsPhotometric)
                    self.assertEqual("epoch" in refSchema, addProperMotion or addParallax)
                    self.assertEqual("coord_raErr" in refSchema, coordErrDim > 0)
                    self.assertEqual("coord_decErr" in refSchema, coordErrDim > 0)
                    self.assertEqual("coord_ra_dec_Cov" in refSchema, coordErrDim == 3)
                    self.assertEqual("pm_ra" in refSchema, addProperMotion)
                    self.assertEqual("pm_dec" in refSchema, addProperMotion)
                    self.assertEqual("pm_raErr" in refSchema, addProperMotion and properMotionErrDim > 0)
                    self.assertEqual("pm_decErr" in refSchema, addProperMotion and properMotionErrDim > 0)
                    self.assertEqual("pm_flag" in refSchema, addProperMotion)
                    self.assertEqual("pm_ra_dec_Cov" in refSchema,
                                     addProperMotion and properMotionErrDim == 3)
                    self.assertEqual("parallax" in refSchema, addParallax)
                    self.assertEqual("parallaxErr" in refSchema, addParallax)
                    self.assertEqual("parallax_flag" in refSchema, addParallax)

    def testFilterAliasMap(self):
        """Make a schema with filter aliases."""
        for filterMap in ({}, {"camr": "r"}):
            config = TrivialLoader.ConfigClass()
            config.filterMap = filterMap
            loader = TrivialLoader(config=config)
            refSchema = TrivialLoader.makeMinimalSchema(filterNameList="r")
            loader._addFluxAliases(refSchema,
                                   anyFilterMapsToThis=config.anyFilterMapsToThis,
                                   filterMap=config.filterMap)

            self.assertIn("r_flux", refSchema)
            self.assertIn("r_fluxErr", refSchema)

            # camera filters aliases are named <filter>_camFlux
            if "camr" in filterMap:
                self.assertEqual(getRefFluxField(refSchema, "camr"), "camr_camFlux")
            else:
                with self.assertRaisesRegex(RuntimeError,
                                            r"Could not find flux field\(s\) camr_camFlux, camr_flux"):
                    getRefFluxField(refSchema, "camr")

            refCat = afwTable.SimpleCatalog(refSchema)
            refObj = refCat.addNew()
            refObj["r_flux"] = 1.23
            self.assertAlmostEqual(refCat[0].get(getRefFluxField(refSchema, "r")), 1.23)
            if "camr" in filterMap:
                self.assertAlmostEqual(refCat[0].get(getRefFluxField(refSchema, "camr")), 1.23)
            refObj["r_fluxErr"] = 0.111
            if "camr" in filterMap:
                self.assertEqual(refCat[0].get("camr_camFluxErr"), 0.111)
            fluxKey, fluxErrKey = getRefFluxKeys(refSchema, "r")
            self.assertEqual(refCat[0].get(fluxKey), 1.23)
            self.assertEqual(refCat[0].get(fluxErrKey), 0.111)
            if "camr" in filterMap:
                fluxKey, fluxErrKey = getRefFluxKeys(refSchema, "camr")
                self.assertEqual(refCat[0].get(fluxErrKey), 0.111)
            else:
                with self.assertRaises(RuntimeError):
                    getRefFluxKeys(refSchema, "camr")

    def testAnyFilterMapsToThisAlias(self):
        # test anyFilterMapsToThis
        config = TrivialLoader.ConfigClass()
        config.anyFilterMapsToThis = "gg"
        loader = TrivialLoader(config=config)
        refSchema = TrivialLoader.makeMinimalSchema(filterNameList=["gg"])
        loader._addFluxAliases(refSchema,
                               anyFilterMapsToThis=config.anyFilterMapsToThis,
                               filterMap=config.filterMap)
        self.assertEqual(getRefFluxField(refSchema, "r"), "gg_flux")
        # raise if "gg" is not in the refcat filter list
        with self.assertRaises(RuntimeError):
            refSchema = TrivialLoader.makeMinimalSchema(filterNameList=["rr"])
            refSchema = loader._addFluxAliases(refSchema,
                                               anyFilterMapsToThis=config.anyFilterMapsToThis,
                                               filterMap=config.filterMap)

    def testCheckFluxUnits(self):
        """Test that we can identify old style fluxes in a schema."""
        schema = ReferenceObjectLoaderBase.makeMinimalSchema(['r', 'z'])
        # the default schema should pass
        self.assertTrue(hasNanojanskyFluxUnits(schema))
        schema.addField('bad_fluxSigma', doc='old flux units', type=float, units='')
        self.assertFalse(hasNanojanskyFluxUnits(schema))

        schema = ReferenceObjectLoaderBase.makeMinimalSchema(['r', 'z'])
        schema.addField('bad_flux', doc='old flux units', type=float, units='')
        self.assertFalse(hasNanojanskyFluxUnits(schema))

        schema = ReferenceObjectLoaderBase.makeMinimalSchema(['r', 'z'])
        schema.addField('bad_flux', doc='old flux units', type=float, units='Jy')
        self.assertFalse(hasNanojanskyFluxUnits(schema))

        schema = ReferenceObjectLoaderBase.makeMinimalSchema(['r', 'z'])
        schema.addField('bad_fluxErr', doc='old flux units', type=float, units='')
        self.assertFalse(hasNanojanskyFluxUnits(schema))

        schema = ReferenceObjectLoaderBase.makeMinimalSchema(['r', 'z'])
        schema.addField('bad_fluxErr', doc='old flux units', type=float, units='Jy')
        self.assertFalse(hasNanojanskyFluxUnits(schema))

        schema = ReferenceObjectLoaderBase.makeMinimalSchema(['r', 'z'])
        schema.addField('bad_fluxSigma', doc='old flux units', type=float, units='')
        self.assertFalse(hasNanojanskyFluxUnits(schema))

    def testConvertOldFluxes(self):
        """Check that we can convert old style fluxes in a catalog."""
        flux = 1.234
        fluxErr = 5.678
        log = lsst.log.Log()

        def make_catalog():
            schema = ReferenceObjectLoaderBase.makeMinimalSchema(['r', 'z'])
            schema.addField('bad_flux', doc='old flux units', type=float, units='')
            schema.addField('bad_fluxErr', doc='old flux units', type=float, units='Jy')
            refCat = afwTable.SimpleCatalog(schema)
            refObj = refCat.addNew()
            refObj["bad_flux"] = flux
            refObj["bad_fluxErr"] = fluxErr
            return refCat

        oldRefCat = make_catalog()
        newRefCat = convertToNanojansky(oldRefCat, log)
        self.assertEqual(newRefCat['bad_flux'], [flux*1e9, ])
        self.assertEqual(newRefCat['bad_fluxErr'], [fluxErr*1e9, ])
        self.assertEqual(newRefCat.schema['bad_flux'].asField().getUnits(), 'nJy')
        self.assertEqual(newRefCat.schema['bad_fluxErr'].asField().getUnits(), 'nJy')

        # check that doConvert=False returns None (it also logs a summary)
        oldRefCat = make_catalog()
        newRefCat = convertToNanojansky(oldRefCat, log, doConvert=False)
        self.assertIsNone(newRefCat)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
