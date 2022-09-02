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

import astropy.time

import lsst.afw.table as afwTable
import lsst.geom
import lsst.log
from lsst.meas.algorithms import ReferenceObjectLoader, getRefFluxField, getRefFluxKeys
from lsst.meas.algorithms.loadReferenceObjects import hasNanojanskyFluxUnits, convertToNanojansky
import lsst.pex.config
import lsst.utils.tests

from ingestIndexTestBase import makeConvertConfig


class ReferenceObjectLoaderTestCase(lsst.utils.tests.TestCase):
    """Test generic parts of loader, but not the actual catalog loading."""
    def testFilterMapVsAnyFilterMapsToThis(self):
        config = ReferenceObjectLoader.ConfigClass()
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
                        ReferenceObjectLoader.makeMinimalSchema(**argDict)
                else:
                    refSchema = ReferenceObjectLoader.makeMinimalSchema(**argDict)
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
            config = ReferenceObjectLoader.ConfigClass()
            config.filterMap = filterMap
            loader = ReferenceObjectLoader(None, None, name=None, config=config)
            refSchema = ReferenceObjectLoader.makeMinimalSchema(filterNameList="r")
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
        config = ReferenceObjectLoader.ConfigClass()
        config.anyFilterMapsToThis = "gg"
        loader = ReferenceObjectLoader(None, None, name=None, config=config)
        refSchema = ReferenceObjectLoader.makeMinimalSchema(filterNameList=["gg"])
        loader._addFluxAliases(refSchema,
                               anyFilterMapsToThis=config.anyFilterMapsToThis,
                               filterMap=config.filterMap)
        self.assertEqual(getRefFluxField(refSchema, "r"), "gg_flux")
        # raise if "gg" is not in the refcat filter list
        with self.assertRaises(RuntimeError):
            refSchema = ReferenceObjectLoader.makeMinimalSchema(filterNameList=["rr"])
            refSchema = loader._addFluxAliases(refSchema,
                                               anyFilterMapsToThis=config.anyFilterMapsToThis,
                                               filterMap=config.filterMap)

    def testCheckFluxUnits(self):
        """Test that we can identify old style fluxes in a schema."""
        schema = ReferenceObjectLoader.makeMinimalSchema(['r', 'z'])
        # the default schema should pass
        self.assertTrue(hasNanojanskyFluxUnits(schema))
        schema.addField('bad_fluxSigma', doc='old flux units', type=float, units='')
        self.assertFalse(hasNanojanskyFluxUnits(schema))

        schema = ReferenceObjectLoader.makeMinimalSchema(['r', 'z'])
        schema.addField('bad_flux', doc='old flux units', type=float, units='')
        self.assertFalse(hasNanojanskyFluxUnits(schema))

        schema = ReferenceObjectLoader.makeMinimalSchema(['r', 'z'])
        schema.addField('bad_flux', doc='old flux units', type=float, units='Jy')
        self.assertFalse(hasNanojanskyFluxUnits(schema))

        schema = ReferenceObjectLoader.makeMinimalSchema(['r', 'z'])
        schema.addField('bad_fluxErr', doc='old flux units', type=float, units='')
        self.assertFalse(hasNanojanskyFluxUnits(schema))

        schema = ReferenceObjectLoader.makeMinimalSchema(['r', 'z'])
        schema.addField('bad_fluxErr', doc='old flux units', type=float, units='Jy')
        self.assertFalse(hasNanojanskyFluxUnits(schema))

        schema = ReferenceObjectLoader.makeMinimalSchema(['r', 'z'])
        schema.addField('bad_fluxSigma', doc='old flux units', type=float, units='')
        self.assertFalse(hasNanojanskyFluxUnits(schema))

    def testConvertOldFluxes(self):
        """Check that we can convert old style fluxes in a catalog."""
        flux = 1.234
        fluxErr = 5.678
        log = lsst.log.Log()

        def make_catalog():
            schema = ReferenceObjectLoader.makeMinimalSchema(['r', 'z'])
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

    def testGetMetadataCircle(self):
        center = lsst.geom.SpherePoint(100*lsst.geom.degrees, 45*lsst.geom.degrees)
        radius = lsst.geom.Angle(1*lsst.geom.degrees)
        loader = ReferenceObjectLoader(None, None, name=None)
        metadata = loader.getMetadataCircle(center, radius, "fakeR")
        self.assertEqual(metadata['RA'], center.getLongitude().asDegrees())
        self.assertEqual(metadata['DEC'], center.getLatitude().asDegrees())
        self.assertEqual(metadata['RADIUS'], radius.asDegrees())
        self.assertEqual(metadata['SMATCHV'], 2)
        self.assertEqual(metadata['FILTER'], 'fakeR')
        self.assertEqual(metadata['JEPOCH'], None)
        self.assertEqual(metadata['TIMESYS'], 'TAI')

        epoch = astropy.time.Time(2023.0, format="jyear", scale="tai")
        metadata = loader.getMetadataCircle(center, radius, "fakeR", epoch=epoch)
        self.assertEqual(metadata['JEPOCH'], epoch.jyear)


class ConvertReferenceCatalogConfigValidateTestCase(lsst.utils.tests.TestCase):
    """Test validation of ConvertReferenceCatalogConfig."""
    def testValidateRaDecMag(self):
        config = makeConvertConfig()
        config.validate()

        for name in ("ra_name", "dec_name", "mag_column_list"):
            with self.subTest(name=name):
                config = makeConvertConfig()
                setattr(config, name, None)
                with self.assertRaises(ValueError):
                    config.validate()

    def testValidateRaDecErr(self):
        # check that a basic config validates
        config = makeConvertConfig(withRaDecErr=True)
        config.validate()

        # check that a config with any of these fields missing does not validate
        for name in ("ra_err_name", "dec_err_name", "coord_err_unit"):
            with self.subTest(name=name):
                config = makeConvertConfig(withRaDecErr=True)
                setattr(config, name, None)
                with self.assertRaises(ValueError):
                    config.validate()

        # check that coord_err_unit must be an astropy unit
        config = makeConvertConfig(withRaDecErr=True)
        config.coord_err_unit = "nonsense unit"
        with self.assertRaisesRegex(ValueError, "is not a valid astropy unit string"):
            config.validate()

    def testValidateMagErr(self):
        config = makeConvertConfig(withMagErr=True)
        config.validate()

        # test for missing names
        for name in config.mag_column_list:
            with self.subTest(name=name):
                config = makeConvertConfig(withMagErr=True)
                del config.mag_err_column_map[name]
                with self.assertRaises(ValueError):
                    config.validate()

        # test for incorrect names
        for name in config.mag_column_list:
            with self.subTest(name=name):
                config = makeConvertConfig(withMagErr=True)
                config.mag_err_column_map["badName"] = config.mag_err_column_map[name]
                del config.mag_err_column_map[name]
                with self.assertRaises(ValueError):
                    config.validate()

    def testValidatePm(self):
        basicNames = ["pm_ra_name", "pm_dec_name", "epoch_name", "epoch_format", "epoch_scale"]

        for withPmErr in (False, True):
            config = makeConvertConfig(withPm=True, withPmErr=withPmErr)
            config.validate()
            del config

            if withPmErr:
                names = basicNames + ["pm_ra_err_name", "pm_dec_err_name"]
            else:
                names = basicNames
                for name in names:
                    with self.subTest(name=name, withPmErr=withPmErr):
                        config = makeConvertConfig(withPm=True, withPmErr=withPmErr)
                        setattr(config, name, None)
                        with self.assertRaises(ValueError):
                            config.validate()

    def testValidateParallax(self):
        """Validation should fail if any parallax-related fields are missing.
        """
        names = ["parallax_name", "epoch_name", "epoch_format", "epoch_scale", "parallax_err_name"]

        config = makeConvertConfig(withParallax=True)
        config.validate()
        del config

        for name in names:
            with self.subTest(name=name):
                config = makeConvertConfig(withParallax=True)
                setattr(config, name, None)
                with self.assertRaises(ValueError, msg=name):
                    config.validate()


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
