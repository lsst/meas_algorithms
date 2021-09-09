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

import os
import unittest
from collections import Counter

import astropy.time
import astropy.units
import numpy as np

import lsst.geom
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
import lsst.daf.persistence as dafPersist
from lsst.meas.algorithms import (IngestIndexedReferenceTask, LoadIndexedReferenceObjectsTask,
                                  LoadIndexedReferenceObjectsConfig, getRefFluxField)
from lsst.meas.algorithms.loadReferenceObjects import hasNanojanskyFluxUnits
import lsst.utils

from ingestIndexTestBase import (makeConvertConfig, ConvertReferenceCatalogTestBase,
                                 make_coord)

REGENERATE_COMPARISON = False  # Regenerate comparison data?


class IngestIndexTaskValidateTestCase(lsst.utils.tests.TestCase):
    """Test validation of IngestIndexReferenceConfig."""
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


class ReferenceCatalogIngestAndLoadTestCase(ConvertReferenceCatalogTestBase, lsst.utils.tests.TestCase):
    """Tests of converting, ingesting, loading and validating an HTM Indexed
    Reference Catalog (gen2 code path).
    """
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.obs_test_dir = lsst.utils.getPackageDir('obs_test')
        cls.input_dir = os.path.join(cls.obs_test_dir, "data", "input")

        # Run the ingest once to create a butler repo we can compare to
        config = makeConvertConfig(withMagErr=True, withRaDecErr=True, withPm=True, withPmErr=True,
                                   withParallax=True)
        # Pregenerated gen2 test refcats have the "cal_ref_cat" name.
        config.dataset_config.ref_dataset_name = "cal_ref_cat"
        config.dataset_config.indexer.active.depth = cls.depth
        config.id_name = 'id'
        config.pm_scale = 1000.0  # arcsec/yr --> mas/yr
        config.parallax_scale = 1e3  # arcsec -> milliarcsec
        # np.savetxt prepends '# ' to the header lines, so use a reader that understands that
        config.file_reader.format = 'ascii.commented_header'
        IngestIndexedReferenceTask.parseAndRun(args=[cls.input_dir, "--output", cls.testRepoPath,
                                                     cls.skyCatalogFile], config=config)
        cls.testButler = dafPersist.Butler(cls.testRepoPath)

    @classmethod
    def tearDownClass(cls):
        del cls.testButler

    def testSanity(self):
        """Sanity-check that compCats contains some entries with sources."""
        numWithSources = 0
        for idList in self.compCats.values():
            if len(idList) > 0:
                numWithSources += 1
        self.assertGreater(numWithSources, 0)

    def testAgainstPersisted(self):
        """Test that we can get a specific shard from a pre-persisted refcat.
        """
        shardId = 2222
        dataset_name = IngestIndexedReferenceTask.ConfigClass().dataset_config.ref_dataset_name
        dataId = self.indexer.makeDataId(shardId, dataset_name)
        self.assertTrue(self.testButler.datasetExists('ref_cat', dataId))
        refCat = self.testButler.get('ref_cat', dataId)
        if REGENERATE_COMPARISON:
            if os.path.exists(self.testCatPath):
                os.unlink(self.testCatPath)
            refCat.writeFits(self.testCatPath)
            self.fail("New comparison data written; unset REGENERATE_COMPARISON in order to proceed")

        ex1 = refCat.extract('*')
        testCat = afwTable.SimpleCatalog.readFits(self.testCatPath)

        ex2 = testCat.extract('*')
        self.assertEqual(set(ex1.keys()), set(ex2.keys()))
        for key in ex1:
            np.testing.assert_array_almost_equal(ex1[key], ex2[key], err_msg=f"{key} values not equal")

    def testIngestSetsVersion(self):
        """Test that newly ingested catalogs get the correct version number set.
        """
        def runTest(withRaDecErr):
            outputPath = os.path.join(self.outPath, "output_setsVersion"
                                      + "_withRaDecErr" if withRaDecErr else "")
            # Test with multiple files and standard config
            config = makeConvertConfig(withRaDecErr=withRaDecErr, withMagErr=True,
                                       withPm=True, withPmErr=True)
            # Pregenerated gen2 test refcats have the "cal_ref_cat" name.
            config.dataset_config.ref_dataset_name = "cal_ref_cat"
            # don't use the default depth, to avoid taking the time to create thousands of file locks
            config.dataset_config.indexer.active.depth = self.depth
            IngestIndexedReferenceTask.parseAndRun(
                args=[self.input_dir, "--output", outputPath, self.skyCatalogFile],
                config=config)
            # A newly-ingested refcat should be marked format_version=1.
            loader = LoadIndexedReferenceObjectsTask(butler=dafPersist.Butler(outputPath))
            self.assertEqual(loader.dataset_config.format_version, 1)

        runTest(withRaDecErr=True)
        runTest(withRaDecErr=False)

    def testIngestConfigOverrides(self):
        """Test IngestIndexedReferenceTask with different configs.
        """
        config2 = makeConvertConfig(withRaDecErr=True, withMagErr=True, withPm=True, withPmErr=True,
                                    withParallax=True)
        config2.ra_name = "ra"
        config2.dec_name = "dec"
        config2.dataset_config.ref_dataset_name = 'myrefcat'
        # Change the indexing depth to prove we can.
        # Smaller is better than larger because it makes fewer files.
        config2.dataset_config.indexer.active.depth = self.depth - 1
        config2.is_photometric_name = 'is_phot'
        config2.is_resolved_name = 'is_res'
        config2.is_variable_name = 'is_var'
        config2.id_name = 'id'
        config2.extra_col_names = ['val1', 'val2', 'val3']
        config2.file_reader.header_lines = 1
        config2.file_reader.colnames = [
            'id', 'ra', 'dec', 'ra_err', 'dec_err', 'a', 'a_err', 'b', 'b_err', 'is_phot',
            'is_res', 'is_var', 'val1', 'val2', 'val3', 'pm_ra', 'pm_dec', 'pm_ra_err',
            'pm_dec_err', 'parallax', 'parallax_err', 'unixtime',
        ]
        config2.file_reader.delimiter = '|'
        # this also tests changing the delimiter
        IngestIndexedReferenceTask.parseAndRun(
            args=[self.input_dir, "--output", self.outPath+"/output_override",
                  self.skyCatalogFileDelim], config=config2)

        # Test if we can get back the catalog with a non-standard dataset name
        butler = dafPersist.Butler(self.outPath+"/output_override")
        loaderConfig = LoadIndexedReferenceObjectsConfig()
        loaderConfig.ref_dataset_name = "myrefcat"
        loader = LoadIndexedReferenceObjectsTask(butler=butler, config=loaderConfig)
        self.checkAllRowsInRefcat(loader, self.skyCatalog, config2)

        # TODO: this test is probably irrelevant in gen3, since the name is now the collection.
        # test that a catalog can be loaded even with a name not used for ingestion
        butler = dafPersist.Butler(self.testRepoPath)
        loaderConfig2 = LoadIndexedReferenceObjectsConfig()
        loaderConfig2.ref_dataset_name = self.testDatasetName
        loader = LoadIndexedReferenceObjectsTask(butler=butler, config=loaderConfig2)
        self.checkAllRowsInRefcat(loader, self.skyCatalog, config2)

    def testLoadIndexedReferenceConfig(self):
        """Make sure LoadIndexedReferenceConfig has needed fields."""
        """
        Including at least one from the base class LoadReferenceObjectsConfig
        """
        config = LoadIndexedReferenceObjectsConfig()
        self.assertEqual(config.ref_dataset_name, "cal_ref_cat")
        self.assertEqual(config.defaultFilter, "")

    def testLoadSkyCircle(self):
        """Test LoadIndexedReferenceObjectsTask.loadSkyCircle with default config."""
        loader = LoadIndexedReferenceObjectsTask(butler=self.testButler)
        for tupl, idList in self.compCats.items():
            cent = make_coord(*tupl)
            lcat = loader.loadSkyCircle(cent, self.searchRadius, filterName='a')
            self.assertTrue(lcat.refCat.isContiguous())
            self.assertFalse("camFlux" in lcat.refCat.schema)
            self.assertEqual(Counter(lcat.refCat['id']), Counter(idList))
            if len(lcat.refCat) > 0:
                # make sure there are no duplicate ids
                self.assertEqual(len(set(Counter(lcat.refCat['id']).values())), 1)
                self.assertEqual(len(set(Counter(idList).values())), 1)
                # A default-loaded sky circle should not have centroids
                self.assertNotIn("centroid_x", lcat.refCat.schema)
                self.assertNotIn("centroid_y", lcat.refCat.schema)
                self.assertNotIn("hasCentroid", lcat.refCat.schema)
            else:
                self.assertEqual(len(idList), 0)

    def testLoadPixelBox(self):
        """Test LoadIndexedReferenceObjectsTask.loadPixelBox with default config."""
        loader = LoadIndexedReferenceObjectsTask(butler=self.testButler)
        numFound = 0
        for tupl, idList in self.compCats.items():
            cent = make_coord(*tupl)
            bbox = lsst.geom.Box2I(lsst.geom.Point2I(30, -5), lsst.geom.Extent2I(1000, 1004))  # arbitrary
            ctr_pix = bbox.getCenter()
            # catalog is sparse, so set pixel scale such that bbox encloses region
            # used to generate compCats
            pixel_scale = 2*self.searchRadius/max(bbox.getHeight(), bbox.getWidth())
            cdMatrix = afwGeom.makeCdMatrix(scale=pixel_scale)
            wcs = afwGeom.makeSkyWcs(crval=cent, crpix=ctr_pix, cdMatrix=cdMatrix)
            result = loader.loadPixelBox(bbox=bbox, wcs=wcs, filterName="a")
            # The following is to ensure the reference catalog coords are
            # getting corrected for proper motion when an epoch is provided.
            # Use an extreme epoch so that differences in corrected coords
            # will be significant.  Note that this simply tests that the coords
            # do indeed change when the epoch is passed.  It makes no attempt
            # at assessing the correctness of the change.  This is left to the
            # explicit testProperMotion() test below.
            resultWithEpoch = loader.loadPixelBox(bbox=bbox, wcs=wcs, filterName="a",
                                                  epoch=astropy.time.Time(30000, format='mjd', scale="tai"))
            self.assertFloatsNotEqual(result.refCat["coord_ra"], resultWithEpoch.refCat["coord_ra"],
                                      rtol=1.0e-4)
            self.assertFloatsNotEqual(result.refCat["coord_dec"], resultWithEpoch.refCat["coord_dec"],
                                      rtol=1.0e-4)
            self.assertFalse("camFlux" in result.refCat.schema)
            self.assertGreaterEqual(len(result.refCat), len(idList))
            numFound += len(result.refCat)
        self.assertGreater(numFound, 0)

    def testDefaultFilterAndFilterMap(self):
        """Test defaultFilter and filterMap parameters of LoadIndexedReferenceObjectsConfig."""
        config = LoadIndexedReferenceObjectsConfig()
        config.defaultFilter = "b"
        config.filterMap = {"aprime": "a"}
        loader = LoadIndexedReferenceObjectsTask(butler=self.testButler, config=config)
        for tupl, idList in self.compCats.items():
            cent = make_coord(*tupl)
            lcat = loader.loadSkyCircle(cent, self.searchRadius)
            self.assertEqual(lcat.fluxField, "camFlux")
            if len(idList) > 0:
                defFluxFieldName = getRefFluxField(lcat.refCat.schema, None)
                self.assertTrue(defFluxFieldName in lcat.refCat.schema)
                aprimeFluxFieldName = getRefFluxField(lcat.refCat.schema, "aprime")
                self.assertTrue(aprimeFluxFieldName in lcat.refCat.schema)
                break  # just need one test

    def testProperMotion(self):
        """Test proper motion correction"""
        center = make_coord(93.0, -90.0)
        loader = LoadIndexedReferenceObjectsTask(butler=self.testButler)
        references = loader.loadSkyCircle(center, self.searchRadius, filterName='a').refCat
        original = references.copy(True)

        # Zero epoch change --> no proper motion correction (except minor numerical effects)
        loader.applyProperMotions(references, self.epoch)
        self.assertFloatsAlmostEqual(references["coord_ra"], original["coord_ra"], rtol=1.0e-14)
        self.assertFloatsAlmostEqual(references["coord_dec"], original["coord_dec"], rtol=1.0e-14)
        self.assertFloatsEqual(references["coord_raErr"], original["coord_raErr"])
        self.assertFloatsEqual(references["coord_decErr"], original["coord_decErr"])

        # One year difference
        loader.applyProperMotions(references, self.epoch + 1.0*astropy.units.yr)
        self.assertFloatsEqual(references["pm_raErr"], original["pm_raErr"])
        self.assertFloatsEqual(references["pm_decErr"], original["pm_decErr"])
        for orig, ref in zip(original, references):
            self.assertAnglesAlmostEqual(orig.getCoord().separation(ref.getCoord()),
                                         self.properMotionAmt, maxDiff=1.0e-6*lsst.geom.arcseconds)
            self.assertAnglesAlmostEqual(orig.getCoord().bearingTo(ref.getCoord()),
                                         self.properMotionDir, maxDiff=1.0e-4*lsst.geom.arcseconds)
        predictedRaErr = np.hypot(original["coord_raErr"], original["pm_raErr"])
        predictedDecErr = np.hypot(original["coord_decErr"], original["pm_decErr"])
        self.assertFloatsAlmostEqual(references["coord_raErr"], predictedRaErr)
        self.assertFloatsAlmostEqual(references["coord_decErr"], predictedDecErr)

    def testRequireProperMotion(self):
        """Tests of the requireProperMotion config field.

        Requiring proper motion corrections for a catalog that does not
        contain valid PM data should result in an exception.

        `data/testHtmIndex-ps1-bad-pm.fits` is a random shard taken from the
        ps1_pv3_3pi_20170110 refcat (that has the unitless PM fields),
        stripped to only 2 rows: we patch it in here to simplify test setup.
        """
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data/testHtmIndex-ps1-bad-pm.fits')
        refcatData = lsst.afw.table.SimpleCatalog.readFits(path)
        center = make_coord(93.0, -90.0)
        epoch = self.epoch + 1.0*astropy.units.yr

        # malformatted catalogs should warn and raise if we require proper motion corrections
        config = LoadIndexedReferenceObjectsConfig()
        config.requireProperMotion = True
        config.anyFilterMapsToThis = "g"  # to use a catalog not made for obs_test
        loader = LoadIndexedReferenceObjectsTask(butler=self.testButler, config=config)
        with unittest.mock.patch.object(self.testButler, 'get', return_value=refcatData):
            msg = "requireProperMotion=True but refcat pm_ra field is not an Angle"
            with self.assertRaisesRegex(RuntimeError, msg):
                loader.loadSkyCircle(center, self.searchRadius, epoch=epoch)

        # not specifying `epoch` with requireProperMotion=True should raise for any catalog
        config = LoadIndexedReferenceObjectsConfig()
        config.requireProperMotion = True
        config.anyFilterMapsToThis = "g"  # to use a catalog not made for obs_test
        loader = LoadIndexedReferenceObjectsTask(butler=self.testButler, config=config)
        msg = "requireProperMotion=True but epoch not provided to loader"
        with self.assertRaisesRegex(RuntimeError, msg):
            loader.loadSkyCircle(center, self.searchRadius, epoch=None)

        # malformatted catalogs should just warn if we do not require proper motion corrections
        config = LoadIndexedReferenceObjectsConfig()
        config.requireProperMotion = False
        config.anyFilterMapsToThis = "g"  # to use a catalog not made for obs_test
        loader = LoadIndexedReferenceObjectsTask(butler=self.testButler, config=config)
        with unittest.mock.patch.object(self.testButler, 'get', return_value=refcatData):
            with lsst.log.UsePythonLogging(), self.assertLogs(level="WARNING") as cm:
                loader.loadSkyCircle(center, self.searchRadius, epoch=epoch)
            warnLog1 = "Reference catalog pm_ra field is not an Angle; cannot apply proper motion."
            self.assertEqual(cm.output, [f"WARNING:LoadIndexedReferenceObjectsTask:{warnLog1}"])

    def testLoadVersion0(self):
        """Test reading a pre-written format_version=0 (Jy flux) catalog.
        It should be converted to have nJy fluxes.
        """
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data/version0')
        loader = LoadIndexedReferenceObjectsTask(butler=dafPersist.Butler(path))
        self.assertEqual(loader.dataset_config.format_version, 0)
        result = loader.loadSkyCircle(make_coord(10, 20),
                                      5*lsst.geom.degrees, filterName='a')
        self.assertTrue(hasNanojanskyFluxUnits(result.refCat.schema))
        catalog = afwTable.SimpleCatalog.readFits(os.path.join(path, 'ref_cats/cal_ref_cat/4022.fits'))
        self.assertFloatsEqual(catalog['a_flux']*1e9, result.refCat['a_flux'])
        self.assertFloatsEqual(catalog['a_fluxSigma']*1e9, result.refCat['a_fluxErr'])
        self.assertFloatsEqual(catalog['b_flux']*1e9, result.refCat['b_flux'])
        self.assertFloatsEqual(catalog['b_fluxSigma']*1e9, result.refCat['b_fluxErr'])

    def testLoadVersion1(self):
        """Test reading a format_version=1 catalog (fluxes unchanged)."""
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data/version1')
        loader = LoadIndexedReferenceObjectsTask(butler=dafPersist.Butler(path))
        self.assertEqual(loader.dataset_config.format_version, 1)
        result = loader.loadSkyCircle(make_coord(10, 20),
                                      5*lsst.geom.degrees, filterName='a')
        self.assertTrue(hasNanojanskyFluxUnits(result.refCat.schema))
        catalog = afwTable.SimpleCatalog.readFits(os.path.join(path, 'ref_cats/cal_ref_cat/4022.fits'))
        self.assertFloatsEqual(catalog['a_flux'], result.refCat['a_flux'])
        self.assertFloatsEqual(catalog['a_fluxErr'], result.refCat['a_fluxErr'])
        self.assertFloatsEqual(catalog['b_flux'], result.refCat['b_flux'])
        self.assertFloatsEqual(catalog['b_fluxErr'], result.refCat['b_fluxErr'])


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
