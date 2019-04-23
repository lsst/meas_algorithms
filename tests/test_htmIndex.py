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

import math
import os
import tempfile
import shutil
import unittest
import string
from collections import Counter

import astropy.time
import numpy as np

import lsst.geom
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
import lsst.daf.persistence as dafPersist
from lsst.meas.algorithms import (IngestIndexedReferenceTask, LoadIndexedReferenceObjectsTask,
                                  LoadIndexedReferenceObjectsConfig, getRefFluxField)
from lsst.meas.algorithms import IndexerRegistry
from lsst.meas.algorithms.loadReferenceObjects import hasNanojanskyFluxUnits
import lsst.utils

OBS_TEST_DIR = lsst.utils.getPackageDir('obs_test')
INPUT_DIR = os.path.join(OBS_TEST_DIR, "data", "input")

REGENERATE_COMPARISON = False  # Regenerate comparison data?


def make_coord(ra, dec):
    """Make an ICRS coord given its RA, Dec in degrees."""
    return lsst.geom.SpherePoint(ra, dec, lsst.geom.degrees)


class HtmIndexTestCase(lsst.utils.tests.TestCase):
    @classmethod
    def make_skyCatalog(cls, outPath, size=1000):
        np.random.seed(123)
        ident = np.arange(1, size+1, dtype=int)
        ra = np.random.random(size)*360.
        dec = np.degrees(np.arccos(2.*np.random.random(size) - 1.))
        dec -= 90.
        ra_err = np.ones(size)*0.1  # arcsec
        dec_err = np.ones(size)*0.1  # arcsec
        a_mag = 16. + np.random.random(size)*4.
        a_mag_err = 0.01 + np.random.random(size)*0.2
        b_mag = 17. + np.random.random(size)*5.
        b_mag_err = 0.02 + np.random.random(size)*0.3
        is_photometric = np.random.randint(2, size=size)
        is_resolved = np.random.randint(2, size=size)
        is_variable = np.random.randint(2, size=size)
        extra_col1 = np.random.normal(size=size)
        extra_col2 = np.random.normal(1000., 100., size=size)
        # compute proper motion and PM error in arcseconds/year
        # and let the ingest task scale them to radians
        pm_amt_arcsec = cls.properMotionAmt.asArcseconds()
        pm_dir_rad = cls.properMotionDir.asRadians()
        pm_ra = np.ones(size)*pm_amt_arcsec*math.cos(pm_dir_rad)
        pm_dec = np.ones(size)*pm_amt_arcsec*math.sin(pm_dir_rad)
        pm_ra_err = np.ones(size)*cls.properMotionErr.asArcseconds()*abs(math.cos(pm_dir_rad))
        pm_dec_err = np.ones(size)*cls.properMotionErr.asArcseconds()*abs(math.sin(pm_dir_rad))
        unixtime = np.ones(size)*cls.epoch.unix

        def get_word(word_len):
            return "".join(np.random.choice([s for s in string.ascii_letters], word_len))
        extra_col3 = np.array([get_word(num) for num in np.random.randint(11, size=size)])

        dtype = np.dtype([('id', float), ('ra_icrs', float), ('dec_icrs', float),
                         ('ra_err', float), ('dec_err', float), ('a', float),
                          ('a_err', float), ('b', float), ('b_err', float), ('is_phot', int),
                          ('is_res', int), ('is_var', int), ('val1', float), ('val2', float),
                          ('val3', '|S11'), ('pm_ra', float), ('pm_dec', float), ('pm_ra_err', float),
                          ('pm_dec_err', float), ('unixtime', float)])

        arr = np.array(list(zip(ident, ra, dec, ra_err, dec_err, a_mag, a_mag_err, b_mag, b_mag_err,
                                is_photometric, is_resolved, is_variable, extra_col1, extra_col2, extra_col3,
                                pm_ra, pm_dec, pm_ra_err, pm_dec_err, unixtime)), dtype=dtype)
        # write the data with full precision; this is not realistic for
        # real catalogs, but simplifies tests based on round tripped data
        saveKwargs = dict(
            header="id,ra_icrs,dec_icrs,ra_err,dec_err,"
                   "a,a_err,b,b_err,is_phot,is_res,is_var,val1,val2,val3,"
                   "pm_ra,pm_dec,pm_ra_err,pm_dec_err,unixtime",
            fmt=["%i", "%.15g", "%.15g", "%.15g", "%.15g",
                 "%.15g", "%.15g", "%.15g", "%.15g", "%i", "%i", "%i", "%.15g", "%.15g", "%s",
                 "%.15g", "%.15g", "%.15g", "%.15g", "%.15g"]
        )

        np.savetxt(outPath+"/ref.txt", arr, delimiter=",", **saveKwargs)
        np.savetxt(outPath+"/ref_test_delim.txt", arr, delimiter="|", **saveKwargs)
        return outPath+"/ref.txt", outPath+"/ref_test_delim.txt", arr

    @classmethod
    def setUpClass(cls):
        cls.outPath = tempfile.mkdtemp()
        cls.testCatPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data",
                                       "testHtmIndex.fits")
        # arbitrary, but reasonable, amount of proper motion (angle/year)
        # and direction of proper motion
        cls.properMotionAmt = 3.0*lsst.geom.arcseconds
        cls.properMotionDir = 45*lsst.geom.degrees
        cls.properMotionErr = 1e-3*lsst.geom.arcseconds
        cls.epoch = astropy.time.Time(58206.861330339219, scale="tai", format="mjd")
        ret = cls.make_skyCatalog(cls.outPath)
        cls.skyCatalogFile, cls.skyCatalogFileDelim, cls.skyCatalog = ret
        cls.testRas = [210., 14.5, 93., 180., 286., 0.]
        cls.testDecs = [-90., -51., -30.1, 0., 27.3, 62., 90.]
        cls.searchRadius = 3. * lsst.geom.degrees
        cls.compCats = {}  # dict of center coord: list of IDs of stars within cls.searchRadius of center
        cls.depth = 4  # gives a mean area of 20 deg^2 per pixel, roughly matching a 3 deg search radius

        config = IndexerRegistry['HTM'].ConfigClass()
        # Match on disk comparison file
        config.depth = cls.depth
        cls.indexer = IndexerRegistry['HTM'](config)
        for ra in cls.testRas:
            for dec in cls.testDecs:
                tupl = (ra, dec)
                cent = make_coord(*tupl)
                cls.compCats[tupl] = []
                for rec in cls.skyCatalog:
                    if make_coord(rec['ra_icrs'], rec['dec_icrs']).separation(cent) < cls.searchRadius:
                        cls.compCats[tupl].append(rec['id'])

        cls.testRepoPath = cls.outPath+"/test_repo"
        config = cls.makeConfig(withMagErr=True, withRaDecErr=True, withPm=True, withPmErr=True)
        # To match on disk test data
        config.dataset_config.indexer.active.depth = cls.depth
        config.id_name = 'id'
        config.pm_scale = 1000.0  # arcsec/yr --> mas/yr
        IngestIndexedReferenceTask.parseAndRun(args=[INPUT_DIR, "--output", cls.testRepoPath,
                                                     cls.skyCatalogFile], config=config)
        cls.defaultDatasetName = config.dataset_config.ref_dataset_name
        cls.testDatasetName = 'diff_ref_name'
        cls.testButler = dafPersist.Butler(cls.testRepoPath)
        os.symlink(os.path.join(cls.testRepoPath, 'ref_cats', cls.defaultDatasetName),
                   os.path.join(cls.testRepoPath, 'ref_cats', cls.testDatasetName))

    @classmethod
    def tearDownClass(cls):
        try:
            shutil.rmtree(cls.outPath)
        except Exception:
            print("WARNING: failed to remove temporary dir %r" % (cls.outPath,))
        del cls.outPath
        del cls.skyCatalogFile
        del cls.skyCatalogFileDelim
        del cls.skyCatalog
        del cls.testRas
        del cls.testDecs
        del cls.searchRadius
        del cls.compCats
        del cls.testButler

    def testSanity(self):
        """Sanity-check that compCats contains some entries with sources."""
        numWithSources = 0
        for idList in self.compCats.values():
            if len(idList) > 0:
                numWithSources += 1
        self.assertGreater(numWithSources, 0)

    def testAgainstPersisted(self):
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
        for kk in ex1:
            np.testing.assert_array_almost_equal(ex1[kk], ex2[kk], )

    @staticmethod
    def makeConfig(withMagErr=False, withRaDecErr=False, withPm=False, withPmErr=False,
                   withParallax=False, withParallaxErr=False):
        """Make a config for IngestIndexedReferenceTask

        This is primarily intended to simplify tests of config validation,
        so fields that are not validated are not set.
        However, it can calso be used to reduce boilerplate in other tests.
        """
        config = IngestIndexedReferenceTask.ConfigClass()
        config.pm_scale = 1000.0
        config.ra_name = 'ra_icrs'
        config.dec_name = 'dec_icrs'
        config.mag_column_list = ['a', 'b']

        if withMagErr:
            config.mag_err_column_map = {'a': 'a_err', 'b': 'b_err'}

        if withRaDecErr:
            config.ra_err_name = "ra_err"
            config.dec_err_name = "dec_err"

        if withPm:
            config.pm_ra_name = "pm_ra"
            config.pm_dec_name = "pm_dec"

        if withPmErr:
            config.pm_ra_err_name = "pm_ra_err"
            config.pm_dec_err_name = "pm_dec_err"

        if withParallax:
            config.parallax_name = "parallax"

        if withParallaxErr:
            config.parallax_err_name = "parallax_err"

        if withPm or withParallax:
            config.epoch_name = "unixtime"
            config.epoch_format = "unix"
            config.epoch_scale = "utc"

        return config

    def testValidateRaDecMag(self):
        config = self.makeConfig()
        config.validate()

        for name in ("ra_name", "dec_name", "mag_column_list"):
            with self.subTest(name=name):
                config = self.makeConfig()
                setattr(config, name, None)
                with self.assertRaises(ValueError):
                    config.validate()

    def testValidateRaDecErr(self):
        config = self.makeConfig(withRaDecErr=True)
        config.validate()

        for name in ("ra_err_name", "dec_err_name"):
            with self.subTest(name=name):
                config = self.makeConfig(withRaDecErr=True)
                setattr(config, name, None)
                with self.assertRaises(ValueError):
                    config.validate()

    def testValidateMagErr(self):
        config = self.makeConfig(withMagErr=True)
        config.validate()

        # test for missing names
        for name in config.mag_column_list:
            with self.subTest(name=name):
                config = self.makeConfig(withMagErr=True)
                del config.mag_err_column_map[name]
                with self.assertRaises(ValueError):
                    config.validate()

        # test for incorrect names
        for name in config.mag_column_list:
            with self.subTest(name=name):
                config = self.makeConfig(withMagErr=True)
                config.mag_err_column_map["badName"] = config.mag_err_column_map[name]
                del config.mag_err_column_map[name]
                with self.assertRaises(ValueError):
                    config.validate()

    def testValidatePm(self):
        basicNames = ["pm_ra_name", "pm_dec_name", "epoch_name", "epoch_format", "epoch_scale"]

        for withPmErr in (False, True):
            config = self.makeConfig(withPm=True, withPmErr=withPmErr)
            config.validate()
            del config

            if withPmErr:
                names = basicNames + ["pm_ra_err_name", "pm_dec_err_name"]
            else:
                names = basicNames
                for name in names:
                    with self.subTest(name=name, withPmErr=withPmErr):
                        config = self.makeConfig(withPm=True, withPmErr=withPmErr)
                        setattr(config, name, None)
                        with self.assertRaises(ValueError):
                            config.validate()

    def testValidateParallax(self):
        basicNames = ["parallax_name", "epoch_name", "epoch_format", "epoch_scale"]

        for withParallaxErr in (False, True):
            config = self.makeConfig(withParallax=True, withParallaxErr=withParallaxErr)
            config.validate()
            del config

            if withParallaxErr:
                names = basicNames + ["parallax_err_name"]
            else:
                names = basicNames
                for name in names:
                    with self.subTest(name=name, withParallaxErr=withParallaxErr):
                        config = self.makeConfig(withParallax=True, withParallaxErr=withParallaxErr)
                        setattr(config, name, None)
                        if name == "parallax_name" and not withParallaxErr:
                            # it is OK to omit parallax_name if no parallax_err_name
                            config.validate()
                        else:
                            with self.assertRaises(ValueError):
                                config.validate()

    def testIngest(self):
        """Test IngestIndexedReferenceTask."""
        # Test with multiple files and standard config
        config = self.makeConfig(withRaDecErr=True, withMagErr=True, withPm=True, withPmErr=True)
        IngestIndexedReferenceTask.parseAndRun(
            args=[INPUT_DIR, "--output", self.outPath+"/output_multifile",
                  self.skyCatalogFile, self.skyCatalogFile],
            config=config)
        # A newly-ingested refcat should be marked format_version=1.
        loader = LoadIndexedReferenceObjectsTask(butler=dafPersist.Butler(self.outPath+"/output_multifile"))
        self.assertEqual(loader.dataset_config.format_version, 1)

        # Test with config overrides
        config2 = self.makeConfig(withRaDecErr=True, withMagErr=True, withPm=True, withPmErr=True)
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
            'pm_dec_err', 'unixtime',
        ]
        config2.file_reader.delimiter = '|'
        # this also tests changing the delimiter
        IngestIndexedReferenceTask.parseAndRun(
            args=[INPUT_DIR, "--output", self.outPath+"/output_override",
                  self.skyCatalogFileDelim], config=config2)

        # This location is known to have objects
        cent = make_coord(93.0, -90.0)

        # Test if we can get back the catalog with a non-standard dataset name
        butler = dafPersist.Butler(self.outPath+"/output_override")
        loaderConfig = LoadIndexedReferenceObjectsConfig()
        loaderConfig.ref_dataset_name = "myrefcat"
        loader = LoadIndexedReferenceObjectsTask(butler=butler, config=loaderConfig)
        cat = loader.loadSkyCircle(cent, self.searchRadius, filterName='a').refCat
        self.assertTrue(len(cat) > 0)
        self.assertTrue(cat.isContiguous())

        # test that a catalog can be loaded even with a name not used for ingestion
        butler = dafPersist.Butler(self.testRepoPath)
        loaderConfig2 = LoadIndexedReferenceObjectsConfig()
        loaderConfig2.ref_dataset_name = self.testDatasetName
        loader = LoadIndexedReferenceObjectsTask(butler=butler, config=loaderConfig2)
        cat = loader.loadSkyCircle(cent, self.searchRadius, filterName='a').refCat
        self.assertTrue(len(cat) > 0)
        self.assertTrue(cat.isContiguous())

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
                for suffix in ("x", "y"):
                    self.assertTrue(np.all(np.isnan(lcat.refCat["centroid_%s" % (suffix,)])))
                self.assertFalse(np.any(lcat.refCat["hasCentroid"]))
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

    def testLoadVersion0(self):
        """Test reading a pre-written format_version=0 (Jy flux) catalog.
        It should be converted to have nJy fluxes.
        """
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data/version0')
        loader = LoadIndexedReferenceObjectsTask(butler=dafPersist.Butler(path))
        self.assertEqual(loader.dataset_config.format_version, 0)
        result = loader.loadSkyCircle(make_coord(10, 20), 5*lsst.geom.degrees, filterName='a')
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
        result = loader.loadSkyCircle(make_coord(10, 20), 5*lsst.geom.degrees, filterName='a')
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
