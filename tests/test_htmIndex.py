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
from builtins import zip

import math
import os
import tempfile
import shutil
import unittest
import string
from collections import Counter

import numpy as np

import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord
import lsst.daf.persistence as dafPersist
from lsst.meas.algorithms import (IngestIndexedReferenceTask, LoadIndexedReferenceObjectsTask,
                                  LoadIndexedReferenceObjectsConfig, getRefFluxField)
from lsst.meas.algorithms import IndexerRegistry
import lsst.utils

obs_test_dir = lsst.utils.getPackageDir('obs_test')
input_dir = os.path.join(obs_test_dir, "data", "input")


def make_coord(ra, dec):
    """Make an ICRS coord given its RA, Dec in degrees."""
    return afwCoord.IcrsCoord(afwGeom.Angle(ra, afwGeom.degrees), afwGeom.Angle(dec, afwGeom.degrees))


def makeWcs(ctr_coord, ctr_pix=afwGeom.Point2D(2036., 2000.),
            pixel_scale=2*afwGeom.arcseconds, pos_angle=afwGeom.Angle(0.0)):
    """Make a simple TAN WCS

    @param[in] ctr_coord  sky coordinate at ctr_pix
    @param[in] ctr_pix  center pixel; an lsst.afw.geom.Point2D; default matches LSST
    @param[in] pixel_scale  desired scale, as sky/pixel; an lsst.afw.geom.Angle; default matches LSST
    @param[in] pos_angle  orientation of CCD w.r.t. ctr_coord, an lsst.afw.geom.Angle
    """
    pos_angleRad = pos_angle.asRadians()
    pixel_scaleDeg = pixel_scale.asDegrees()
    cdMat = np.array([[math.cos(pos_angleRad), math.sin(pos_angleRad)],
                      [-math.sin(pos_angleRad), math.cos(pos_angleRad)]], dtype=float) * pixel_scaleDeg
    return lsst.afw.image.makeWcs(ctr_coord, ctr_pix, cdMat[0, 0], cdMat[0, 1], cdMat[1, 0], cdMat[1, 1])


class HtmIndexTestCase(lsst.utils.tests.TestCase):

    @staticmethod
    def make_sky_catalog(out_path, size=1000):
        np.random.seed(123)
        ident = np.arange(1, size+1, dtype=int)
        ra = np.random.random(size)*360.
        dec = np.degrees(np.arccos(2.*np.random.random(size) - 1.))
        dec -= 90.
        a_mag = 16. + np.random.random(size)*4.
        a_mag_err = 0.01 + np.random.random(size)*0.2
        b_mag = 17. + np.random.random(size)*5.
        b_mag_err = 0.02 + np.random.random(size)*0.3
        is_photometric = np.random.randint(2, size=size)
        is_resolved = np.random.randint(2, size=size)
        is_variable = np.random.randint(2, size=size)
        extra_col1 = np.random.normal(size=size)
        extra_col2 = np.random.normal(1000., 100., size=size)

        def get_word(word_len):
            return "".join(np.random.choice([s for s in string.ascii_letters], word_len))
        extra_col3 = np.array([get_word(num) for num in np.random.randint(11, size=size)])

        dtype = np.dtype([('id', float), ('ra_icrs', float), ('dec_icrs', float), ('a', float),
                          ('a_err', float), ('b', float), ('b_err', float), ('is_phot', int),
                          ('is_res', int), ('is_var', int), ('val1', float), ('val2', float),
                          ('val3', '|S11')])

        arr = np.array(list(zip(ident, ra, dec, a_mag, a_mag_err, b_mag, b_mag_err, is_photometric,
                                is_resolved, is_variable, extra_col1, extra_col2, extra_col3)), dtype=dtype)
        np.savetxt(out_path+"/ref.txt", arr, delimiter=",",
                   header="id,ra_icrs,dec_icrs,a,a_err,b,b_err,is_phot,is_res,is_var,val1,val2,val3",
                   fmt=["%i", "%.6g", "%.6g", "%.4g", "%.4g", "%.4g", "%.4g", "%i",
                        "%i", "%i", "%.2g", "%.2g", "%s"])
        np.savetxt(out_path+"/ref_test_delim.txt", arr, delimiter="|",
                   header="id,ra_icrs,dec_icrs,a,a_err,b,b_err,is_phot,is_res,is_var,val1,val2,val3",
                   fmt=["%i", "%.6g", "%.6g", "%.4g", "%.4g", "%.4g", "%.4g", "%i",
                        "%i", "%i", "%.2g", "%.2g", "%s"])
        return out_path+"/ref.txt", out_path+"/ref_test_delim.txt", arr

    @classmethod
    def setUpClass(cls):
        cls.out_path = tempfile.mkdtemp()
        test_cat_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data", "testHtmIndex.fits")
        cls.test_cat = afwTable.SourceCatalog.readFits(test_cat_path)
        ret = cls.make_sky_catalog(cls.out_path)
        cls.sky_catalog_file, cls.sky_catalog_file_delim, cls.sky_catalog = ret
        cls.test_ras = [210., 14.5, 93., 180., 286., 0.]
        cls.test_decs = [-90., -51., -30.1, 0., 27.3, 62., 90.]
        cls.search_radius = 3. * afwGeom.degrees
        cls.comp_cats = {}  # dict of center coord: list of IDs of stars within cls.search_radius of center
        config = IndexerRegistry['HTM'].ConfigClass()
        # Match on disk comparison file
        config.depth = 8
        cls.indexer = IndexerRegistry['HTM'](config)
        for ra in cls.test_ras:
            for dec in cls.test_decs:
                tupl = (ra, dec)
                cent = make_coord(*tupl)
                cls.comp_cats[tupl] = []
                for rec in cls.sky_catalog:
                    if make_coord(rec['ra_icrs'], rec['dec_icrs']).angularSeparation(cent) \
                       < cls.search_radius:
                        cls.comp_cats[tupl].append(rec['id'])

        cls.test_repo_path = cls.out_path+"/test_repo"
        config = IngestIndexedReferenceTask.ConfigClass()
        # To match on disk test data
        config.dataset_config.indexer.active.depth = 8
        config.ra_name = 'ra_icrs'
        config.dec_name = 'dec_icrs'
        config.mag_column_list = ['a', 'b']
        config.id_name = 'id'
        config.mag_err_column_map = {'a': 'a_err', 'b': 'b_err'}
        IngestIndexedReferenceTask.parseAndRun(args=[input_dir, "--output", cls.test_repo_path,
                                                     cls.sky_catalog_file], config=config)
        cls.default_dataset_name = config.dataset_config.ref_dataset_name
        cls.test_dataset_name = 'diff_ref_name'
        cls.test_butler = dafPersist.Butler(cls.test_repo_path)
        os.symlink(os.path.join(cls.test_repo_path, 'ref_cats', cls.default_dataset_name),
                   os.path.join(cls.test_repo_path, 'ref_cats', cls.test_dataset_name))

    @classmethod
    def tearDownClass(cls):
        try:
            shutil.rmtree(cls.out_path)
        except Exception:
            print("WARNING: failed to remove temporary dir %r" % (cls.out_path,))
        del cls.out_path
        del cls.sky_catalog_file
        del cls.sky_catalog_file_delim
        del cls.sky_catalog
        del cls.test_ras
        del cls.test_decs
        del cls.search_radius
        del cls.comp_cats
        del cls.test_butler
        del cls.test_cat

    def testSanity(self):
        """Sanity-check that comp_cats contains some entries with sources."""
        numWithSources = 0
        for idList in self.comp_cats.values():
            if len(idList) > 0:
                numWithSources += 1
        self.assertGreater(numWithSources, 0)

    def testAgainstPersisted(self):
        pix_id = 671901
        dataset_name = IngestIndexedReferenceTask.ConfigClass().dataset_config.ref_dataset_name
        data_id = self.indexer.make_data_id(pix_id, dataset_name)
        self.assertTrue(self.test_butler.datasetExists('ref_cat', data_id))
        ref_cat = self.test_butler.get('ref_cat', data_id)
        ex1 = ref_cat.extract('*')
        ex2 = self.test_cat.extract('*')
        # compare sets as the order may be different
        self.assertDictEqual(ex1, ex2)

    def testIngest(self):
        """Test IngestIndexedReferenceTask."""
        default_config = IngestIndexedReferenceTask.ConfigClass()
        # test ingest with default config
        # This should raise since I haven't specified the ra/dec/mag columns.
        with self.assertRaises(ValueError):
            IngestIndexedReferenceTask.parseAndRun(
                args=[input_dir, "--output", self.out_path+"/output", self.sky_catalog_file],
                config=default_config)
        # test with ~minimum config.  Mag errors are not technically necessary, but might as well test here
        default_config.ra_name = 'ra_icrs'
        default_config.dec_name = 'dec_icrs'
        default_config.mag_column_list = ['a', 'b']
        default_config.mag_err_column_map = {'a': 'a_err'}
        # should raise since all columns need an error column if any do
        with self.assertRaises(ValueError):
            IngestIndexedReferenceTask.parseAndRun(
                args=[input_dir, "--output", self.out_path+"/output", self.sky_catalog_file],
                config=default_config)
        # test with multiple files and correct config
        default_config.mag_err_column_map = {'a': 'a_err', 'b': 'b_err'}
        IngestIndexedReferenceTask.parseAndRun(
            args=[input_dir, "--output", self.out_path+"/output_multifile",
                  self.sky_catalog_file, self.sky_catalog_file],
            config=default_config)
        # test with config overrides
        default_config = IngestIndexedReferenceTask.ConfigClass()
        default_config.ra_name = 'ra'
        default_config.dec_name = 'dec'
        default_config.mag_column_list = ['a', 'b']
        default_config.mag_err_column_map = {'a': 'a_err', 'b': 'b_err'}
        default_config.dataset_config.ref_dataset_name = 'myrefcat'
        default_config.dataset_config.indexer.active.depth = 10
        default_config.is_photometric_name = 'is_phot'
        default_config.is_resolved_name = 'is_res'
        default_config.is_variable_name = 'is_var'
        default_config.id_name = 'id'
        default_config.extra_col_names = ['val1', 'val2', 'val3']
        default_config.file_reader.header_lines = 1
        default_config.file_reader.colnames = ['id', 'ra', 'dec', 'a', 'a_err', 'b', 'b_err', 'is_phot',
                                               'is_res', 'is_var', 'val1', 'val2', 'val3']
        default_config.file_reader.delimiter = '|'
        # this also tests changing the delimiter
        IngestIndexedReferenceTask.parseAndRun(
            args=[input_dir, "--output", self.out_path+"/output_override",
                  self.sky_catalog_file_delim], config=default_config)

        # This location is known to have objects
        cent = make_coord(93.0, -90.0)

        # Test if we can get back the catalog with a non-standard dataset name
        butler = dafPersist.Butler(self.out_path+"/output_override")
        config = LoadIndexedReferenceObjectsConfig()
        config.ref_dataset_name = "myrefcat"
        loader = LoadIndexedReferenceObjectsTask(butler=butler, config=config)
        cat = loader.loadSkyCircle(cent, self.search_radius, filterName='a')
        self.assertTrue(len(cat) > 0)

        # test that a catalog can be loaded even with a name not used for ingestion
        butler = dafPersist.Butler(self.test_repo_path)
        config = LoadIndexedReferenceObjectsConfig()
        config.ref_dataset_name = self.test_dataset_name
        loader = LoadIndexedReferenceObjectsTask(butler=butler, config=config)
        cat = loader.loadSkyCircle(cent, self.search_radius, filterName='a')
        self.assertTrue(len(cat) > 0)

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
        loader = LoadIndexedReferenceObjectsTask(butler=self.test_butler)
        for tupl, idList in self.comp_cats.items():
            cent = make_coord(*tupl)
            lcat = loader.loadSkyCircle(cent, self.search_radius, filterName='a')
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
        loader = LoadIndexedReferenceObjectsTask(butler=self.test_butler)
        numFound = 0
        for tupl, idList in self.comp_cats.items():
            cent = make_coord(*tupl)
            bbox = afwGeom.Box2I(afwGeom.Point2I(30, -5), afwGeom.Extent2I(1000, 1004))  # arbitrary
            ctr_pix = afwGeom.Box2D(bbox).getCenter()
            # catalog is sparse, so set pixel scale such that bbox encloses region
            # used to generate comp_cats
            pixel_scale = 2*self.search_radius/max(bbox.getHeight(), bbox.getWidth())
            wcs = makeWcs(ctr_coord=cent, ctr_pix=ctr_pix, pixel_scale=pixel_scale)
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
        loader = LoadIndexedReferenceObjectsTask(butler=self.test_butler, config=config)
        for tupl, idList in self.comp_cats.items():
            cent = make_coord(*tupl)
            lcat = loader.loadSkyCircle(cent, self.search_radius)
            self.assertEqual(lcat.fluxField, "camFlux")
            if len(idList) > 0:
                defFluxFieldName = getRefFluxField(lcat.refCat.schema, None)
                self.assertTrue(defFluxFieldName in lcat.refCat.schema)
                aprimeFluxFieldName = getRefFluxField(lcat.refCat.schema, "aprime")
                self.assertTrue(aprimeFluxFieldName in lcat.refCat.schema)
                break  # just need one test


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
