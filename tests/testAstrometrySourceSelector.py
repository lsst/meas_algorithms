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
from __future__ import division, absolute_import, print_function
from builtins import range

import unittest
import numpy as np

import lsst.afw.table as afwTable
from lsst.meas.algorithms import sourceSelector
import lsst.meas.base.tests
import lsst.utils.tests

# badFlags are flags that, if the bad flag is set the objecct should be
# imediately thrown out. This is as opposed to combinations of flags
# signifying good or bad data. Such flags should then be a part of the
# isGood test.
badFlags = ["base_PixelFlags_flag_edge",
            "base_PixelFlags_flag_interpolatedCenter",
            "base_PixelFlags_flag_saturated",
            "base_PixelFlags_flag_saturatedCenter",
            "base_PixelFlags_flag_crCenter",
            "base_PixelFlags_flag_bad",
            "base_PixelFlags_flag_interpolated",
            "slot_Centroid_flag",
            "slot_ApFlux_flag",
            ]


def add_good_source(src, num=0):
    """Insert a likely-good source into the catalog."""
    """
    num is added to various values to distinguish them in catalogs with multiple objects.
    """
    src.addNew()
    src['coord_ra'][-1] = 1. + num
    src['coord_dec'][-1] = 2. + num
    src['slot_Centroid_x'][-1] = 10. + num
    src['slot_Centroid_y'][-1] = 20. + num
    src['slot_ApFlux_flux'][-1] = 100. + num
    src['slot_ApFlux_fluxSigma'][-1] = 1.


class TestAstrometrySourceSelector(lsst.utils.tests.TestCase):

    def setUp(self):
        schema = lsst.meas.base.tests.TestDataset.makeMinimalSchema()
        schema.addField("slot_ApFlux_flux", type=np.float64)
        schema.addField("slot_ApFlux_fluxSigma", type=np.float64)
        for flag in badFlags:
            schema.addField(flag, type="Flag")

        self.src = afwTable.SourceCatalog(schema)
        self.sourceSelector = sourceSelector.sourceSelectorRegistry['astrometry']()

    def tearDown(self):
        del self.src
        del self.sourceSelector

    def testSelectSources_good(self):
        for i in range(5):
            add_good_source(self.src, i)
        result = self.sourceSelector.selectSources(self.src)
        # TODO: assertEqual doesn't work correctly on source catalogs.
        # self.assertEqual(result.sourceCat, self.src)
        for x in self.src['id']:
            self.assertIn(x, result.sourceCat['id'])

    def testSelectSources_bad(self):
        for i, flag in enumerate(badFlags):
            add_good_source(self.src, i)
            self.src[i].set(flag, True)
        result = self.sourceSelector.selectSources(self.src)
        for i, x in enumerate(self.src['id']):
            self.assertNotIn(x, result.sourceCat['id'], "should not have found %s" % badFlags[i])

    def testSelectSources_bad_centroid(self):
        add_good_source(self.src, 1)
        self.src[0].set('slot_Centroid_x', np.nan)
        result = self.sourceSelector.selectSources(self.src)
        self.assertNotIn(self.src['id'][0], result.sourceCat['id'])

    def testSelectSources_is_parent(self):
        add_good_source(self.src, 1)
        self.src[0].set('parent', 1)
        result = self.sourceSelector.selectSources(self.src)
        self.assertNotIn(self.src['id'][0], result.sourceCat['id'])

    def testSelectSources_has_children(self):
        add_good_source(self.src, 1)
        self.src[0].set('deblend_nChild', 1)
        result = self.sourceSelector.selectSources(self.src)
        self.assertNotIn(self.src['id'][0], result.sourceCat['id'])

    def testSelectSources_highSN_cut(self):
        add_good_source(self.src, 1)
        add_good_source(self.src, 2)
        self.src['slot_ApFlux_flux'][0] = 20.
        self.src['slot_ApFlux_flux'][1] = 1000.

        self.sourceSelector.config.minSnr = 100
        result = self.sourceSelector.selectSources(self.src)
        self.assertNotIn(self.src[0]['id'], result.sourceCat['id'])
        self.assertIn(self.src[1]['id'], result.sourceCat['id'])

    def testSelectSources_no_SN_cut(self):
        self.sourceSelector.config.minSnr = 0
        add_good_source(self.src, 1)
        self.src['slot_ApFlux_flux'][0] = 0
        result = self.sourceSelector.selectSources(self.src)
        self.assertIn(self.src[0]['id'], result.sourceCat['id'])

    def testSelectSources_non_contiguous(self):
        """Should raise Pex:RuntimeError if sourceSelector fails on non-contiguous catalogs."""
        for i in range(3):
            add_good_source(self.src, i)
        del self.src[1]  # take one out of the middle to make it non-contiguous.
        self.assertFalse(self.src.isContiguous(), "Catalog is contiguous: the test won't work.")

        result = self.sourceSelector.selectSources(self.src)
        # NOTE: have to use find() to search non-contiguous catalogs.
        for x in self.src:
            self.assertTrue(result.sourceCat.find(x['id']))


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
