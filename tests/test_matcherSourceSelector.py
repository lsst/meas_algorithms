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

import lsst.afw.table as afwTable
from lsst.meas.algorithms import sourceSelector
import lsst.meas.base.tests
import lsst.utils.tests

badFlags = [
    "slot_Centroid_flag",
    "slot_ApFlux_flag",
    "base_PixelFlags_flag_edge",
    "base_PixelFlags_flag_interpolatedCenter",
    "base_PixelFlags_flag_saturated",
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
    src['slot_ApFlux_instFlux'][-1] = 100. + num
    src['slot_ApFlux_instFluxErr'][-1] = 1.


class TestMatcherSourceSelector(lsst.utils.tests.TestCase):

    def setUp(self):
        schema = lsst.meas.base.tests.TestDataset.makeMinimalSchema()
        schema.addField("slot_ApFlux_instFlux", type=float)
        schema.addField("slot_ApFlux_instFluxErr", type=float)
        for flag in badFlags:
            schema.addField(flag, type="Flag")

        self.src = afwTable.SourceCatalog(schema)
        self.sourceSelector = sourceSelector.sourceSelectorRegistry['matcher']()

    def tearDown(self):
        del self.src
        del self.sourceSelector

    def testSelectSources_good(self):
        for i in range(5):
            add_good_source(self.src, i)
        result = self.sourceSelector.run(self.src)
        # TODO: assertEqual doesn't work correctly on source catalogs.
        # self.assertEqual(result.sourceCat, self.src)
        for x in self.src['id']:
            self.assertIn(x, result.sourceCat['id'])
        self.assertTrue(result.sourceCat.isContiguous())

    def testSelectSources_bad_centroid(self):
        add_good_source(self.src, 1)
        self.src[0].set('slot_Centroid_x', np.nan)
        result = self.sourceSelector.run(self.src)
        self.assertNotIn(self.src['id'][0], result.sourceCat['id'])
        self.assertTrue(result.sourceCat.isContiguous())

    def testSelectSources_is_parent(self):
        add_good_source(self.src, 1)
        self.src[0].set('parent', 1)
        result = self.sourceSelector.run(self.src)
        self.assertNotIn(self.src['id'][0], result.sourceCat['id'])
        self.assertTrue(result.sourceCat.isContiguous())

    def testSelectSources_highSN_cut(self):
        add_good_source(self.src, 1)
        add_good_source(self.src, 2)
        self.src['slot_ApFlux_instFlux'][0] = 20.
        self.src['slot_ApFlux_instFlux'][1] = 1000.

        self.sourceSelector.config.minSnr = 100
        result = self.sourceSelector.run(self.src)
        self.assertNotIn(self.src[0]['id'], result.sourceCat['id'])
        self.assertIn(self.src[1]['id'], result.sourceCat['id'])
        self.assertTrue(result.sourceCat.isContiguous())

    def testSelectSources_no_SN_cut(self):
        self.sourceSelector.config.minSnr = 0
        add_good_source(self.src, 1)
        self.src['slot_ApFlux_instFlux'][0] = 0
        result = self.sourceSelector.run(self.src)
        self.assertIn(self.src[0]['id'], result.sourceCat['id'])
        self.assertTrue(result.sourceCat.isContiguous())

    def testSelectSources_is_edge(self):
        add_good_source(self.src, 1)
        self.src[0].set('base_PixelFlags_flag_edge', 1)
        result = self.sourceSelector.run(self.src)
        self.assertNotIn(self.src['id'][0], result.sourceCat['id'])
        self.assertTrue(result.sourceCat.isContiguous())

    def testSelectSources_is_interpolated(self):
        add_good_source(self.src, 1)
        self.src[0].set('base_PixelFlags_flag_interpolatedCenter', 1)
        result = self.sourceSelector.run(self.src)
        self.assertNotIn(self.src['id'][0], result.sourceCat['id'])
        self.assertTrue(result.sourceCat.isContiguous())

    def testSelectSources_is_saturated(self):
        add_good_source(self.src, 1)
        self.src[0].set('base_PixelFlags_flag_saturated', 1)
        result = self.sourceSelector.run(self.src)
        self.assertNotIn(self.src['id'][0], result.sourceCat['id'])
        self.assertTrue(result.sourceCat.isContiguous())


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
