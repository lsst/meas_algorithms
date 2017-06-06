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

import lsst.afw.table as afwTable
from lsst.meas.algorithms import sourceSelector
import lsst.meas.base.tests
from lsst.pex.exceptions import RuntimeError
import lsst.utils.tests


def add_good_source(src, num=0):
    """Insert a likely-good source into the catalog. Num is added to various
    values to distinguish them in catalogs with multiple objects.
    """
    src.addNew()
    src['coord_ra'][-1] = 1. + num
    src['coord_dec'][-1] = 2. + num
    src['slot_Centroid_x'][-1] = 10. + num
    src['slot_Centroid_y'][-1] = 20. + num
    src['slot_ApFlux_flux'][-1] = 100. + num
    src['slot_ApFlux_fluxSigma'][-1] = 1.
    src[-1].set("calib_psfUsed", True)


class TestFlaggedSourceSelector(lsst.utils.tests.TestCase):

    def setUp(self):
        schema = lsst.meas.base.tests.TestDataset.makeMinimalSchema()
        schema.addField("slot_ApFlux_flux", type=float)
        schema.addField("slot_ApFlux_fluxSigma", type=float)
        schema.addField("calib_psfUsed", type="Flag")
        schema.addField("is_selected", type="Flag")

        self.selected_key = "is_selected"

        self.src = afwTable.SourceCatalog(schema)
        self.sourceSelector = \
            sourceSelector.sourceSelectorRegistry['flagged'](schema=schema)

    def tearDown(self):
        del self.src
        del self.sourceSelector

    def testSelectSources_good(self):
        """Insert sources that pass our criteria and test that they indeed do
        so.
        """
        for i in range(5):
            add_good_source(self.src, i)
        result = self.sourceSelector.run(self.src)
        # TODO: assertEqual doesn't work correctly on source catalogs.
        # self.assertEqual(result.sourceCat, self.src)
        for x in self.src['id']:
            self.assertIn(x, result.source_cat['id'])

    def testSelectSources_selected_field(self):
        """Test the behavior of source_selected_field in run.

        This test asserts that the field will specified in
        source_selected_field is properly set by the source selector. We test
        this both for sources that fail and pass our cuts.
        """
        for i in range(5):
            add_good_source(self.src, i)
        self.src[0].set("calib_psfUsed", False)
        result = self.sourceSelector.run(
            self.src, source_selected_field=self.selected_key)
        for src_idx in range(5):
            if src_idx == 0:
                self.assertFalse(self.src[src_idx].get("is_selected"))
            else:
                self.assertTrue(self.src[src_idx].get("is_selected"))

    def testSelectSources_bad(self):
        """Add a source that fails the source selector test and check
        that our output array is indeed empty.
        """
        add_good_source(self.src, 1)
        self.src[0].set('calib_psfUsed', False)
        result = self.sourceSelector.run(self.src)
        self.assertNotIn(self.src['id'][0], result.source_cat['id'])

    def testSelectSources_non_contiguous(self):
        """Should raise Pex:RuntimeError if sourceSelector fails on
        non-contiguous catalogs.
        """
        for i in range(3):
            add_good_source(self.src, i)
            self.src[i].set("calib_psfUsed", True)
        # take one out of the middle to make it non-contiguous.
        del self.src[1]
        self.assertFalse(self.src.isContiguous(),
                         "Catalog is contiguous: the test won't work.")

        with self.assertRaises(RuntimeError):
            result = self.sourceSelector.run(self.src)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
