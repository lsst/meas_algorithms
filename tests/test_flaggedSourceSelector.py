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
import unittest

import lsst.geom
import lsst.afw.table as afwTable
from lsst.meas.algorithms import sourceSelector
import lsst.meas.base.tests
import lsst.utils.tests


def addGoodSource(src, num=0):
    """Insert a likely-good source into the catalog. Num is added to various
    values to distinguish them in catalogs with multiple objects.
    """
    record = src.addNew()
    record['coord_ra'] = (1. + num) * lsst.geom.degrees
    record['coord_dec'] = (2. + num) * lsst.geom.degrees
    record['slot_Centroid_x'] = 10. + num
    record['slot_Centroid_y'] = 20. + num
    record['slot_ApFlux_flux'] = 100. + num
    record['slot_ApFlux_fluxErr'] = 1.
    record.set("calib_psfUsed", True)


class TestFlaggedSourceSelector(lsst.utils.tests.TestCase):

    def setUp(self):
        schema = lsst.meas.base.tests.TestDataset.makeMinimalSchema()
        schema.addField("slot_ApFlux_flux", type=float)
        schema.addField("slot_ApFlux_fluxErr", type=float)
        schema.addField("calib_psfUsed", type="Flag")

        self.src = afwTable.SourceCatalog(schema)
        self.sourceSelector = sourceSelector.sourceSelectorRegistry['flagged']()

    def tearDown(self):
        del self.src
        del self.sourceSelector

    def testSelectSourcesGood(self):
        """Insert good sources and check that they were selected."""
        for i in range(5):
            addGoodSource(self.src, i)
        result = self.sourceSelector.run(self.src)
        for x in self.src['id']:
            self.assertIn(x, result.sourceCat['id'])

    def testSelectSourcesBad(self):
        """Add a source that fails the source selector test and check
        that our output array is indeed empty.
        """
        addGoodSource(self.src, 1)
        self.src[0].set('calib_psfUsed', False)
        result = self.sourceSelector.run(self.src)
        self.assertNotIn(self.src['id'][0], result.sourceCat['id'])


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
