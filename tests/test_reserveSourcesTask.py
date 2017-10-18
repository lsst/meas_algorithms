#
# LSST Data Management System
#
# Copyright 2008-2017  AURA/LSST.
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

import unittest
import numpy as np

import lsst.afw.table
import lsst.meas.algorithms
import lsst.utils.tests

from lsst.pipe.base import Struct


class ReserveSourcesTaskTest(lsst.utils.tests.TestCase):
    """TestCase for the ReserveSourcesTask"""
    def setUp(self):
        self.num = 100  # Number of sources
        self.longMessage = True

    def construct(self, name, fraction):
        """Construct the test environment

        This isn't called 'setUp' because we want to vary the `fraction`.

        Parameters
        ----------
        name : `str`
            Name of column for flagging reserved sources (without "_reserved").
        fraction : `float`
            Fraction of sources to reserve.

        Return struct elements
        ----------------------
        catalog : `lsst.afw.table.SourceCatalog`
            Catalog of sources.
        task : `lsst.meas.algorithms.ReserveSourcesTask`
            Task to do the reservations.
        key : `lsst.afw.table.Key`
            Key to the flag column.
        """
        schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        config = lsst.meas.algorithms.ReserveSourcesConfig()
        config.fraction = fraction
        task = lsst.meas.algorithms.ReserveSourcesTask(columnName=name, schema=schema,
                                                       doc="Documentation is good", config=config)
        key = schema[name + "_reserved"].asKey()
        catalog = lsst.afw.table.SourceCatalog(schema)
        catalog.reserve(self.num)
        for _ in range(self.num):
            catalog.addNew()
        return Struct(catalog=catalog, task=task, key=key)

    def check(self, sources, task, key, fraction):
        """Check that source reservation is working

        We test that the source reservation works, that it depends
        on the RNG seed and that things behave as we expect when there's
        a prior selection.

        Parameters
        ----------
        catalog : `lsst.afw.table.Catalog` or `list` of `lsst.afw.table.Record`
            List of sources.
        task : `lsst.meas.algorithms.ReserveSourcesTask`
            Task to do the reservations.
        key : `lsst.afw.table.Key`
            Key to the flag column.
        fraction : `float`
            Fraction of sources to reserve.
        """
        message = "fraction=%s, key=%s" % (fraction,key)
        numExpect = int(fraction*len(sources))

        # No prior
        results1 = task.run(sources, expId=1)
        flagged1 = np.array([ss.get(key) for ss in sources])
        self.assertEqual(flagged1.sum(), numExpect, message)
        np.testing.assert_array_equal(results1.reserved, flagged1, message)
        np.testing.assert_array_equal(results1.use, np.logical_not(flagged1), message)

        # Second run with different seed; clear the flag first
        for ss in sources:
            ss.set(key, False)
        results2 = task.run(sources, expId=2)
        flagged2 = np.array([ss.get(key) for ss in sources])
        self.assertEqual(flagged1.sum(), numExpect, message)
        np.testing.assert_array_equal(results2.reserved, flagged2, message)
        np.testing.assert_array_equal(results2.use, np.logical_not(flagged2), message)
        if numExpect > 0:
            self.assertFalse(np.all(flagged1 == flagged2),
                             "Pretty unlikely since different seeds\n" + message)

        # Run with prior selection; clear the flag first
        for ss in sources:
            ss.set(key, False)
        prior = np.arange(0, self.num, 1, dtype=int) % 2 == 0
        results3 = task.run(sources, prior, expId=1)
        flagged3 = np.array([ss.get(key) for ss in sources])
        self.assertEqual(flagged3.sum(), fraction*prior.sum(), message)
        if numExpect > 0:
            self.assertFalse(np.all(flagged1 == flagged3),
                             "Flags should change, despite same see\n" + message)
        np.testing.assert_array_equal(results3.reserved, flagged3, message)
        self.assertEqual((results3.use & flagged3).sum(), 0,
                         "No sources should be both reserved and used\n" + message)  #
        self.assertEqual(results3.use.sum(), int((1.0 - fraction)*prior.sum()), message)
        self.assertEqual(results3.reserved.sum(), int(fraction*prior.sum()), message)
        np.testing.assert_array_equal(results3.use, prior & ~flagged3,
                                      "Actual definition of 'use'\n" + message)

    def testCatalog(self):
        """Test source reservation with a Catalog

        We test multiple reservation fractions.
        """
        for fraction in (0.0, 0.1, 0.2, 0.5):
            data = self.construct("this_table_is", fraction)
            self.check(data.catalog, data.task, data.key, fraction)

    def testSources(self):
        """Test source reservation with a list of sources"""
        fraction = 0.2
        data = self.construct("this_table_is", fraction)
        sources = [ss for ss in data.catalog]
        self.check(sources, data.task, data.key, fraction)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    import sys
    setup_module(sys.modules[__name__])
    unittest.main()
