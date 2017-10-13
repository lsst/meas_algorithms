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

from lsst.meas.algorithms import ScienceSourceSelectorTask, ReferenceSourceSelectorTask, ColorLimit


class SourceSelectorTester(object):
    """Mixin for testing

    This provides a base class for doing tests common to the
    ScienceSourceSelectorTask and ReferenceSourceSelectorTask.
    """
    Task = None

    def setUp(self):
        schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        schema.addField("flux", float, "Flux value")
        schema.addField("flux_flag", "Flag", "Bad flux?")
        schema.addField("other_flux", float, "Flux value 2")
        schema.addField("other_flux_flag", "Flag", "Bad flux 2?")
        schema.addField("goodFlag", "Flag", "Flagged if good")
        schema.addField("badFlag", "Flag", "Flagged if bad")
        schema.addField("starGalaxy", float, "0=star, 1=galaxy")
        self.catalog = lsst.afw.table.SourceCatalog(schema)
        self.catalog.reserve(10)
        self.config = self.Task.ConfigClass()

    def tearDown(self):
        del self.catalog

    def check(self, expected):
        task = self.Task(config=self.config)
        results = task.selectSources(self.catalog)
        self.assertListEqual(results.selection.tolist(), expected)
        self.assertListEqual([src.getId() for src in results.sourceCat],
                             [src.getId() for src, ok in zip(self.catalog, expected) if ok])

    def testFlags(self):
        bad1 = self.catalog.addNew()
        bad1.set("goodFlag", False)
        bad1.set("badFlag", False)
        bad2 = self.catalog.addNew()
        bad2.set("goodFlag", True)
        bad2.set("badFlag", True)
        bad3 = self.catalog.addNew()
        bad3.set("goodFlag", False)
        bad3.set("badFlag", True)
        good = self.catalog.addNew()
        good.set("goodFlag", True)
        good.set("badFlag", False)
        self.catalog["flux"] = 1.0
        self.catalog["other_flux"] = 1.0
        self.config.flags.good = ["goodFlag"]
        self.config.flags.bad = ["badFlag"]
        self.check([False, False, False, True])


class ScienceSourceSelectorTaskTest(SourceSelectorTester, lsst.utils.tests.TestCase):
    Task = ScienceSourceSelectorTask

    def setUp(self):
        SourceSelectorTester.setUp(self)
        self.config.fluxLimit.fluxField = "flux"
        self.config.flags.bad = []
        self.config.doFluxLimit = True
        self.config.doFlags = True
        self.config.doUnresolved = False

    def testFluxLimit(self):
        tooBright = self.catalog.addNew()
        tooBright.set("flux", 1.0e10)
        tooBright.set("flux_flag", False)
        good = self.catalog.addNew()
        good.set("flux", 1000.0)
        good.set("flux_flag", False)
        bad = self.catalog.addNew()
        bad.set("flux", good.get("flux"))
        bad.set("flux_flag", True)
        tooFaint = self.catalog.addNew()
        tooFaint.set("flux", 1.0)
        tooFaint.set("flux_flag", False)
        self.config.fluxLimit.minimum = 10.0
        self.config.fluxLimit.maximum = 1.0e6
        self.config.fluxLimit.fluxField = "flux"
        self.check([False, True, False, False])

        # Works with no maximum set?
        self.config.fluxLimit.maximum = None
        self.check([True, True, False, False])

        # Works with no minimum set?
        self.config.fluxLimit.minimum = None
        self.check([True, True, False, True])

    def testUnresolved(self):
        num = 5
        for _ in range(num):
            self.catalog.addNew()
        self.catalog["flux"] = 1.0
        starGalaxy = np.linspace(0.0, 1.0, num, False)
        self.catalog["starGalaxy"] = starGalaxy
        self.config.doUnresolved = True
        self.config.unresolved.name = "starGalaxy"
        minimum, maximum = 0.3, 0.7
        self.config.unresolved.minimum = minimum
        self.config.unresolved.maximum = maximum
        self.check(((starGalaxy > minimum) & (starGalaxy < maximum)).tolist())

        # Works with no minimum set?
        self.config.unresolved.minimum = None
        self.check((starGalaxy < maximum).tolist())

        # Works with no maximum set?
        self.config.unresolved.minimum = minimum
        self.config.unresolved.maximum = None
        self.check((starGalaxy > minimum).tolist())


class ReferenceSourceSelectorTaskTest(SourceSelectorTester, lsst.utils.tests.TestCase):
    Task = ReferenceSourceSelectorTask

    def setUp(self):
        SourceSelectorTester.setUp(self)
        self.config.magLimit.fluxField = "flux"
        self.config.doMagLimit = True
        self.config.doFlags = True

    def testMagnitudeLimit(self):
        tooBright = self.catalog.addNew()
        tooBright.set("flux", 1.0e10)
        tooBright.set("flux_flag", False)
        good = self.catalog.addNew()
        good.set("flux", 1000.0)
        good.set("flux_flag", False)
        bad = self.catalog.addNew()
        bad.set("flux", good.get("flux"))
        bad.set("flux_flag", True)
        tooFaint = self.catalog.addNew()
        tooFaint.set("flux", 1.0)
        tooFaint.set("flux_flag", False)
        # Note: magnitudes are backwards, so the minimum flux is the maximum magnitude
        self.config.magLimit.minimum = lsst.afw.image.abMagFromFlux(1.0e6)
        self.config.magLimit.maximum = lsst.afw.image.abMagFromFlux(10.0)
        self.config.magLimit.fluxField = "flux"
        self.check([False, True, False, False])

        # Works with no minimum set?
        self.config.magLimit.minimum = None
        self.check([True, True, False, False])

        # Works with no maximum set?
        self.config.magLimit.maximum = None
        self.check([True, True, False, True])

    def testColorLimits(self):
        num = 10
        for _ in range(num):
            self.catalog.addNew()
        color = np.linspace(-0.5, 0.5, num, True)
        flux = 1000.0
        # Definition: color = mag(flux) - mag(otherFlux)
        otherFlux = lsst.afw.image.fluxFromABMag(lsst.afw.image.abMagFromFlux(flux) - color)
        self.catalog["flux"] = flux
        self.catalog["other_flux"] = otherFlux
        minimum, maximum = -0.1, 0.2
        self.config.colorLimits = {"test": ColorLimit(primary="flux", secondary="other_flux",
                                                      minimum=minimum, maximum=maximum)}
        self.check(((color > minimum) & (color < maximum)).tolist())

        # Works with no minimum set?
        self.config.colorLimits["test"].minimum = None
        self.check((color < maximum).tolist())

        # Works with no maximum set?
        self.config.colorLimits["test"].maximum = None
        self.config.colorLimits["test"].minimum = minimum
        self.check((color > minimum).tolist())

        # Multiple limits
        self.config.colorLimits = {"test": ColorLimit(primary="flux", secondary="other_flux",
                                                      minimum=minimum),
                                   "other": ColorLimit(primary="flux", secondary="other_flux",
                                                       maximum=maximum)}
        assert maximum > minimum  # To be non-mutually-exclusive
        self.check(((color > minimum) & (color < maximum)).tolist())

        # Multiple mutually-exclusive limits
        self.config.colorLimits["test"] = ColorLimit(primary="flux", secondary="other_flux", maximum=-0.1)
        self.config.colorLimits["other"] = ColorLimit(primary="flux", secondary="other_flux", minimum=0.1)
        self.check([False]*num)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    import sys
    setup_module(sys.modules[__name__])
    unittest.main()
