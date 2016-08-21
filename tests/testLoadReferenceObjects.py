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

import itertools
import unittest

import lsst.afw.table as afwTable
from lsst.meas.algorithms import LoadReferenceObjectsTask, getRefFluxField, getRefFluxKeys
import lsst.utils.tests


class TrivialLoader(LoadReferenceObjectsTask):
    """Minimal subclass of LoadReferenceObjectsTask to allow instantiation
    """

    def loadSkyCircle(self, ctrCoord, radius, filterName):
        pass


class TestLoadReferenceObjects(lsst.utils.tests.TestCase):
    """Test case for LoadReferenceObjectsTask abstract base class

    Only methods with concrete implementations are tested (hence not loadSkyCircle)
    """

    def testMakeMinimalSchema(self):
        """Make a schema and check it."""
        for filterNameList in (["r"], ["foo", "_bar"]):
            for addFluxSigma, addIsPhotometric, addIsResolved, addIsVariable in itertools.product(
                    (False, True), (False, True), (False, True), (False, True)):
                refSchema = LoadReferenceObjectsTask.makeMinimalSchema(
                    filterNameList=filterNameList,
                    addFluxSigma=addFluxSigma,
                    addIsPhotometric=addIsPhotometric,
                    addIsResolved=addIsResolved,
                    addIsVariable=addIsVariable,
                )
                self.assertEqual("resolved" in refSchema, addIsResolved)
                self.assertEqual("variable" in refSchema, addIsVariable)
                self.assertEqual("photometric" in refSchema, addIsPhotometric)
                for filterName in filterNameList:
                    fluxField = filterName + "_flux"
                    self.assertIn(fluxField, refSchema)
                    self.assertNotIn("x" + fluxField, refSchema)
                    fluxSigmaField = fluxField + "Sigma"
                    self.assertEqual(fluxSigmaField in refSchema, addFluxSigma)
                    self.assertEqual(getRefFluxField(refSchema, filterName), filterName + "_flux")

    def testFilterAliasMap(self):
        """Make a schema with filter aliases."""
        for defaultFilter in ("", "r", "camr"):
            for filterMap in ({}, {"camr": "r"}):
                for addFluxSigma in (False, True):
                    config = TrivialLoader.ConfigClass()
                    config.defaultFilter = defaultFilter
                    config.filterMap = filterMap
                    loader = TrivialLoader(config=config)
                    refSchema = TrivialLoader.makeMinimalSchema(filterNameList="r", addFluxSigma=addFluxSigma)
                    try:
                        loader._addFluxAliases(refSchema)
                        self.assertNotEqual(defaultFilter, "camr")
                    except Exception:
                        # only reference filters are allowed as default filters
                        self.assertEqual(defaultFilter, "camr")
                        continue

                    self.assertIn("r_flux", refSchema)
                    self.assertEqual("r_fluxSigma" in refSchema, addFluxSigma)

                    # camera filters aliases are named <filter>_camFlux
                    if "camr" in filterMap:
                        self.assertEqual(getRefFluxField(refSchema, "camr"), "camr_camFlux")
                    else:
                        with self.assertRaises(RuntimeError):
                            getRefFluxField(refSchema, "camr")

                    # if a non-empty default filter is specified then camFlux
                    # and camFluxSigma (if addFluxSigma) should be present
                    hasDefault = bool(defaultFilter)
                    self.assertEqual("camFlux" in refSchema, hasDefault)
                    self.assertEqual("camFluxSigma" in refSchema, hasDefault and addFluxSigma)

                    refCat = afwTable.SimpleCatalog(refSchema)
                    refObj = refCat.addNew()
                    refObj["r_flux"] = 1.23
                    self.assertAlmostEqual(refCat[0].get(getRefFluxField(refSchema, "r")), 1.23)
                    if "camr" in filterMap:
                        self.assertAlmostEqual(refCat[0].get(getRefFluxField(refSchema, "camr")), 1.23)
                    if hasDefault:
                        self.assertEqual(getRefFluxField(refSchema, ""), "camFlux")
                        self.assertAlmostEqual(refCat[0].get(getRefFluxField(refSchema, "")), 1.23)
                    if addFluxSigma:
                        refObj["r_fluxSigma"] = 0.111
                        if "camr" in filterMap:
                            self.assertEqual(refCat[0].get("camr_camFluxSigma"), 0.111)
                    fluxKey, fluxSigmaKey = getRefFluxKeys(refSchema, "r")
                    self.assertEqual(refCat[0].get(fluxKey), 1.23)
                    if addFluxSigma:
                        self.assertEqual(refCat[0].get(fluxSigmaKey), 0.111)
                    else:
                        self.assertEqual(fluxSigmaKey, None)
                    if "camr" in filterMap:
                        fluxKey, fluxSigmaKey = getRefFluxKeys(refSchema, "camr")
                        if addFluxSigma:
                            self.assertEqual(refCat[0].get(fluxSigmaKey), 0.111)
                        else:
                            self.assertEqual(fluxSigmaKey, None)
                    else:
                        with self.assertRaises(RuntimeError):
                            getRefFluxKeys(refSchema, "camr")


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
