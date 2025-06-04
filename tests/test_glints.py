# This file is part of meas_algorithms.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import unittest

import numpy as np

from lsst.meas.algorithms import findGlints
import lsst.geom
from lsst.meas.base.tests import TestDataset
import lsst.utils.tests


class TestFindGlints(lsst.utils.tests.TestCase):
    def setUp(self):
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(5, 4), lsst.geom.Point2I(1005, 1084))
        dataset = TestDataset(bbox)
        x0 = 200
        y0 = 300
        scale = 50
        rng = np.random.default_rng()
        # the glint trail
        for x in range(8):
            dataset.addSource(instFlux=10000,
                              centroid=lsst.geom.Point2D(x*scale + x0 + rng.random(),
                                                         x*scale + y0 + rng.random()))
        # a source that shouldn't appear in the found trail
        dataset.addSource(instFlux=10000, centroid=lsst.geom.Point2D(500, 200))

        schema = dataset.makeMinimalSchema()
        self.exposure, self.catalog = dataset.realize(10.0, schema=schema)
        # So that the catalog ids don't look like simple indices when debugging.
        self.catalog["id"] += 100

    def test_simple(self):
        import lsst.afw.display
        display = lsst.afw.display.Display()
        display.frame = 1
        display.image(self.exposure, title="something")
        display.centroids(self.catalog, size=20, ctype="red", symbol="x")
        # import ipdb; ipdb.set_trace();

        config = findGlints.FindGlintsTask.ConfigClass()
        config.radius = 300
        task = findGlints.FindGlintsTask(config=config)
        task.run(self.catalog)


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()

