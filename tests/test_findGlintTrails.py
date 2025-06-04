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

from lsst.meas.algorithms import findGlintTrails
import lsst.geom
import lsst.meas.base.tests
import lsst.utils.tests

# Set this and have display_ds9 setup to see the fit on the image.
display = False
if display:
    import lsst.afw.display
    display = lsst.afw.display.Display()


class TestFindGlintTrails(lsst.utils.tests.TestCase):
    """Generate two glint trails and test that they are both found, while
    outlier points are not included.
    """

    def _make_image_3_trails(self):
        """Make an image with 3 trails on it, one of them negative.
        """
        rng = np.random.default_rng(10)
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(-5, -4), lsst.geom.Point2I(1005, 1084))
        dataset = lsst.meas.base.tests.TestDataset(bbox)
        flux = 10000

        x0 = 200
        y0 = 300
        scale = 50
        # one outlier source near the first glint trail.
        dataset.addSource(instFlux=flux, centroid=lsst.geom.Point2D(x0 - 31, y0 - 27))
        # two glint trails
        for i in range(8):
            dataset.addSource(instFlux=flux,
                              centroid=lsst.geom.Point2D(i*scale + x0 + rng.random(),
                                                         i*scale + y0 + rng.random()))

        # A negative trail (e.g. came from the template) that shouldn't be
        # included in any of the fits.
        x0 = 170
        y0 = 690
        step = 70
        for i in range(6):
            dataset.addSource(instFlux=-flux,
                              centroid=lsst.geom.Point2D(i*step + x0 + rng.random(),
                                                         -i*step + y0 + rng.random()),
                              negative=True)

        # a source that shouldn't appear in the found trails
        dataset.addSource(instFlux=flux, centroid=lsst.geom.Point2D(500, 200))

        # Perpendicular trail that has a point lying on the same line as the
        # first trail.
        x0 = 380
        y0 = 760
        step = 70
        for i in range(5):
            dataset.addSource(instFlux=flux,
                              centroid=lsst.geom.Point2D(i*step + x0 + rng.random(),
                                                         -i*step + y0 + rng.random()))

        schema = dataset.makeMinimalSchema()
        schema.addField("ip_diffim_DipoleFit_classification", type="Flag")
        exposure, catalog = dataset.realize(10.0, schema=schema)
        # So that the catalog ids don't look like simple indices when debugging.
        catalog["id"] += 100
        return exposure, catalog

    def test_simple(self):
        """Test the basic operation of the glint finder on an image with two
        glint trails.
        """
        exposure, catalog = self._make_image_3_trails()

        config = findGlintTrails.FindGlintTrailsTask.ConfigClass()
        # Limit the search radius so that we can test the line-extension
        # part of the fitter.
        config.radius = 300
        # Use a tighter threshold for this simulated data, which has better
        # controlled centroids than real data.
        config.threshold = 1.5
        task = findGlintTrails.FindGlintTrailsTask(config=config)
        result = task.run(catalog)

        if display:
            display.frame = 1
            display.image(exposure, title="something")
            for i, trail in enumerate(result.trails):
                display.centroids(trail, size=5, ctype="cyan", symbol=i)

        self.assertEqual(len(result.trails), 2)
        # Note that if you add sources to the simulated catalog, the ids in
        # the fitted trails may change.
        self.assertSetEqual(set(result.trails[0]["id"]),
                            {102, 103, 104, 105, 106, 107, 108, 109, 119})
        self.assertSetEqual(set(result.trails[1]["id"]),
                            {117, 118, 119, 120, 121})
        self.assertSetEqual(result.trailed_ids, {102, 103, 104, 105, 106, 107, 108,
                            109, 117, 118, 119, 120, 121})
        # Expected lengths from the step size of the simulated trails.
        self.assertFloatsAlmostEqual(result.parameters[0].length,
                                     np.sqrt((7*50)**2 + (7*50)**2),
                                     rtol=1e-3)
        self.assertFloatsAlmostEqual(result.parameters[1].length,
                                     np.sqrt((4*70)**2 + (4*70)**2),
                                     rtol=1e-3)
        self.assertFloatsAlmostEqual(np.degrees(result.parameters[0].angle), 45.0, rtol=2e-3)
        self.assertFloatsAlmostEqual(np.degrees(result.parameters[1].angle), -45.0, rtol=2e-3)

    def test_empty_image(self):
        """Test that no trails are found on an empty image."""
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(-5, -4), lsst.geom.Point2I(1005, 1084))
        dataset = lsst.meas.base.tests.TestDataset(bbox)
        schema = dataset.makeMinimalSchema()
        schema.addField("ip_diffim_DipoleFit_classification", type="Flag")
        exposure, catalog = dataset.realize(10.0, schema=schema)

        task = findGlintTrails.FindGlintTrailsTask()
        result = task.run(catalog)

        self.assertEqual(len(result.trails), 0)
        self.assertEqual(result.trailed_ids, set())

    def test_many_points_no_trail(self):
        """Test that no trails are found on an image with many points, but no
        sets of three that could lie on the same line.
        """
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(-5, -4), lsst.geom.Point2I(1005, 1084))
        dataset = lsst.meas.base.tests.TestDataset(bbox)
        flux = 10000
        dataset.addSource(instFlux=flux, centroid=lsst.geom.Point2D(500, 500))
        dataset.addSource(instFlux=flux, centroid=lsst.geom.Point2D(480, 610))
        dataset.addSource(instFlux=flux, centroid=lsst.geom.Point2D(520, 680))
        dataset.addSource(instFlux=flux, centroid=lsst.geom.Point2D(570, 500))
        dataset.addSource(instFlux=flux, centroid=lsst.geom.Point2D(640, 630))
        dataset.addSource(instFlux=flux, centroid=lsst.geom.Point2D(600, 700))
        schema = dataset.makeMinimalSchema()
        schema.addField("ip_diffim_DipoleFit_classification", type="Flag")
        exposure, catalog = dataset.realize(10.0, schema=schema)

        task = findGlintTrails.FindGlintTrailsTask()
        result = task.run(catalog)

        self.assertEqual(len(result.trails), 0)
        self.assertEqual(result.trailed_ids, set())

    def test_rejected_trail_min_points(self):
        """Test that no trails are found on an image with an initially
        reasonable trail that is rejected during fitting for having too many
        outliers.
        """
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(-5, -4), lsst.geom.Point2I(1005, 1084))
        dataset = lsst.meas.base.tests.TestDataset(bbox)
        flux = 10000
        dataset.addSource(instFlux=flux, centroid=lsst.geom.Point2D(100, 100))
        dataset.addSource(instFlux=flux, centroid=lsst.geom.Point2D(210, 200))
        dataset.addSource(instFlux=flux, centroid=lsst.geom.Point2D(300, 310))
        dataset.addSource(instFlux=flux, centroid=lsst.geom.Point2D(400, 400))
        dataset.addSource(instFlux=flux, centroid=lsst.geom.Point2D(500, 500))
        dataset.addSource(instFlux=flux, centroid=lsst.geom.Point2D(600, 600))
        schema = dataset.makeMinimalSchema()
        schema.addField("ip_diffim_DipoleFit_classification", type="Flag")
        exposure, catalog = dataset.realize(10.0, schema=schema)

        task = findGlintTrails.FindGlintTrailsTask()
        result = task.run(catalog)

        self.assertEqual(len(result.trails), 0)
        self.assertEqual(result.trailed_ids, set())


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
