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

import lsst.geom
import lsst.afw.geom
import lsst.afw.detection
from lsst.meas.algorithms import MaskBrightStarsTask
import lsst.meas.base
import lsst.meas.base.tests
import lsst.utils.tests


import lsst.afw.display
display = lsst.afw.display.Display()
display.frame = 1


class MaskBrightStarsTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(-50, -50), lsst.geom.Point2I(200, 300))
        self.small_bbox = lsst.geom.Box2I(lsst.geom.Point2I(10, 20), lsst.geom.Point2I(150, 200))
        self.sky_center = lsst.geom.SpherePoint(245.0, -45.0, lsst.geom.degrees)
        self.photo_calib = 12.3
        dataset = lsst.meas.base.tests.TestDataset(bbox, crval=self.sky_center,
                                                   calibration=self.photo_calib, psfSigma=4.0)
        dataset.addSource(1e9, lsst.geom.Point2D(100, 150))
        dataset.addSource(1e9, lsst.geom.Point2D(-5.6, 10.2))
        noise = 2.0  # stddev of noise per pixel
        schema = dataset.makeMinimalSchema()
        schema.addField("truth_flux", float, "True source flux.", "nJy")
        schema.addField("centroid_x", float, "Centroid position, as in a loaded refcat")
        schema.addField("centroid_y", float, "Centroid position, as in a loaded refcat")
        self.exposure, self.truth = dataset.realize(noise=noise, schema=schema)
        for source in self.truth:
            source["truth_flux"] = self.exposure.photoCalib.instFluxToNanojansky(source["truth_instFlux"],
                                                                                 source.getCentroid())
            source["centroid_x"] = source["truth_x"]
            source["centroid_y"] = source["truth_y"]
        # Loaded refcats are SimpleCatalogs, and don't have e.g. centroid slots.
        self.catalog = lsst.afw.table.SimpleCatalog(self.truth.schema)
        self.catalog.extend(self.truth)

    def testRun(self):
        display.image(self.exposure, "original")
        # display.centroids(catalog, size=10, ctype="red", symbol="x")
        # for x in catalog:
        #     display.dot(x['id'], x.getX(), x.getY(), size=10, ctype="cyan")

        config = MaskBrightStarsTask.ConfigClass()
        config.scale = 5
        task = MaskBrightStarsTask(config=config)
        result = task.run(self.exposure, self.catalog, "truth_flux")

        display.frame += 1
        display.image(self.exposure.mask, "mask")
        import os; print(os.getpid()); import ipdb; ipdb.set_trace();

        self.assertEqual(len(result.masked_footprints), 2)
        bright = self.exposure.mask.getPlaneBitMask("BRIGHT")
        detected = self.exposure.mask.getPlaneBitMask("DETECTED")
        self.assertEqual(self.exposure.mask.array & bright, self.exposure.mask.array & detected)


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
