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

from astropy.table import Table

import lsst.afw.detection as afwDet
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.geom
import lsst.pex.config
import lsst.utils.tests
from lsst.meas.algorithms import MultiResolutionDetectionConfig, MultiResolutionDetectionTask
from lsst.meas.algorithms.testUtils import plantSources


class MultiResolutionDetectionTaskTestCase(lsst.utils.tests.TestCase):

    def setUp(self) -> None:
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(128, 128))
        self.kwid = 11
        self.sky = 100.0
        self.starSigma = 2.0
        # (x, y, counts, sigma) for each planted star.
        self.coords = [(30, 30, 5000, self.starSigma),
                       (80, 90, 8000, self.starSigma),
                       (60, 20, 3000, self.starSigma)]

    def _makeExposure(self, offset: int = 0) -> afwImage.Exposure:
        coords = [(x, y, counts + offset, sigma) for x, y, counts, sigma in self.coords]
        return plantSources(self.bbox, self.kwid, self.sky, coords, addPoissonNoise=True)

    def testConfigValidation(self) -> None:
        """The config must default to no background re-estimation and reject
        turning it back on."""
        config = MultiResolutionDetectionConfig()
        self.assertFalse(config.reEstimateBackground)
        self.assertFalse(config.doTempLocalBackground)
        self.assertFalse(config.doTempWideBackground)
        config.reEstimateBackground = True
        with self.assertRaises(lsst.pex.config.FieldValidationError):
            config.validate()

    def testSingleBand(self) -> None:
        """Detection on a single-band exposure recovers the planted sources."""
        exposure = self._makeExposure()
        schema = afwTable.SourceTable.makeMinimalSchema()
        task = MultiResolutionDetectionTask(schema=schema)
        table = afwTable.SourceTable.make(schema)
        results = task.run(table, exposure)

        self.assertGreaterEqual(len(results.sources), len(self.coords))
        self.assertEqual(results.numPos, len(results.positive.getFootprints()))
        self.assertIsInstance(results.positive, afwDet.FootprintSet)
        self.assertIsNone(results.negative)
        self.assertIsNone(results.background)

        # The scarlet intermediate products are attached to the result as
        # astropy tables.
        for attr in ("peaks", "candidates", "positions"):
            table = getattr(results, attr)
            self.assertIsInstance(table, Table)
            self.assertIn("polarity", table.colnames)
            self.assertGreater(len(table), 0)
        self.assertIsNotNone(results.significanceMap)
        self.assertIsNotNone(results.starlets)

        # The DETECTED plane is set on the exposure.
        mask = exposure.mask
        detected = (mask.array & mask.getPlaneBitMask("DETECTED")) > 0
        self.assertGreater(detected.sum(), 0)

    def testMultiBand(self) -> None:
        """Detection on a MultibandExposure produces a single merged catalog
        and sets the DETECTED plane in every band."""
        bands = ["g", "r", "i"]
        exposures = [self._makeExposure(offset=1000 * i) for i in range(len(bands))]
        mExposure = afwImage.MultibandExposure.fromExposures(bands, exposures)

        schema = afwTable.SourceTable.makeMinimalSchema()
        task = MultiResolutionDetectionTask(schema=schema)
        table = afwTable.SourceTable.make(schema)
        results = task.run(table, mExposure)

        self.assertGreaterEqual(len(results.sources), len(self.coords))
        self.assertEqual(results.significanceMap.shape[1], len(bands) + 1)
        self.assertEqual(results.starlets.shape[1], len(bands))
        for single in mExposure.singles:
            mask = single.mask
            detected = (mask.array & mask.getPlaneBitMask("DETECTED")) > 0
            self.assertGreater(detected.sum(), 0)

    def testBothPolarities(self) -> None:
        """With thresholdPolarity='both', negative footprints are returned and
        the peak table records both polarities."""
        exposure = self._makeExposure()
        config = MultiResolutionDetectionConfig()
        config.thresholdPolarity = "both"
        schema = afwTable.SourceTable.makeMinimalSchema()
        task = MultiResolutionDetectionTask(schema=schema, config=config)
        table = afwTable.SourceTable.make(schema)
        results = task.run(table, exposure)

        self.assertIsInstance(results.positive, afwDet.FootprintSet)
        self.assertIsInstance(results.negative, afwDet.FootprintSet)
        self.assertIsNotNone(task.negativeFlagKey)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
