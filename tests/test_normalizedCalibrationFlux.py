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
import logging

import lsst.afw.image
import lsst.afw.table
import lsst.utils.tests
from lsst.meas.algorithms import NormalizedCalibrationFluxTask


class NormalizedCalibrationFluxTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        self.ap_name = "base_CircularApertureFlux_12_0"
        self.cg_name = "base_CompensatedTophatFlux_12"
        self.exposure = lsst.afw.image.ExposureF(1000, 1000)

    def _make_task(self, apply_only=False):
        """
        Make a normalization task for testing.

        Parameters
        ----------
        apply_only : `bool`, optional
            Configure task in apply_only mode?

        Returns
        -------
        norm_task : `lsst.meas.algorithms.NormalizedCalibrationFluxTask`
        """
        schema = lsst.afw.table.SourceTable.makeMinimalSchema()

        config = NormalizedCalibrationFluxTask.ConfigClass()
        if apply_only:
            config.do_measure_ap_corr = False

        for name in [self.ap_name, self.cg_name]:
            schema.addField(name + "_instFlux", type=float)
            schema.addField(name + "_instFluxErr", type=float)
            schema.addField(name + "_flag", type="Flag")
        schema.addField("base_Centroid_x", type=float)
        schema.addField("base_Centroid_y", type=float)
        schema.getAliasMap().set("slot_Centroid", "base_Centroid")
        schema.addField("calib_psf_used", type="Flag")
        config.measure_ap_corr.refFluxName = self.ap_name
        flux_field = self.cg_name + "_instFlux"
        err_field = flux_field + "Err"
        config.measure_ap_corr.sourceSelector["science"].signalToNoise.fluxField = flux_field
        config.measure_ap_corr.sourceSelector["science"].signalToNoise.errField = err_field

        norm_task = NormalizedCalibrationFluxTask(schema=schema, config=config)

        return norm_task

    def _make_catalog(self, schema, num_sources=100, ap_flux_offset=0.0, cg_scale=0.5):
        """Make a source catalog with optional flux offset and scaling.

        Parameters
        ----------
        schema : `lsst.afw.table.Schema`
            The table schema.
        num_sources : `int`
            Number of sources to put into the catalog.
        ap_flux_offset : `float`
            Constant offset in aperture fluxes due to "background" issues.
        cg_scale : `float`
            Flux ratio between unnormalized flux and aperture (reference) flux.
        """
        source_cat = lsst.afw.table.SourceCatalog(schema)

        x = np.random.rand(num_sources)*self.exposure.getWidth() + self.exposure.getX0()
        y = np.random.rand(num_sources)*self.exposure.getHeight() + self.exposure.getY0()
        flux = np.random.uniform(low=10000.0, high=100000.0, size=num_sources)

        source_cat.resize(num_sources)
        source_cat["slot_Centroid_x"] = x
        source_cat["slot_Centroid_y"] = y
        source_cat["calib_psf_used"] = True

        # Make a very simple error model.
        noise_per_pix = 2.0
        ap_flux_err = np.sqrt((np.pi*12**2.*noise_per_pix)**2. + (flux + ap_flux_offset))
        cg_flux_err = np.sqrt((np.pi*4**2.*noise_per_pix)**2. + flux)

        source_cat[self.cg_name + "_instFlux"] = cg_scale*flux
        source_cat[self.cg_name + "_instFluxErr"] = cg_scale*cg_flux_err
        source_cat[self.cg_name + "_flag"] = np.zeros(num_sources, dtype=bool)
        source_cat[self.ap_name + "_instFlux"] = flux + ap_flux_offset
        source_cat[self.ap_name + "_instFluxErr"] = ap_flux_err
        source_cat[self.ap_name + "_flag"] = np.zeros(num_sources, dtype=bool)

        return source_cat

    def tearDown(self):
        del self.exposure

    def testNormalizedCalibrationFlux(self):
        np.random.seed(12345)
        norm_task = self._make_task()
        catalog = self._make_catalog(norm_task.schema)

        ap_corr_map = norm_task.run(catalog=catalog, exposure=self.exposure).ap_corr_map

        self.assertEqual(
            catalog.schema.getAliasMap().get("slot_CalibFlux"),
            norm_task.config.normalized_calibflux_name,
        )

        self.assertIn("base_CompensatedTophatFlux_12_instFlux", ap_corr_map)
        self.assertIn("base_CompensatedTophatFlux_12_instFluxErr", ap_corr_map)

        # The full set should have a 1.0 ratio when the
        # aperture flux offset is 0.0
        ratio = np.mean(catalog["slot_CalibFlux_instFlux"]/catalog[self.ap_name + "_instFlux"])
        self.assertFloatsAlmostEqual(ratio, 1.0, rtol=1e-10)

        # The subset that was used should always have a 1.0 ratio.
        used = catalog["apcorr_base_CompensatedTophatFlux_12_used"]
        ratio_used = np.mean(
            catalog["slot_CalibFlux_instFlux"][used]/catalog[self.ap_name + "_instFlux"][used]
        )
        self.assertFloatsAlmostEqual(ratio_used, 1.0, rtol=1e-10)

        # The error ratios should match the input and output.
        self.assertFloatsAlmostEqual(
            catalog["slot_CalibFlux_instFluxErr"]/catalog["slot_CalibFlux_instFlux"],
            catalog[self.cg_name + "_instFluxErr"]/catalog[self.cg_name + "_instFlux"],
        )

    def testNormalizedCalibrationFluxOffset(self):
        np.random.seed(12345)
        norm_task = self._make_task()

        for offset in [-10.0, 10.0]:
            catalog = self._make_catalog(norm_task.schema, ap_flux_offset=offset)

            norm_task.run(catalog=catalog, exposure=self.exposure)

            # The full set should not have a 1.0 ratio when the
            # aperture flux offset is not 0.0
            ratio = np.mean(catalog["slot_CalibFlux_instFlux"]/catalog[self.ap_name + "_instFlux"])
            self.assertFloatsNotEqual(ratio, 1.0)
            # Whether the full set is less than or greater than 1.0 depends on
            # the sign of the background offset.
            if offset < 0.0:
                self.assertGreater(ratio, 1.0)
            else:
                self.assertLess(ratio, 1.0)

            # The subset that was used should always have a 1.0 ratio, though
            # this may be not quite zero because of the trend in the ratio
            # vs flux even at the bright end.
            used = catalog["apcorr_base_CompensatedTophatFlux_12_used"]
            ratio_used = np.median(
                catalog["slot_CalibFlux_instFlux"][used]/catalog[self.ap_name + "_instFlux"][used]
            )
            self.assertFloatsAlmostEqual(ratio_used, 1.0, rtol=1e-10)

    def testNormalizedCalibrationFluxTooFew(self):
        np.random.seed(12345)
        norm_task = self._make_task()
        catalog = self._make_catalog(norm_task.schema)

        flags = np.ones(len(catalog), dtype=bool)
        flags[0] = False
        catalog[self.cg_name + "_flag"] = flags

        ap_corr_map = norm_task.run(catalog=catalog, exposure=self.exposure).ap_corr_map

        self.assertIn("base_CompensatedTophatFlux_12_instFlux", ap_corr_map)
        self.assertIn("base_CompensatedTophatFlux_12_instFluxErr", ap_corr_map)

        # The full set should have a 1.0 ratio when the
        # aperture flux offset is 0.0
        ratio = np.mean(catalog["slot_CalibFlux_instFlux"]/catalog[self.ap_name + "_instFlux"])
        self.assertFloatsAlmostEqual(ratio, 1.0, rtol=1e-10)

        self.assertTrue(np.all(~flags == catalog["apcorr_base_CompensatedTophatFlux_12_used"]))

    def testNormalizedCalibrationFluxApplyOnly(self):
        # Run the regular task in default mode first.
        np.random.seed(12345)
        norm_task = self._make_task()
        catalog_run1 = self._make_catalog(norm_task.schema)
        exposure_run1 = self.exposure.clone()

        ap_corr_map = norm_task.run(catalog=catalog_run1, exposure=exposure_run1).ap_corr_map

        exposure_run1.info.setApCorrMap(ap_corr_map)

        # Rerun the task; we need to make sure we have the same input so re-seed.
        np.random.seed(12345)
        norm_task2 = self._make_task(apply_only=True)
        catalog_run2 = self._make_catalog(norm_task.schema)

        ap_corr_map2 = norm_task2.run(catalog=catalog_run2, exposure=exposure_run1).ap_corr_map

        # Check that the ap_corr_map and ap_corr_map2 are the same.
        self.assertEqual(set(ap_corr_map2.keys()), set(ap_corr_map.keys()))
        for key in ap_corr_map.keys():
            self.assertEqual(ap_corr_map2[key], ap_corr_map[key])

        # Check that the slot is set correctly
        self.assertEqual(
            catalog_run2.schema.getAliasMap().get("slot_CalibFlux"),
            norm_task2.config.normalized_calibflux_name,
        )

        # Check that the final normalized catalog values are the same.
        self.assertFloatsAlmostEqual(
            catalog_run2["slot_CalibFlux_instFlux"],
            catalog_run1["slot_CalibFlux_instFlux"],
        )

    def testNormalizedCalibrationFluxApplyOnlyFail(self):
        np.random.seed(12345)
        norm_task = self._make_task()
        catalog_run1 = self._make_catalog(norm_task.schema)
        exposure_run1 = self.exposure.clone()

        norm_task.run(catalog=catalog_run1, exposure=exposure_run1).ap_corr_map

        np.random.seed(12345)
        norm_task2 = self._make_task(apply_only=True)
        catalog_run2 = self._make_catalog(norm_task.schema)

        # Try without setting an aperture correction map at all.
        with self.assertLogs(level=logging.WARNING) as cm:
            _ = norm_task2.run(catalog=catalog_run2, exposure=exposure_run1)
        warnings = '\n'.join(cm.output)
        self.assertIn("does not have a valid normalization", warnings)

        # Try again after setting an incomplete aperture correction map.
        ap_corr_map_blank = lsst.afw.image.ApCorrMap()
        exposure_run1.info.setApCorrMap(ap_corr_map_blank)

        with self.assertLogs(level=logging.WARNING) as cm:
            _ = norm_task2.run(catalog=catalog_run2, exposure=exposure_run1)
        warnings = '\n'.join(cm.output)
        self.assertIn("aperture correction map is missing base_CompensatedTophatFlux_12_instFlux", warnings)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
