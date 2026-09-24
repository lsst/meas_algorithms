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
import unittest.mock

import numpy as np
from astropy.table import Table

import lsst.afw.detection as afwDet
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.geom
import lsst.pex.config
import lsst.scarlet.lite as scl
import lsst.utils.tests
from lsst.meas.algorithms import MultiResolutionDetectionConfig, MultiResolutionDetectionTask
from lsst.meas.algorithms import multiResolutionDetection
from lsst.meas.algorithms.multiResolutionDetection import SIGMA_TO_FWHM
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

    def _makeExtendedExposure(self) -> afwImage.Exposure:
        """Plant one compact and one broad source, at a non-zero origin."""
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(1000, 2000), lsst.geom.Extent2I(180, 180))
        coords = [(1040, 2040, 6000, self.starSigma),
                  (1110, 2110, 60000, 9.0)]
        return plantSources(bbox, 15, self.sky, coords, addPoissonNoise=True)

    def testConfigValidation(self) -> None:
        """The config must reject the background steps this task does not
        implement, and allow the base task's background re-estimation."""
        config = MultiResolutionDetectionConfig()
        self.assertFalse(config.doTempLocalBackground)
        self.assertFalse(config.doTempWideBackground)
        config.reEstimateBackground = True
        config.validate()
        for name in ("doTempLocalBackground", "doTempWideBackground"):
            config = MultiResolutionDetectionConfig()
            setattr(config, name, True)
            with self.assertRaises(lsst.pex.config.FieldValidationError):
                config.validate()

    def testExtendedSourceFootprints(self) -> None:
        """Footprints must follow the extended flux that the starlet
        coefficients suppress, and stay in the exposure's coordinate system."""
        exposure = self._makeExtendedExposure()
        config = MultiResolutionDetectionConfig()
        config.reEstimateBackground = False
        schema = afwTable.SourceTable.makeMinimalSchema()
        task = MultiResolutionDetectionTask(schema=schema, config=config)
        results = task.run(afwTable.SourceTable.make(schema), exposure)

        self.assertEqual(results.positive.getRegion(), exposure.getBBox())
        footprints = results.positive.getFootprints()
        self.assertEqual(len(footprints), 2)
        for footprint in footprints:
            self.assertTrue(exposure.getBBox().contains(footprint.getBBox()))

        starletResult = scl.detect.detect_peaks(
            images=exposure.image.array[None],
            variance=exposure.variance.array[None],
            scales=config.starletScales,
            generation=config.generation,
            first_scale=config.firstScale,
            origin=(exposure.getY0(), exposure.getX0()),
            min_separation=0,
            min_area=config.minPixels,
            peak_thresh=config.peakThreshold,
            footprint_thresh=config.filterThreshold,
            psf_fwhm=self.starSigma*SIGMA_TO_FWHM,
            kappa=config.saddleThreshold,
        )
        # The starlet coefficient footprints stop inside the source, so the
        # broadest of them is still smaller than the smallest chi^2 footprint.
        widestStarlet = max(int(fp.data.sum()) for fp in starletResult.footprints)
        self.assertGreater(min(fp.getArea() for fp in footprints), widestStarlet)

    def testPeaksAssignedToFootprints(self) -> None:
        """Every footprint must hold at least one peak, and the peak table must
        index the footprint each peak was placed in."""
        exposure = self._makeExtendedExposure()
        config = MultiResolutionDetectionConfig()
        config.reEstimateBackground = False
        schema = afwTable.SourceTable.makeMinimalSchema()
        task = MultiResolutionDetectionTask(schema=schema, config=config)
        results = task.run(afwTable.SourceTable.make(schema), exposure)

        footprints = results.positive.getFootprints()
        for footprint in footprints:
            self.assertGreater(len(footprint.getPeaks()), 0)
            peakValues = [peak["peakValue"] for peak in footprint.getPeaks()]
            self.assertEqual(peakValues, sorted(peakValues, reverse=True))

        self.assertIn("footprint", results.peaks.colnames)
        for row in results.peaks:
            if row["footprint"] < 0:
                continue
            footprint = footprints[row["footprint"]]
            self.assertTrue(footprint.contains(lsst.geom.Point2I(row["x"], row["y"])))
        self.assertEqual(sum(len(fp.getPeaks()) for fp in footprints),
                         int((results.peaks["footprint"] >= 0).sum()))

    def testPerBandPsfWidths(self) -> None:
        """Bands with different PSF widths are each correlated with their own
        PSF, so the convolution trims them by different amounts."""
        bands = ["g", "r", "i"]
        psfSigmas = [1.5, 2.5, 4.0]
        exposures = []
        for bandSigma in psfSigmas:
            coords = [(x, y, counts, bandSigma) for x, y, counts, _ in self.coords]
            exposure = plantSources(self.bbox, self.kwid, self.sky, coords, addPoissonNoise=True)
            exposure.setPsf(afwDet.GaussianPsf(self.kwid, self.kwid, bandSigma))
            exposures.append(exposure)
        mExposure = afwImage.MultibandExposure.fromExposures(bands, exposures)

        config = MultiResolutionDetectionConfig()
        config.reEstimateBackground = False
        schema = afwTable.SourceTable.makeMinimalSchema()
        task = MultiResolutionDetectionTask(schema=schema, config=config)

        singles = list(mExposure.singles)
        psfs = task._getDetectionPsfs(singles, None)
        radii = [psf.computeShape(psf.getAveragePosition()).getDeterminantRadius() for psf in psfs]
        self.assertFloatsAlmostEqual(np.array(radii), np.array(psfSigmas), atol=1e-6)
        self.assertGreater(len(set(psf.getDimensions().getX() for psf in psfs)), 1)

        results = task.run(afwTable.SourceTable.make(schema), mExposure)
        self.assertGreater(len(results.positive.getFootprints()), 0)
        for footprint in results.positive.getFootprints():
            self.assertGreater(len(footprint.getPeaks()), 0)

    def testPeakValueFromDetectionImage(self) -> None:
        """A peak's value must be read off the image its footprint came from,
        not carried over from the starlet significance map."""
        exposure = self._makeExtendedExposure()
        config = MultiResolutionDetectionConfig()
        config.reEstimateBackground = False
        config.nSigmaToGrow = 0
        schema = afwTable.SourceTable.makeMinimalSchema()
        task = MultiResolutionDetectionTask(schema=schema, config=config)
        results = task.run(afwTable.SourceTable.make(schema), exposure)

        peakRows = {(row["x"], row["y"]): row for row in results.peaks if row["footprint"] >= 0}
        self.assertGreater(len(peakRows), 0)
        differs = False
        for footprint in results.positive.getFootprints():
            for peak in footprint.getPeaks():
                row = peakRows[(peak.getIx(), peak.getIy())]
                # Every peak is above the footprint threshold on that image.
                self.assertGreater(peak["peakValue"], config.thresholdValue)
                if abs(peak["peakValue"] - row["peak_sigma"]) > 1e-4:
                    differs = True
        self.assertTrue(differs, "peakValue should not simply mirror peak_sigma")

    def testNoRedundantThresholding(self) -> None:
        """Each polarity is thresholded once, on the coadd built for it, even
        when background re-estimation asks for both."""
        exposure = self._makeExposure()
        cases = (("positive", True, 2), ("positive", False, 1),
                 ("both", True, 2), ("both", False, 2),
                 ("negative", True, 2), ("negative", False, 1))
        for polarity, reEstimateBackground, expected in cases:
            config = MultiResolutionDetectionConfig()
            config.thresholdPolarity = polarity
            config.reEstimateBackground = reEstimateBackground
            # Growing also builds a FootprintSet, which would be counted below.
            config.nSigmaToGrow = 0
            schema = afwTable.SourceTable.makeMinimalSchema()
            task = MultiResolutionDetectionTask(schema=schema, config=config)

            original = multiResolutionDetection.afwDet.FootprintSet
            scans = []

            def counting(*args, **kwargs):
                # The thresholding constructor takes an image and a threshold;
                # FootprintSet(bbox) takes neither.
                if len(args) > 1:
                    scans.append(args[1])
                return original(*args, **kwargs)

            with unittest.mock.patch.object(multiResolutionDetection.afwDet, "FootprintSet", counting):
                task.run(afwTable.SourceTable.make(schema), exposure.clone())
            self.assertEqual(len(scans), expected,
                             f"thresholdPolarity={polarity}, "
                             f"reEstimateBackground={reEstimateBackground}")

    def testThresholdTypeValidation(self) -> None:
        """thresholdType selects the noise the bands are divided by, so only
        the two standard deviation choices mean anything here."""
        self.assertEqual(MultiResolutionDetectionConfig().thresholdType, "stdev")
        for thresholdType in ("pixel_stdev", "stdev"):
            config = MultiResolutionDetectionConfig()
            config.thresholdType = thresholdType
            config.validate()
        for thresholdType in ("value", "variance"):
            config = MultiResolutionDetectionConfig()
            config.thresholdType = thresholdType
            with self.assertRaises(lsst.pex.config.FieldValidationError):
                config.validate()

    def testBandNoise(self) -> None:
        """The per-pixel noise follows the variance plane, and the single
        value skips the pixels with no usable variance."""
        convolved = afwImage.MaskedImageF(self.bbox)
        convolved.variance.array[:] = 4.0
        convolved.variance.array[:40] = 0.0
        convolved.mask.array[:40] |= convolved.mask.getPlaneBitMask("NO_DATA")
        schema = afwTable.SourceTable.makeMinimalSchema()

        config = MultiResolutionDetectionConfig()
        config.thresholdType = "pixel_stdev"
        task = MultiResolutionDetectionTask(schema=schema, config=config)
        perPixel = task._bandNoise(convolved)
        self.assertEqual(np.shape(perPixel), convolved.image.array.shape)
        self.assertFloatsAlmostEqual(perPixel[40:], 2.0)

        config = MultiResolutionDetectionConfig()
        config.thresholdType = "stdev"
        task = MultiResolutionDetectionTask(schema=schema, config=config)
        single = task._bandNoise(convolved)
        self.assertEqual(np.shape(single), ())
        # A raw median over the whole plane would be dragged to zero.
        self.assertFloatsAlmostEqual(float(single), 2.0)

    def testThresholdTypeReachesPeakDetection(self) -> None:
        """The starlet coefficients are standardized by the same noise as the
        footprints, so the choice reaches scarlet lite too."""
        exposure = self._makeExtendedExposure()
        expected = {"pixel_stdev": "pixel", "stdev": "median"}
        for thresholdType, varianceMode in expected.items():
            config = MultiResolutionDetectionConfig()
            config.thresholdType = thresholdType
            config.reEstimateBackground = False
            schema = afwTable.SourceTable.makeMinimalSchema()
            task = MultiResolutionDetectionTask(schema=schema, config=config)
            seen = []
            original = multiResolutionDetection.scl.detect.detect_peaks

            def recording(*args, **kwargs):
                seen.append(kwargs["variance_mode"])
                return original(*args, **kwargs)

            with unittest.mock.patch.object(multiResolutionDetection.scl.detect,
                                            "detect_peaks", recording):
                results = task.run(afwTable.SourceTable.make(schema), exposure.clone())
            self.assertEqual(seen, [varianceMode])
            # The noise the coefficients were divided by follows the choice.
            self.assertEqual(results.sigma.ndim, 4)
            perPixel = varianceMode == "pixel"
            self.assertEqual(results.sigma.shape[-2:] != (1, 1), perPixel)

    def testThresholdTypeChangesFootprints(self) -> None:
        """A depth step across the image is treated differently by the two
        noise choices, on the same pixels."""
        base = self._makeExtendedExposure()
        base.variance.array[:, :90] *= 0.25
        detected = {}
        for thresholdType in ("pixel_stdev", "stdev"):
            config = MultiResolutionDetectionConfig()
            config.thresholdType = thresholdType
            config.reEstimateBackground = False
            schema = afwTable.SourceTable.makeMinimalSchema()
            task = MultiResolutionDetectionTask(schema=schema, config=config)
            exposure = base.clone()
            results = task.run(afwTable.SourceTable.make(schema), exposure)
            self.assertGreater(len(results.positive.getFootprints()), 0)
            mask = exposure.mask
            detected[thresholdType] = int(((mask.array & mask.getPlaneBitMask("DETECTED")) > 0).sum())
        self.assertNotEqual(detected["pixel_stdev"], detected["stdev"])

    def testImageUnchanged(self) -> None:
        """Removing the bad pixels must act on a copy, leaving the exposure's
        image plane alone."""
        exposure = self._makeExposure()
        exposure.mask.array[10:20, 10:20] |= exposure.mask.getPlaneBitMask("BAD")
        original = exposure.image.array.copy()
        config = MultiResolutionDetectionConfig()
        config.reEstimateBackground = False
        config.excludeMaskPlanes = ["BAD"]
        schema = afwTable.SourceTable.makeMinimalSchema()
        task = MultiResolutionDetectionTask(schema=schema, config=config)
        task.run(afwTable.SourceTable.make(schema), exposure)

        self.assertFloatsEqual(exposure.image.array, original)

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
        self.assertIsInstance(results.background, afwMath.BackgroundList)

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
        # One background per band, so that each can be added back to its own.
        self.assertEqual(len(results.background), len(bands))
        for background in results.background:
            self.assertIsInstance(background, afwMath.BackgroundList)
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
