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

from __future__ import annotations

__all__ = ("MultiResolutionDetectionConfig", "MultiResolutionDetectionTask")

from contextlib import contextmanager
from typing import Any

import numpy as np
from astropy.table import Table
from numpy.lib import recfunctions as rfn

import lsst.afw.detection as afwDet
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.geom as geom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.scarlet.lite as scl
from lsst.utils.timer import timeMethod

from .detection import SourceDetectionConfig, SourceDetectionTask


# Conversion from Gaussian sigma to full width at half maximum.
# This is also the default minimum FWHM when the exposure has no usable PSF.
SIGMA_TO_FWHM = 2.0*np.sqrt(2.0*np.log(2.0))


def scarletFootprintToAfw(footprint: scl.detect.Footprint) -> afwDet.Footprint:
    """Convert a scarlet lite footprint into an afw footprint.

    Parameters
    ----------
    footprint:
        The scarlet footprint to convert.

    Returns
    -------
    newFootprint:
        The converted afw footprint, with the scarlet peaks copied over.

    Notes
    -----
    This duplicates ``scarletFootprintToAfw`` from meas_extensions_scarlet.
    That package depends on meas_algorithms, so importing it here would create
    a circular dependency.
    """
    xy0 = geom.Point2I(int(footprint.bbox.origin[1]), int(footprint.bbox.origin[0]))
    data = afwImage.Mask(footprint.data.astype(np.int32), xy0=xy0)
    spans = afwGeom.SpanSet.fromMask(data)
    newFootprint = afwDet.Footprint(spans)
    for peak in footprint.peaks:
        newFootprint.addPeak(peak.x, peak.y, peak.flux)
    return newFootprint


class MultiResolutionDetectionConfig(SourceDetectionConfig):
    """Configuration parameters for the MultiResolutionDetectionTask.
    """
    starletScales = pexConfig.Field(
        doc="The maximum number of starlet scales to use in the multi-resolution detection.",
        dtype=int, default=3,
    )
    firstScale = pexConfig.Field(
        doc="The first starlet scale to use in the multi-resolution detection.",
        dtype=int, default=1,
    )
    generation = pexConfig.Field(
        doc="The generation of the starlet transform to use (1 or 2).",
        dtype=int, default=2,
    )
    minPixels = pexConfig.Field(
        doc="Detected footprints with fewer than the specified number of pixels will be ignored.",
        dtype=int, default=8,
    )
    peakThreshold = pexConfig.Field(
        doc="Threshold for detecting peaks within footprints, in sigma.",
        dtype=float, default=5.0,
    )
    thresholdValue = pexConfig.Field(
        doc="Threshold for detecting footprints, in sigma.",
        dtype=float, default=2.0,
    )
    saddleThreshold = pexConfig.Field(
        doc="Minimum prominence of a peak above the saddle to a brighter peak, in sigma.",
        dtype=float, default=3.0,
    )

    def setDefaults(self) -> None:
        super().setDefaults()
        # Starlets act as a compensated filter, so none of the background
        # subtraction machinery of the base task is used or wanted.
        self.reEstimateBackground = False
        self.doTempLocalBackground = False
        self.doTempWideBackground = False

    def validate(self) -> None:
        super().validate()
        if self.reEstimateBackground:
            raise pexConfig.FieldValidationError(
                self.__class__.reEstimateBackground, self,
                "reEstimateBackground must be False for MultiResolutionDetectionTask"
            )


class MultiResolutionDetectionTask(SourceDetectionTask):
    """Detect peaks and footprints of sources in an image.

    Parameters
    ----------
    schema:
        Schema object used to create the output `lsst.afw.table.SourceCatalog`.
    **kwds:
        Keyword arguments passed to `lsst.pipe.base.Task.__init__`.

    Notes
    -----
    This task uses the starlet transform as a compensated filter to detect
    sources at multiple scales, and (optionally) multiple bands.
    Once the starlet coefficients are built, a significance map is made that
    converts images into units of sigma and creates a chi^2 image
    (if multiple bands are used). Candidate peaks are then detected in each
    significance map using a watershed algorithm. Next the peaks are merged
    across bands and scales to produce a final catalog of detected sources.
    Finally the peaks are added to footprints and are optionally grown
    according to the configuration.

    Because the starlet detection accounts for the noise in the exposure
    itself, the exposure must not have its background re-estimated or its
    variance rescaled before detection.
    """
    ConfigClass = MultiResolutionDetectionConfig
    _DefaultName = "multiResolutionDetection"

    def __init__(self, schema: afwTable.Schema | None, **kwds: Any) -> None:
        # Skip SourceDetectionTask.__init__ so that none of the background
        # subtasks are constructed; only the negative-detection flag is needed.
        pipeBase.Task.__init__(self, **kwds)
        if schema is not None and self.config.thresholdPolarity == "both":
            self.negativeFlagKey = schema.addField(
                "is_negative", type="Flag",
                doc="Set if source peak was detected as negative."
            )
        else:
            if self.config.thresholdPolarity == "both":
                self.log.warning("Detection polarity set to 'both', but no flag will be "
                                 "set to distinguish between positive and negative detections")
            self.negativeFlagKey = None

    @timeMethod
    def run(
        self,
        table: afwTable.SourceTable,
        exposure: afwImage.Exposure | afwImage.MultibandExposure,
        doSmooth: bool = True,
        sigma: float | None = None,
        clearMask: bool = True,
        expId: int | None = None,
        background: afwMath.BackgroundList | None = None,
        backgroundToPhotometricRatio: afwImage.Image | None = None,
    ) -> pipeBase.Struct:
        r"""Detect sources and return a catalog of detections.

        Parameters
        ----------
        table:
            Table object that will be used to create the SourceCatalog.
        exposure:
            Exposure to process; the DETECTED mask plane will be set in-place.
            May be a single-band ``Exposure`` or a multi-band
            ``MultibandExposure``.
        doSmooth:
            Unused; retained for compatibility with the base
            ``SourceDetectionTask``.
        sigma:
            Gaussian sigma of the PSF (pixels), used to grow detections; if
            `None` it is measured from the PSF of the ``exposure``.
        clearMask:
            Clear DETECTED{,_NEGATIVE} planes before running detection.
        expId:
            Unused; retained for compatibility with the base task.
        background:
            Unused; retained for compatibility with the base task.
        backgroundToPhotometricRatio:
            Unused; retained for compatibility with the base task.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            See `detectFootprints` for the contents, plus:

            ``sources``
                Detected sources on the exposure.
                (`lsst.afw.table.SourceCatalog`)

        Raises
        ------
        ValueError
            Raised if flags.negative is needed, but isn't in table's schema.

        Notes
        -----
        This override exists only to accept a `MultibandExposure` and to
        document the parameters that this task actually uses; the base
        implementation already delegates the detection to the overridden
        `detectFootprints`.
        """
        return super().run(
            table=table,
            exposure=exposure,
            doSmooth=doSmooth,
            sigma=sigma,
            clearMask=clearMask,
            expId=expId,
            background=background,
            backgroundToPhotometricRatio=backgroundToPhotometricRatio,
        )

    @timeMethod
    def detectFootprints(
        self,
        exposure: afwImage.Exposure | afwImage.MultibandExposure,
        doSmooth: bool = False,
        sigma: float | None = None,
        clearMask: bool = True,
        expId: int | None = None,
        background: afwMath.BackgroundList | None = None,
        backgroundToPhotometricRatio: afwImage.Image | None = None,
    ) -> pipeBase.Struct:
        """Detect footprints on an exposure.

        Parameters
        ----------
        exposure:
            Exposure to process; the DETECTED{,_NEGATIVE} mask plane will be
            set in-place. May be a single-band ``Exposure`` or a multi-band
            ``MultibandExposure``.
        doSmooth:
            Unused; retained for compatibility with the base task.
        sigma:
            Gaussian sigma of the PSF (pixels), used to grow detections; if
            `None` it is measured from the PSF of the ``exposure``.
        clearMask:
            Clear both DETECTED and DETECTED_NEGATIVE planes before running
            detection.
        expId:
            Unused; retained for compatibility with the base task.
        background:
            Unused; retained for compatibility with the base task.
        backgroundToPhotometricRatio:
            Unused; retained for compatibility with the base task.

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            See `scarletToAfwResults` for the contents. In addition
            `finalizeFootprints` adds the ``numPos``, ``numPosPeaks``,
            ``numNeg`` and ``numNegPeaks`` counts.
        """
        multiband = hasattr(exposure, "bands")
        x0, y0 = exposure.getBBox().getMin()

        if multiband:
            singles = list(exposure.singles)
            if clearMask:
                for single in singles:
                    self.clearMask(single.mask)
                    self.removeBadPixels(single.maskedImage)
            fwhm = self._getFwhm(singles, sigma)
            images = exposure.image.array
            variance = exposure.variance.array
            mask = singles[0].mask
        else:
            maskedImage = exposure.maskedImage
            if clearMask:
                self.clearMask(maskedImage.getMask())
            self.removeBadPixels(maskedImage)
            fwhm = self._getFwhm([exposure], sigma)
            images = exposure.image.array[None]
            variance = exposure.variance.array[None]
            mask = maskedImage.getMask()

        positive = None
        negative = None
        if self.config.thresholdPolarity in ("positive", "both"):
            positive = self._detectPeaks(images, variance, (y0, x0), fwhm)
        if self.config.thresholdPolarity in ("negative", "both"):
            negative = self._detectPeaks(-images, variance, (y0, x0), fwhm)

        results = self.scarletToAfwResults(exposure, positive, negative)

        # Grow the footprints and set the DETECTED mask planes. The footprints
        # live in a single coordinate system, so growing once is enough; the
        # mask plane is then set on every band.
        sigmaToGrow = fwhm/SIGMA_TO_FWHM
        self.finalizeFootprints(mask, results, sigmaToGrow, factor=1.0, factorNeg=1.0)
        if multiband:
            for single in singles[1:]:
                for polarity, maskName in (("positive", "DETECTED"), ("negative", "DETECTED_NEGATIVE")):
                    fpSet = getattr(results, polarity)
                    if fpSet is not None:
                        fpSet.setMask(single.mask, maskName)
                self.clearUnwantedResults(single.mask, results)

        self.clearUnwantedResults(mask, results)

        return results

    def _getFwhm(self, singles: list[afwImage.Exposure], sigma: float | None) -> float:
        """Measure the PSF FWHM to use for growing and linking detections.

        Parameters
        ----------
        singles:
            The single-band exposures to measure the PSF from. For a multi-band
            exposure the widest PSF is used, since it sets the most conservative
            linking radius across bands.
        sigma:
            Gaussian sigma of the PSF (pixels); if `None` it is measured from
            each exposure's PSF.

        Returns
        -------
        fwhm:
            The PSF full width at half maximum, in pixels.
        """
        if sigma is not None:
            return sigma*SIGMA_TO_FWHM
        sigmas = []
        for single in singles:
            try:
                psf = self.getPsf(single, sigma=None)
                sigmas.append(psf.computeShape(psf.getAveragePosition()).getDeterminantRadius())
            except Exception:
                continue
        if not sigmas:
            return SIGMA_TO_FWHM
        return max(max(sigmas)*SIGMA_TO_FWHM, SIGMA_TO_FWHM)

    def _detectPeaks(
        self,
        images: np.ndarray,
        variance: np.ndarray,
        origin: tuple[int, int],
        fwhm: float,
    ) -> scl.detect.PeakDetectionResult:
        """Run the scarlet lite peak detection on a set of images.

        Parameters
        ----------
        images:
            The images to detect on, with shape ``(n_bands, Ny, Nx)``.
        variance:
            The per-pixel variance, with the same shape as ``images``.
        origin:
            The ``(y, x)`` location of the lower corner of the images.
        fwhm:
            The PSF full width at half maximum, in pixels.

        Returns
        -------
        result:
            The detected peaks and the intermediate detection products.
        """
        return scl.detect.detect_peaks(
            images=images,
            variance=variance,
            scales=self.config.starletScales,
            generation=self.config.generation,
            first_scale=self.config.firstScale,
            origin=origin,
            min_separation=0,
            min_area=self.config.minPixels,
            peak_thresh=self.config.peakThreshold,
            footprint_thresh=self.config.thresholdValue,
            psf_fwhm=fwhm,
            kappa=self.config.saddleThreshold,
        )

    def scarletToAfwResults(
        self,
        exposure: afwImage.Exposure | afwImage.MultibandExposure,
        positive: scl.detect.PeakDetectionResult | None = None,
        negative: scl.detect.PeakDetectionResult | None = None,
    ) -> pipeBase.Struct:
        """Convert scarlet detection results into an afw detection Struct.

        Parameters
        ----------
        exposure:
            The exposure detection was run on, used for the footprint region.
        positive:
            The positive polarity detection result.
        negative:
            The negative polarity detection result.

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            A `~lsst.pipe.base.Struct` containing:

            ``positive``
                Positive polarity footprints.
                (`lsst.afw.detection.FootprintSet` or `None`)
            ``negative``
                Negative polarity footprints.
                (`lsst.afw.detection.FootprintSet` or `None`)
            ``factor``
                Multiplication factor applied to the positive threshold; always
                1.0. (`float`)
            ``factorNeg``
                Multiplication factor applied to the negative threshold; always
                1.0. (`float`)
            ``background``
                Always `None`; this task does not estimate a background.
            ``peaks``
                The final detections, one row per source, with a ``polarity``
                column. (`astropy.table.Table`)
            ``candidates``
                The peak candidates, before grouping, with a ``polarity``
                column. (`astropy.table.Table`)
            ``positions``
                The unique candidate positions, with a ``polarity`` column.
                (`astropy.table.Table`)
            ``significanceMap``
                The per-scale significance map of the positive detection (or the
                negative detection if only negative was run). (`numpy.ndarray`
                or `None`)
            ``starlets``
                The multiband starlet coefficients, as for ``significanceMap``.
                (`numpy.ndarray` or `None`)
            ``sigma``
                The per-scale, per-band coefficient noise, as for
                ``significanceMap``. (`numpy.ndarray` or `None`)
        """
        bbox = exposure.getBBox()
        results = pipeBase.Struct(
            positive=self._makeFootprintSet(positive, bbox),
            negative=self._makeFootprintSet(negative, bbox),
            factor=1.0,
            factorNeg=1.0,
            background=None,
            peaks=self._joinPolarities(positive, negative, "peaks"),
            candidates=self._joinPolarities(positive, negative, "candidates"),
            positions=self._joinPolarities(positive, negative, "positions"),
        )
        primary = positive if positive is not None else negative
        if primary is not None:
            results.significanceMap = primary.significance_map
            results.starlets = primary.starlets
            results.sigma = primary.sigma
        else:
            results.significanceMap = None
            results.starlets = None
            results.sigma = None
        return results

    @staticmethod
    def _makeFootprintSet(
        result: scl.detect.PeakDetectionResult | None,
        bbox: geom.Box2I,
    ) -> afwDet.FootprintSet | None:
        """Build an afw FootprintSet from a scarlet detection result.

        Parameters
        ----------
        result:
            The detection result to convert.
        bbox:
            The region the footprints were detected in.

        Returns
        -------
        footprintSet:
            The footprints, or `None` if ``result`` is `None`.
        """
        if result is None:
            return None
        footprintSet = afwDet.FootprintSet(bbox)
        footprintSet.setFootprints([scarletFootprintToAfw(fp) for fp in result.footprints])
        return footprintSet

    @staticmethod
    def _joinPolarities(
        positive: scl.detect.PeakDetectionResult | None,
        negative: scl.detect.PeakDetectionResult | None,
        attr: str,
    ) -> Table:
        """Join a structured-array product across detection polarities.

        Parameters
        ----------
        positive:
            The positive polarity detection result.
        negative:
            The negative polarity detection result.
        attr:
            The name of the structured-array attribute to join (``peaks``,
            ``candidates`` or ``positions``).

        Returns
        -------
        joined:
            The concatenated rows as an `astropy.table.Table` with an added
            ``polarity`` column (``1`` for positive, ``-1`` for negative). The
            table is empty if neither result is set.
        """
        arrays = []
        for result, polarity in ((positive, 1), (negative, -1)):
            if result is None:
                continue
            array = getattr(result, attr)
            array = rfn.append_fields(
                array, "polarity", np.full(len(array), polarity, dtype=np.int32), usemask=False
            )
            arrays.append(array)
        if not arrays:
            return Table()
        if len(arrays) == 1:
            return Table(arrays[0])
        return Table(np.concatenate(arrays))

    def applyTempLocalBackground(self, exposure, middle, results):
        raise NotImplementedError("applyTempLocalBackground is not used in MultiResolutionDetectionTask")

    def convolveImage(self, maskedImage, psf, doSmooth=True):
        raise NotImplementedError("convolveImage is not used in MultiResolutionDetectionTask")

    def applyThreshold(self, middle, bbox, factor=1.0, factorNeg=None):
        raise NotImplementedError("applyThreshold is not used in MultiResolutionDetectionTask")

    def reEstimateBackground(self, maskedImage, backgrounds, backgroundToPhotometricRatio=None):
        raise NotImplementedError("reEstimateBackground is not used in MultiResolutionDetectionTask")

    def setPeakSignificance(self, exposure, footprints, threshold, negative=False):
        raise NotImplementedError("setPeakSignificance is not used in MultiResolutionDetectionTask")

    def makeThreshold(self, image, thresholdParity, factor=1.0):
        raise NotImplementedError("makeThreshold is not used in MultiResolutionDetectionTask")

    def updatePeaks(self, fpSet, image, threshold):
        raise NotImplementedError("updatePeaks is not used in MultiResolutionDetectionTask")

    @contextmanager
    def tempWideBackgroundContext(self, exposure):
        raise NotImplementedError("tempWideBackgroundContext is not used in MultiResolutionDetectionTask")
