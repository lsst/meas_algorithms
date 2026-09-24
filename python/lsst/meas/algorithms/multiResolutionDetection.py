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
from typing import cast

import numpy as np
from astropy.table import Table
from numpy.lib import recfunctions as rfn

import lsst.afw.detection as afwDet
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

# The scarlet lite name for each threshold type, so that the peaks and the
# footprints are standardized by the same noise.
VARIANCE_MODE = {"pixel_stdev": "pixel", "stdev": "median"}


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
        doc="Threshold for detecting footprints on the chi^2 detection image, in sigma. Whether that "
            "is the per-pixel or the per-band sigma is set by thresholdType.",
        dtype=float, default=2.0,
    )
    filterThreshold = pexConfig.Field(
        doc="Threshold for the starlet coefficient footprints that peaks are detected on, in sigma.",
        dtype=float, default=2.0,
    )
    saddleThreshold = pexConfig.Field(
        doc="Minimum prominence of a peak above the saddle to a brighter peak, in sigma.",
        dtype=float, default=3.0,
    )

    def setDefaults(self) -> None:
        super().setDefaults()
        # Peaks come from the starlet coefficients, which are unaffected by a
        # low frequency background, so neither local nor wide temporary
        # background subtraction is implemented here. However, the background
        # re-estimation of the base task is used, since footprints are
        # detected on the image itself.
        self.doTempLocalBackground = False
        self.doTempWideBackground = False

        # Unlike SourceDetectionTask, we want the default thresholdType to
        # be 'stdev', since the detection gains on starlet coefficients by
        # using `pixel_stdev` is minimal at the cost of significantly more
        # memory usage.
        self.thresholdType = "stdev"

    def validate(self) -> None:
        super().validate()
        for name in ("doTempLocalBackground", "doTempWideBackground"):
            if getattr(self, name):
                raise pexConfig.FieldValidationError(
                    getattr(self.__class__, name), self,
                    f"{name} must be False for MultiResolutionDetectionTask"
                )
        if self.thresholdType not in ("pixel_stdev", "stdev"):
            raise pexConfig.FieldValidationError(
                self.__class__.thresholdType, self,
                "thresholdType must be 'pixel_stdev' or 'stdev' for MultiResolutionDetectionTask: "
                "the bands are divided by their noise before they are coadded, so the chi^2 detection "
                "image is in sigma and there is no image in counts to threshold"
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

    The footprints come from the image itself rather than from the starlet
    coefficients. Footprints built off of the compensated filter intentionally
    remove signal from lower scales, so the footprints do not capture the
    full extent of the source. Instead every band is correlated with its own
    PSF, the bands are coadded into a chi^2 detection image in units of sigma,
    and that image is thresholded as in the base task.
    Each peak is then inserted into the footprint that contains it, and
    footprints with no peak are dropped.

    References
    ----------
    .. [1] Lupton, R. H., "chi^2 coadds", 2021, section 1.2.2.

    .. [2] Szalay, A. S., Connolly, A. J., and Szokoly, G. P., "Simultaneous
    Multicolor Detection of Faint Galaxies in the Hubble Deep Field",
    The Astronomical Journal, vol. 117, no. 1, pp. 68-74, 1999.
    doi:10.1086/300689.
    """
    ConfigClass = MultiResolutionDetectionConfig
    _DefaultName = "multiResolutionDetection"

    @timeMethod
    def run(
        self,
        table: afwTable.SourceTable,
        exposure: afwImage.Exposure | afwImage.MultibandExposure,
        doSmooth: bool = True,
        sigma: float | None = None,
        clearMask: bool = True,
        expId: int | None = None,
        background: afwMath.BackgroundList | list[afwMath.BackgroundList] | None = None,
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
            If True, convolve each band with a Gaussian of width ``sigma``, or
            of the measured PSF width of ``exposure``, before the footprints
            are detected. The peaks are unaffected either way.
        sigma:
            Gaussian sigma of the PSF (pixels), used for smoothing and to grow
            footprint detections; if `None` it is measured from the PSF of the
            ``exposure``.
        clearMask:
            Clear DETECTED{,_NEGATIVE} planes before running detection.
        expId:
            Unused; retained for compatibility with the base task.
        background:
            Background that was already subtracted from the exposure; will be
            modified in-place if ``reEstimateBackground=True``. One per band
            for a ``MultibandExposure``.
        backgroundToPhotometricRatio:
            Image to convert photometric-flattened image to
            background-flattened image if ``reEstimateBackground=True`` and
            exposure has been photometric-flattened.

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
        doSmooth: bool = True,
        sigma: float | None = None,
        clearMask: bool = True,
        expId: int | None = None,
        background: afwMath.BackgroundList | list[afwMath.BackgroundList] | None = None,
        backgroundToPhotometricRatio: afwImage.Image | None = None,
        factor: float = 1.0,
        factorNeg: float | None = None,
    ) -> pipeBase.Struct:
        """Detect footprints on an exposure.

        Parameters
        ----------
        exposure:
            Exposure to process; the DETECTED{,_NEGATIVE} mask plane will be
            set in-place. May be a single-band ``Exposure`` or a multi-band
            ``MultibandExposure``.
        doSmooth:
            If True, smooth the image before detection using a Gaussian
            of width ``sigma``, or the measured PSF width of ``exposure``.
            Set to False when running on e.g. a pre-convolved image, or a mask
            plane.
        sigma:
            Gaussian sigma of the PSF (pixels), used to grow detections; if
            `None` it is measured from the PSF of the ``exposure``.
        clearMask:
            Clear both DETECTED and DETECTED_NEGATIVE planes before running
            detection.
        expId:
            Unused; retained for compatibility with the base task.
        background:
            Background that was already subtracted from the exposure; will be
            modified in-place if ``reEstimateBackground=True``.
        backgroundToPhotometricRatio:
            Image to convert photometric-flattened image to
            background-flattened image if ``reEstimateBackground=True`` and
            exposure has been photometric-flattened.
        factor:
            Multiplier for the configured threshold for positive detection
            polarity.
        factorNeg:
            Multiplier for the configured threshold for negative detection
            polarity. If `None`, will be set equal to ``factor`` (i.e. equal
            to the factor used for positive detection polarity).

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            See `_mergePeaksAndFootprints` for the contents. In addition
            `finalizeFootprints` adds the ``numPos``, ``numPosPeaks``,
            ``numNeg`` and ``numNegPeaks`` counts.
        """
        isMultiband = hasattr(exposure, "bands")
        x0, y0 = exposure.getBBox().getMin()

        # The masked images are views, so setting their mask planes sets them
        # on the caller's exposure. Everything that modifies pixels below works
        # on a copy instead.
        if isMultiband:
            singles = list(cast(afwImage.MultibandExposure, exposure).singles)
        else:
            singles = [cast(afwImage.Exposure, exposure)]
        maskedImages = [single.maskedImage for single in singles]

        if clearMask:
            for maskedImage in maskedImages:
                self.clearMask(maskedImage.getMask())

        psfs = self._getDetectionPsfs(singles, sigma)
        # The widest band sets the growth radius and the radius that links
        # peaks across scales, which is the most conservative choice.
        sigma = max(psf.computeShape(psf.getAveragePosition()).getDeterminantRadius() for psf in psfs)
        fwhm = sigma*SIGMA_TO_FWHM

        # Instead of calling ``SourceDetectionTask.removeBadPixels``
        # individually for each band, here it is applied to every band at
        # once. Only the image is copied, since it is the only plane the peak
        # detection modifies.
        if isMultiband:
            images = exposure.image.array.copy()
            variance = exposure.variance.array
            masks = exposure.mask.array
        else:
            images = exposure.image.array[None].copy()
            variance = exposure.variance.array[None]
            masks = exposure.mask.array[None]
        badPixelMask = maskedImages[0].mask.getPlaneBitMask(self.config.excludeMaskPlanes)
        if badPixelMask:
            images[(masks & badPixelMask) > 0] = 0

        # Find peaks using starlets as compensated filters
        positive = None
        negative = None
        if self.config.thresholdPolarity in ("positive", "both"):
            positive = self._detectPeaks(images, variance, (y0, x0), fwhm)
        if self.config.thresholdPolarity in ("negative", "both"):
            negative = self._detectPeaks(-images, variance, (y0, x0), fwhm)

        # Find footprints on the image itself, where extended sources keep the
        # flux that the starlet coefficients suppress.
        footprintResults = self._detectChi2Footprints(
            maskedImages, psfs, doSmooth, factor, factorNeg, background
        )

        results = self._mergePeaksAndFootprints(exposure, footprintResults, positive, negative)

        # Grow the footprints and set the DETECTED mask planes. The footprints
        # live in a single coordinate system, so growing once is enough; the
        # mask plane is then set on every band.
        self.finalizeFootprints(maskedImages[0].mask, results, sigma, factor=factor, factorNeg=factorNeg)
        for maskedImage in maskedImages[1:]:
            for polarity, maskName in (("positive", "DETECTED"), ("negative", "DETECTED_NEGATIVE")):
                fpSet = getattr(results, polarity)
                if fpSet is not None:
                    fpSet.setMask(maskedImage.mask, maskName)

        if self.config.reEstimateBackground:
            for maskedImage, backgrounds in zip(maskedImages, results.background):
                self.reEstimateBackground(
                    maskedImage,
                    backgrounds,
                    backgroundToPhotometricRatio=backgroundToPhotometricRatio,
                )

        # Drop the unwanted polarity only once every band has had its mask
        # plane set, since the first call clears that polarity from ``results``.
        for maskedImage in maskedImages:
            self.clearUnwantedResults(maskedImage.getMask(), results)

        if not isMultiband:
            results.background = results.background[0]

        return results

    def _getDetectionPsfs(
        self,
        singles: list[afwImage.Exposure],
        sigma: float | None,
    ) -> list[afwDet.GaussianPsf]:
        """Build the Gaussian PSF to convolve each band with.

        Parameters
        ----------
        singles:
            The single-band exposures to measure the PSF from.
        sigma:
            Gaussian sigma of the PSF (pixels); if `None` it is measured from
            each exposure's own PSF.

        Returns
        -------
        psfs:
            One PSF per band, in band order.

        Notes
        -----
        The chi^2 coadd weights each band by the variance of its own matched
        filter, so each band is correlated with its own PSF.
        See Lupton, chi^2 coadds, equation 1.16.
        """
        psfs = []
        for single in singles:
            bandSigma = sigma
            if bandSigma is None:
                psf = single.getPsf()
                if psf is not None:
                    bandSigma = psf.computeShape(psf.getAveragePosition()).getDeterminantRadius()
                if bandSigma is None or not np.isfinite(bandSigma):
                    bandSigma = 1.0
                # A Gaussian narrower than a pixel is not representable on the
                # pixel grid.
                bandSigma = max(bandSigma, 1.0)
            psfs.append(self.getPsf(single, sigma=bandSigma))
        return psfs

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
            footprint_thresh=self.config.filterThreshold,
            psf_fwhm=fwhm,
            kappa=self.config.saddleThreshold,
            variance_mode=VARIANCE_MODE[self.config.thresholdType],
        )

    def _bandNoise(self, convolvedImage: afwImage.MaskedImage) -> np.ndarray | float:
        """Measure the noise of one convolved band.

        Parameters
        ----------
        convolvedImage:
            The convolved masked image of a single band.

        Returns
        -------
        noise:
            The standard deviation to divide the band by, per pixel for
            ``thresholdType="pixel_stdev"`` and a single value for
            ``thresholdType="stdev"``.
        """
        if self.config.thresholdType == "pixel_stdev":
            return np.sqrt(convolvedImage.variance.array)
        # The median has to skip the pixels with no usable variance, which
        # would otherwise drag it to zero over a region with no data.
        statsControl = afwMath.StatisticsControl()
        statsControl.setAndMask(convolvedImage.mask.getPlaneBitMask(self.config.statsMask))
        statistics = afwMath.makeStatistics(convolvedImage.variance, convolvedImage.mask,
                                            afwMath.MEDIAN, statsControl)
        return np.sqrt(statistics.getValue(afwMath.MEDIAN))

    def _standardizeBands(self, convolvedImages: list[afwImage.MaskedImage]) -> np.ndarray:
        """Divide each convolved band by the noise of its own matched filter.

        Parameters
        ----------
        convolvedImages:
            The convolved masked images, one per band, all sharing a bounding
            box.

        Returns
        -------
        standardized:
            The per-band planes, with shape ``(nBands, Ny, Nx)``, each in
            units of its own noise.
        """
        standardized = np.zeros((len(convolvedImages),) + convolvedImages[0].image.array.shape,
                                dtype=np.float32)
        for index, convolvedImage in enumerate(convolvedImages):
            # A pixel with no variance carries no information, and dividing
            # by it would spread NaN across the whole coadd.
            noise = self._bandNoise(convolvedImage)
            np.divide(convolvedImage.image.array, noise, out=standardized[index], where=noise > 0)
        return standardized

    @staticmethod
    def _assembleChi2Image(
        convolvedImages: list[afwImage.MaskedImage],
        standardized: np.ndarray,
    ) -> afwImage.MaskedImage:
        """Coadd the standardized bands into a chi^2 detection image.

        Parameters
        ----------
        convolvedImages:
            The convolved masked images the standardized bands came from, which
            supply the bounding box and the mask planes.
        standardized:
            The per-band planes from `_standardizeBands`, negated to detect
            negative sources.

        Returns
        -------
        chi2Image:
            The coadd, in units of Gaussian sigma, with the union of the
            per-band masks.

        Notes
        -----
        The variance plane of the coadd is left empty. Clipping the per-band
        planes at zero means the chi^2 coadd has no meaningful per-pixel
        variance; only the upper tail is calibrated, which is what
        `build_chi2_significance` maps to sigma and is all that a positive
        threshold uses.
        """
        chi2Image = afwImage.MaskedImageF(convolvedImages[0].getBBox())
        for convolvedImage in convolvedImages:
            chi2Image.mask.array |= convolvedImage.mask.array
        chi2Image.image.array[:] = scl.detect.build_chi2_significance(standardized)
        return chi2Image

    def makeThreshold(
        self,
        image: afwImage.MaskedImage,
        thresholdParity: str,
        factor: float = 1.0,
    ) -> afwDet.Threshold:
        """Make the threshold to apply to a chi^2 detection image.

        Parameters
        ----------
        image:
            Unused; the threshold does not depend on the image.
        thresholdParity:
            One of "positive" or "negative", to set the kind of fluctuations
            the threshold will detect.
        factor:
            Factor by which to multiply the configured detection threshold.

        Returns
        -------
        threshold:
            Detection threshold.

        Notes
        -----
        The bands are divided by their noise in `_standardizeBands`, so by the
        time the chi^2 coadd exists it is already in the units that
        ``config.thresholdType`` names and the threshold is applied to it
        directly. The base task instead divides here, because it thresholds an
        image in counts.
        """
        threshold = afwDet.createThreshold(self.config.thresholdValue*factor, "value",
                                           thresholdParity != "negative")
        threshold.setIncludeMultiplier(self.config.includeThresholdMultiplier)
        self.log.debug("Detection threshold: %s", threshold)
        return threshold

    def applyThreshold(
        self,
        middle: afwImage.MaskedImage | list[afwImage.MaskedImage],
        bbox: geom.Box2I,
        factor: float = 1.0,
        factorNeg: float | None = None,
    ) -> pipeBase.Struct:
        r"""Threshold the chi^2 coadd of the convolved bands.

        Parameters
        ----------
        middle:
            The convolved masked images to coadd, one per band, all sharing a
            bounding box. A single masked image is treated as one band.
        bbox:
            Bounding box of the unconvolved image, which the returned
            `~lsst.afw.detection.FootprintSet`\ s take as their region.
        factor:
            Multiplier for the configured threshold.
        factorNeg:
            Multiplier for the configured threshold for negative detection
            polarity. If `None`, will be set equal to ``factor``.

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
                Multiplication factor applied to the configured threshold
                for positive detection polarity. (`float`)
            ``factorNeg``
                Multiplication factor applied to the configured threshold
                for negative detection polarity. (`float`)
            ``negativeThreshold``
                Threshold used for negative detection polarity. (`float`)
            ``positiveImage``
                Coadded image used for positive detection polarity.
                (`lsst.afw.image.MaskedImage` or `None`)
            ``negativeImage``
                Coadded image used for negative detection polarity.
                (`lsst.afw.image.MaskedImage` or `None`)

        Notes
        -----
        A chi^2 coadd is clipped at zero, so it has no negative side and each
        polarity needs a coadd of its own: the negative detections are the
        positive detections of the coadd built from the negated bands. Both
        are therefore thresholded on their positive side.
        """
        if factorNeg is None:
            factorNeg = factor
        self.log.info("Threshold scaling factor for positive detections is: %.3f. For negative "
                      "detections it is: %.3f", factor, factorNeg)
        convolvedImages = list(middle) if isinstance(middle, (list, tuple)) else [middle]
        results = pipeBase.Struct(
            positive=None,
            negative=None,
            factor=factor,
            factorNeg=factorNeg,
            positiveThreshold=None,
            negativeThreshold=None,
            positiveImage=None,
            negativeImage=None,
        )

        # Both polarities coadd the same standardized bands, so the negative
        # pass negates them in place rather than building them again.
        standardized = self._standardizeBands(convolvedImages)
        for polarity, maskName, polarityFactor in (("positive", "DETECTED", factor),
                                                   ("negative", "DETECTED_NEGATIVE", factorNeg)):
            if (not self.config.reEstimateBackground
                    and self.config.thresholdPolarity not in (polarity, "both")):
                continue
            if polarity == "negative":
                np.negative(standardized, out=standardized)
            image = self._assembleChi2Image(convolvedImages, standardized)
            threshold = self.makeThreshold(image, "positive", factor=polarityFactor)
            footprints = afwDet.FootprintSet(image, threshold, maskName, self.config.minPixels)
            footprints.setRegion(bbox)
            setattr(results, polarity, footprints)
            setattr(results, polarity + "Threshold", threshold)
            setattr(results, polarity + "Image", image)
        return results

    def _detectChi2Footprints(
        self,
        maskedImages: list[afwImage.MaskedImage],
        psfs: list[afwDet.GaussianPsf],
        doSmooth: bool = True,
        factor: float = 1.0,
        factorNeg: float | None = None,
        background: afwMath.BackgroundList | list[afwMath.BackgroundList] | None = None,
    ) -> pipeBase.Struct:
        r"""Convolve each band and threshold their chi^2 coadd.

        Parameters
        ----------
        maskedImages:
            The masked images to detect footprints on, one per band.
        psfs:
            The PSF to convolve each band with, in band order.
        doSmooth:
            Whether to convolve before detection.
        factor:
            The multiplication factor applied to the positive threshold.
        factorNeg:
            The multiplication factor applied to the negative threshold.
        background:
            Background that was already subtracted from the exposure; will be
            modified in-place if ``reEstimateBackground=True``. One per band,
            or a single one for a single-band exposure.

        Returns
        -------
        results: `lsst.pipe.base.Struct`
            What `applyThreshold` returns, plus:

            ``background``
                One `lsst.afw.math.BackgroundList` per band, to be filled by
                `reEstimateBackground`. (`list`)
        """
        bbox = maskedImages[0].getBBox()
        convolvedImages = []
        for maskedImage, psf in zip(maskedImages, psfs):
            middle = self.convolveImage(maskedImage, psf, doSmooth=doSmooth).middle
            self.removeBadPixels(middle)
            convolvedImages.append(middle)

        # A wider PSF loses a wider border to the convolution, so the bands can
        # come back different sizes. The coadd is defined where all of them are.
        common = geom.Box2I(convolvedImages[0].getBBox())
        for middle in convolvedImages[1:]:
            common.clip(middle.getBBox())
        if any(middle.getBBox() != common for middle in convolvedImages):
            convolvedImages = [middle.Factory(middle, common, afwImage.PARENT, False)
                               for middle in convolvedImages]

        results = self.applyThreshold(convolvedImages, bbox, factor=factor, factorNeg=factorNeg)

        if background is None:
            background = [afwMath.BackgroundList() for _ in maskedImages]
        elif isinstance(background, afwMath.BackgroundList):
            background = [background]
        results.background = background
        return results

    def _mergePeaksAndFootprints(
        self,
        exposure: afwImage.Exposure | afwImage.MultibandExposure,
        footprintResults: pipeBase.Struct,
        positive: scl.detect.PeakDetectionResult | None = None,
        negative: scl.detect.PeakDetectionResult | None = None,
    ) -> pipeBase.Struct:
        """Convert scarlet detection results into an afw detection Struct.

        Parameters
        ----------
        exposure:
            The exposure detection was run on, used for the footprint region.
        footprintResults:
            The footprints detected on the chi^2 image, from
            `_detectChi2Footprints`.
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
                Multiplication factor applied to the positive threshold.
                (`float`)
            ``factorNeg``
                Multiplication factor applied to the negative threshold.
                (`float`)
            ``background``
                One `lsst.afw.math.BackgroundList` per band, to be filled by
                `reEstimateBackground`. (`list`)
            ``peaks``
                The final detections, one row per source, with a ``polarity``
                column and a ``footprint`` column indexing the footprint the
                peak was placed in, ``-1`` for a peak that fell outside every
                footprint. (`astropy.table.Table`)
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
        # The peak tables record which footprint each peak landed in, so the
        # footprint sets have to be built before the tables are joined.
        positiveFootprints = self._insertPeaks(positive, footprintResults.positive,
                                               footprintResults.positiveImage, bbox)
        negativeFootprints = self._insertPeaks(negative, footprintResults.negative,
                                               footprintResults.negativeImage, bbox)
        results = pipeBase.Struct(
            positive=positiveFootprints,
            negative=negativeFootprints,
            factor=footprintResults.factor,
            factorNeg=footprintResults.factorNeg,
            background=footprintResults.background,
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

    def _insertPeaks(
        self,
        result: scl.detect.PeakDetectionResult | None,
        footprintSet: afwDet.FootprintSet | None,
        detectionImage: afwImage.MaskedImage | None,
        bbox: geom.Box2I,
    ) -> afwDet.FootprintSet | None:
        """Place the starlet peaks into the footprints that contain them.

        Parameters
        ----------
        result:
            The scarlet detection result holding the peaks, or `None` if this
            polarity was not detected.
        footprintSet:
            The footprints detected on the chi^2 image for this polarity. Their
            own peaks, found by `lsst.afw.detection.FootprintSet` on the way
            in, are discarded in favor of the starlet peaks.
        detectionImage:
            The chi^2 image the footprints were detected on, which supplies
            each peak's ``peakValue``.
        bbox:
            The region the footprints were detected in.

        Returns
        -------
        footprintSet:
            The footprints that hold at least one peak, brightest peak first,
            or `None` if any input is `None`.

        Notes
        -----
        A peak's ``peakValue`` is read off the detection image, as
        `lsst.afw.detection.FootprintSet` does, rather than carried over from
        the starlet significance map the peak was found in. The two are
        different images, and the peak ordering within a footprint has to
        follow the image the footprint came from. The starlet significance
        stays available as ``peak_sigma`` in the output peak table.

        ``result.peaks["footprint"]`` is overwritten with the index of the
        footprint each peak was placed in, and set to ``-1`` for a peak that
        fell outside every footprint. Those peaks reach the output table but no
        footprint, so a detection is never silently dropped.

        A footprint with no peak in it is dropped. The chi^2 threshold sits
        well below the peak threshold, so such footprints are mostly noise, and
        a parent with no peak is of no use to the deblender in any case.
        """
        if result is None or footprintSet is None or detectionImage is None:
            return None

        peaks = result.peaks
        footprints = footprintSet.getFootprints()
        merged = afwDet.FootprintSet(bbox)
        if len(peaks) == 0 or len(footprints) == 0:
            peaks["footprint"] = -1
            return merged

        # Label every footprint's pixels with its index, so that each peak can
        # be looked up in one pass rather than tested against every footprint.
        # The peaks come from the full exposure, but the detection image loses
        # a border to the convolution, so they are not all inside it.
        imageBBox = detectionImage.getBBox()
        x0, y0 = imageBBox.getMin()
        width, height = imageBBox.getDimensions()
        column = peaks["x"] - x0
        row = peaks["y"] - y0
        inside = (column >= 0) & (column < width) & (row >= 0) & (row < height)

        labels = afwImage.ImageI(imageBBox)
        for index, footprint in enumerate(footprints):
            footprint.spans.setImage(labels, index + 1, doClip=True)
        labelOf = np.full(len(peaks), -1, dtype=int)
        labelOf[inside] = labels.array[row[inside], column[inside]] - 1
        peakValue = np.zeros(len(peaks), dtype=float)
        peakValue[inside] = detectionImage.image.array[row[inside], column[inside]]

        # Brightest first within a footprint, matching the afw convention.
        order = np.lexsort((-peakValue, labelOf))
        kept: list[afwDet.Footprint] = []
        footprintOf = np.full(len(peaks), -1, dtype=peaks["footprint"].dtype)
        current = -1
        for index in order:
            label = labelOf[index]
            if label < 0:
                continue
            if label != current:
                footprints[label].getPeaks().clear()
                kept.append(footprints[label])
                current = label
            footprintOf[index] = len(kept) - 1
            kept[-1].addPeak(int(peaks["x"][index]), int(peaks["y"][index]), float(peakValue[index]))
        peaks["footprint"] = footprintOf

        merged.setFootprints(kept)
        self.log.info("Merged %d of %d peaks into %d of %d chi^2 footprints",
                      int(np.sum(footprintOf >= 0)), len(peaks), len(kept), len(footprints))
        return merged

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

    def setPeakSignificance(self, exposure, footprints, threshold, negative=False):
        raise NotImplementedError("setPeakSignificance is not used in MultiResolutionDetectionTask")

    def updatePeaks(self, fpSet, image, threshold):
        raise NotImplementedError("updatePeaks is not used in MultiResolutionDetectionTask")

    @contextmanager
    def tempWideBackgroundContext(self, exposure):
        raise NotImplementedError("tempWideBackgroundContext is not used in MultiResolutionDetectionTask")
