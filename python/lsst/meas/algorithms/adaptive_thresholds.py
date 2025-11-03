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

__all__ = [
    "AdaptiveThresholdDetectionConfig",
    "AdaptiveThresholdDetectionTask",
    "AdaptiveThresholdBackgroundConfig",
    "AdaptiveThresholdBackgroundTask",
]

from contextlib import contextmanager

import numpy as np

from lsst.pex.config import Field, Config, ConfigField, DictField, FieldValidationError, ListField
from lsst.pipe.base import Task

from lsst.afw.geom import SpanSet
from lsst.afw.image import Mask
from lsst.afw.math import BackgroundList
from .detection import SourceDetectionConfig, SourceDetectionTask
from .subtractBackground import SubtractBackgroundConfig, SubtractBackgroundTask


class AdaptiveThresholdDetectionConfig(Config):
    """Configuration for AdaptiveThresholdDetectionTask
    """
    maxAdaptiveDetIter = Field(dtype=int, default=20,
                               doc="Maximum number of adaptive threshold detection iterations.")
    maxNumPeakPerBand = DictField(
        keytype=str,
        itemtype=int,
        default={
            "u": 3000,
            "g": 3000,
            "r": 3000,
            "i": 5000,
            "z": 5000,
            "y": 3000,
            "fallback": 5000,
        },
        doc="Maximum number of peaks per band. If the current band for the exposure "
        "is not included as a key in this dict, the value associated with the "
        "\"fallback\" key will be used.",
    )
    minFootprint = Field(dtype=int, default=15,
                         doc="Minimum number of footprints considered sufficient to exit the "
                         "iteration loop. This should be larger than or equal to ``minIsolated`` as "
                         "it does not take into account whether any given footprint is multi-peaked "
                         "(i.e. blended thus not providing any isolated sources for PSF modeling.)")
    maxPeakToFootRatio = Field(dtype=float, default=800.0,
                               doc="Maximum ratio of the number of peaks per footprint considered "
                               "sufficient to exit the iteration loop. This helps guard against "
                               "large contiguous footprints leaving no isolated sources for "
                               "inclusion in the PSF modeling.")
    minIsolated = Field(dtype=int, default=6,
                        doc="Minimum number of single-peaked (i.e. isolated) footprints "
                        "considered sufficient to exit the iteration loop. This should be "
                        "larger than the minimum number of sources desired for PSF modeling "
                        "since some may be rejected by the source selector.")
    sufficientIsolated = Field(dtype=int, default=100,
                               doc="Number of single-peaked (isolated) footprints considered "
                               "sufficient to exit the iteration loop. Must be larger than "
                               "``minIsolated``.")
    sufficientFractionIsolated = Field(dtype=float, default=0.45,
                                       doc="Fraction of single-peaked (isolated) footprints considered "
                                       "sufficient to exit the iteration loop.")
    baseline = ConfigField(
        "Baseline configuration for SourceDetectionTask in the absence of any iteration. "
        "All options other than thresholdPolarity, thresholdValue, and includeThresholdMultiplier "
        "are held fixed at these values.",
        SourceDetectionConfig,
    )

    def setDefaults(self):
        self.baseline.reEstimateBackground = False
        self.baseline.doTempWideBackground = True
        self.baseline.tempWideBackground.binSize = 512
        self.baseline.thresholdPolarity = "positive"  # for schema and final run.
        self.baseline.includeThresholdMultiplier = 2.0

    def validate(self):
        super().validate()
        if "fallback" not in self.maxNumPeakPerBand:
            msg = ("Must include a \"fallback\" key in the config.maxNumPeakPerBand config dict. "
                   f"It is currently: {self.maxNumPeakPerBand}.")
            raise FieldValidationError(self.__class__.maxNumPeakPerBand, self, msg)
        if self.minFootprint < self.minIsolated:
            msg = (f"The config.minFootprint (= {self.minFootprint}) must be >= that of "
                   f"config.minIsolated (= {self.minIsolated}).")
            raise FieldValidationError(self.__class__.minFootprint, self, msg)
        if self.sufficientIsolated < self.minIsolated:
            msg = (f"The config.sufficientIsolated (= {self.sufficientIsolated}) must be >= that of "
                   f"config.minIsolated (= {self.minIsolated}).")
            raise FieldValidationError(self.__class__.sufficientIsolated, self, msg)
        if self.baseline.reEstimateBackground:
            raise FieldValidationError(
                self.__class__.baseline, self,
                "Baseline detection configuration must not include background re-estimation."
            )


class AdaptiveThresholdDetectionTask(Task):
    """Detection of sources on an image using an adaptive scheme for
    the detection threshold.
    """
    ConfigClass = AdaptiveThresholdDetectionConfig
    _DefaultName = "detection"

    def __init__(self, schema=None, **kwargs):
        super().__init__(**kwargs)
        # We make a baseline SourceDetectionTask only to set up the schema.
        if schema is not None:
            SourceDetectionTask(config=self.config.baseline, schema=schema)

    def run(self, table, exposure, **kwargs):
        """Perform detection with an adaptive threshold detection scheme
        conditioned to maximize the likelihood of a successful PSF model fit
        for any given "scene".

        In particular, we'd like to be able to handle different scenes, from
        sparsely populated ones through very crowded ones, and possibly high
        fill-factor nebulosity along the way, with a single pipeline (and set
        of configs). This requires some flexibility in setting the detection
        thresholds in order to detect enough sources suitable for PSF modeling
        (e.g. crowded fields require higher thresholds to ensure the detections
        don't end up overlapping into a single or very small number of blended
        footprints.

        We first detect sources using the default threshold and multiplier.
        Then, cycling through a series of criteria based on the DETECTED mask
        planes (number of footprints, number of peaks, number of isolated
        footprints, number of peaks-per-footrpint, etc.) conditioned to identify
        a "Goldilocks Zone" where we have enough isolated peaks from which to
        measure the PSF, we iterate while adjusting the detection thresholds
        in the appropriate direction until all criteria are met (or the maximum
        number of iterations is reached).

        Parameters
        ----------
        table : `lsst.afw.table.SourceTable`
            Table object that will be used to create the SourceCatalog.
        exposure : `lsst.afw.image.Exposure`
            Exposure to process; DETECTED mask plane will be set in-place.
        **kwargs
            Forwarded to internal runs of `SourceDetectionTask`.

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            The adaptive threshold detection results.  Most fields are directly
            produced by `SourceDetectionTask.run`.  Additional fields include:

            ``thresholdValue``
                The final threshold value used to the configure the final round
                of detection (`float`).
            ``includeThresholdMultiplier``
                The final multiplication factor applied to the configured detection
                threshold. (`float`).
        """
        band = "fallback"
        if exposure.filter is not None:
            if exposure.filter.hasBandLabel():
                band = exposure.filter.bandLabel
        if band in self.config.maxNumPeakPerBand:
            maxNumPeak = self.config.maxNumPeakPerBand[band]
        else:
            maxNumPeak = self.config.maxNumPeakPerBand["fallback"]

        # Set up and configure the adaptive detection task on first iteration.
        thresholdFactor = 1.0
        adaptiveDetectionConfig = self.config.baseline.copy()
        adaptiveDetectionConfig.thresholdPolarity = "both"
        self.log.info("Using adaptive detection with thresholdValue = %.2f and multiplier = %.1f",
                      adaptiveDetectionConfig.thresholdValue,
                      adaptiveDetectionConfig.includeThresholdMultiplier)
        adaptiveDetectionTask = SourceDetectionTask(config=adaptiveDetectionConfig)

        maxNumNegFactor = 1.0
        # We use a 1-indexed iteration variable just to make logs intuitive.
        for nAdaptiveDetIter in range(1, self.config.maxAdaptiveDetIter + 1):
            detRes = adaptiveDetectionTask.run(table=table, exposure=exposure, doSmooth=True, **kwargs)
            sourceCat = detRes.sources
            nFootprint = len(sourceCat)
            nPeak = 0
            nPosPeak = detRes.numPosPeaks
            nNegPeak = detRes.numNegPeaks  # Often detected in high nebulosity scenes.
            maxNumNegPeak = max(15, int(0.025*nPosPeak))*maxNumNegFactor
            nIsolated = 0
            nPeakPerSrcMax = 0
            maxNumPeakPerSrcMax = 0.2*maxNumPeak
            for src in sourceCat:
                nPeakSrc = len(src.getFootprint().getPeaks())
                if nPeakSrc == 1:
                    nIsolated += 1
                nPeak += nPeakSrc
                nPeakPerSrcMax = max(nPeakPerSrcMax, nPeakSrc)
            if nFootprint > 0.0:
                fractionIsolated = nIsolated/nFootprint
                avgPeakPerFoot = nPeak/nFootprint
            else:
                fractionIsolated = float("nan")
                avgPeakPerFoot = float("nan")
            self.log.info("In adaptive detection iter %d: nFootprints = %d, nPosPeak = %d (max is %d), "
                          "nNegPeak = %d (max is %d), nPeak/nFoot = %.1f (max is %.1f), "
                          "nPeakPerkSrcMax = %d (max is %d), nIsolated = %d, fractionIsolated = %.2f",
                          nAdaptiveDetIter, nFootprint, nPeak, maxNumPeak, nNegPeak, maxNumNegPeak,
                          avgPeakPerFoot, self.config.maxPeakToFootRatio, nPeakPerSrcMax,
                          maxNumPeakPerSrcMax, nIsolated, fractionIsolated)
            if (nIsolated > self.config.sufficientIsolated
                    and fractionIsolated > self.config.sufficientFractionIsolated
                    and (nAdaptiveDetIter > 1 or self.config.maxAdaptiveDetIter == 1)):
                if ((nIsolated > 5.0*self.config.sufficientIsolated and nPeak < 2.5*maxNumPeak
                    and nNegPeak < 100.0*maxNumNegPeak)
                        or (nNegPeak < 5.0*maxNumNegPeak and nPeak < 1.2*maxNumPeak)):
                    self.log.info("Sufficient isolated footprints (%d > %d) and fraction of isolated "
                                  "footprints (%.2f > %.2f) for PSF modeling.  Exiting adaptive detection "
                                  "at iter: %d.",
                                  nIsolated, self.config.sufficientIsolated, fractionIsolated,
                                  self.config.sufficientFractionIsolated, nAdaptiveDetIter)
                    break

            if nFootprint == 0 or nPosPeak == 0:
                thresholdFactor = 0.25
                maxNumNegFactor *= 10
                self.log.warning("Adaptive threshold increase went too far (nFootprint = 0). "
                                 "Decreasing threshold to %.2f and rerunning.",
                                 thresholdFactor*adaptiveDetectionConfig.thresholdValue)
                adaptiveDetectionConfig.thresholdValue = (
                    thresholdFactor*adaptiveDetectionConfig.thresholdValue)
                adaptiveDetectionTask = SourceDetectionTask(config=adaptiveDetectionConfig)
                continue

            if ((nPeak/nFootprint > self.config.maxPeakToFootRatio and nIsolated < self.config.minIsolated)
                    or nNegPeak > maxNumNegPeak):
                if nNegPeak > 2*maxNumNegPeak:
                    thresholdFactor = 1.5
                else:
                    thresholdFactor = 1.25
                thresholdFactor *= adaptiveDetectionConfig.includeThresholdMultiplier
                adaptiveDetectionConfig.includeThresholdMultiplier = 1.0
                self.log.warning("Adaptive detection iter %d: catalog had nPeak/nFootprint = "
                                 "%.1f (max is %.1f) and %d negative peaks (max is %d). "
                                 "Increasing threshold to %.2f and setting multiplier "
                                 "to %.1f and rerunning.",
                                 nAdaptiveDetIter, nPeak/nFootprint, self.config.maxPeakToFootRatio,
                                 nNegPeak, maxNumNegPeak,
                                 thresholdFactor*adaptiveDetectionConfig.thresholdValue,
                                 adaptiveDetectionConfig.includeThresholdMultiplier)
                adaptiveDetectionConfig.thresholdValue = (
                    thresholdFactor*adaptiveDetectionConfig.thresholdValue)
                adaptiveDetectionTask = SourceDetectionTask(config=adaptiveDetectionConfig)
                continue

            if (nPeak > maxNumPeak or nPeakPerSrcMax > maxNumPeakPerSrcMax
                    or nFootprint <= self.config.minFootprint):
                if nPeak > maxNumPeak or nPeakPerSrcMax > 0.25*maxNumPeak:
                    if nAdaptiveDetIter < 0.5*self.config.maxAdaptiveDetIter:
                        if nPeak > 3*maxNumPeak or nPeakPerSrcMax > maxNumPeak:
                            thresholdFactor = 1.7
                        elif nPeak > 2*maxNumPeak or nPeakPerSrcMax > 0.5*maxNumPeak:
                            thresholdFactor = 1.4
                        else:
                            thresholdFactor = 1.2
                    else:
                        thresholdFactor = 1.2
                    thresholdFactor *= adaptiveDetectionConfig.includeThresholdMultiplier
                    newThresholdMultiplier = max(1.0, 0.5*adaptiveDetectionConfig.includeThresholdMultiplier)
                    adaptiveDetectionConfig.includeThresholdMultiplier = newThresholdMultiplier
                    adaptiveDetectionConfig.thresholdValue = (
                        thresholdFactor*adaptiveDetectionConfig.thresholdValue)
                    self.log.warning("Adaptive detection iter %d catalog had nPeak = %d (max = %d) "
                                     "and nPeakPerSrcMax = %d (max = %d). Increasing threshold to %.2f "
                                     "and setting multiplier to %.1f and rerunning.",
                                     nAdaptiveDetIter, nPeak, maxNumPeak, nPeakPerSrcMax, maxNumPeakPerSrcMax,
                                     adaptiveDetectionConfig.thresholdValue,
                                     adaptiveDetectionConfig.includeThresholdMultiplier)
                    adaptiveDetectionTask = SourceDetectionTask(config=adaptiveDetectionConfig)
                    continue

                if nFootprint <= self.config.minFootprint:
                    maxNumNegFactor *= 10  # Allow more -ve detections at this point.
                    thresholdFactor = min(0.85, 0.4*np.log10(10*(nFootprint + 1)))
                    self.log.warning("Adaptive detection iter %d catalog had only %d footprints. "
                                     "Lowering threshold to %.2f and increasing the allowance for "
                                     "negative detections and rerunning.",
                                     nAdaptiveDetIter, nFootprint,
                                     thresholdFactor*adaptiveDetectionConfig.thresholdValue)
                    adaptiveDetectionTask = SourceDetectionTask(config=adaptiveDetectionConfig)
                adaptiveDetectionConfig.thresholdValue = (
                    thresholdFactor*adaptiveDetectionConfig.thresholdValue
                )
                adaptiveDetectionTask = SourceDetectionTask(config=adaptiveDetectionConfig)
            else:
                break
        # Final round of detection with positive polarity
        adaptiveDetectionConfig.thresholdPolarity = "positive"
        adaptiveDetectionTask = SourceDetectionTask(config=adaptiveDetectionConfig)
        self.log.info("Perfomring final round of detection with threshold %.2f and multiplier %.1f",
                      adaptiveDetectionConfig.thresholdValue,
                      adaptiveDetectionConfig.includeThresholdMultiplier)
        detections = adaptiveDetectionTask.run(table=table, exposure=exposure, doSmooth=True, **kwargs)
        detections.thresholdValue = adaptiveDetectionConfig.thresholdValue
        detections.includeThresholdMultiplier = adaptiveDetectionConfig.includeThresholdMultiplier
        return detections



class AdaptiveThresholdBackgroundConfig(SubtractBackgroundConfig):
    detectedFractionBadMaskPlanes = ListField(
        "Mask planes to ignore when computing the detected fraction.", dtype=str,
        default=["BAD", "EDGE", "NO_DATA"]
    )
    minDetFracForFinalBg = Field(
        "Minimum detected fraction for the final background.",
        dtype=float, default=0.02
    )
    maxDetFracForFinalBg = Field(
        "Maximum detected fraction for the final background.",
        dtype=float, default=0.93
    )


class AdaptiveThresholdBackgroundTask(SubtractBackgroundTask):
    """A background subtraction task that does its own masking of detected
    sources, using an adaptive scheme that iterates until bounds on the mask
    fraction are satisfied.

    Notes
    -----
    This task is only designed for use on detector images, as it is aware of
    amplifier geometry (to deal with the fact that some amps have much higher
    noise than others, and hence very different detected-mask fractions for the
    same detection threshold.
    """

    ConfigClass = AdaptiveThresholdBackgroundConfig
    _DETECTED_MASK_PLANES = ("DETECTED", "DETECTED_NEGATIVE")

    def run(self, exposure, background=None, stats=True, statsKeys=None, backgroundToPhotometricRatio=None):
        # Restore the previously measured background and remeasure it
        # using an adaptive threshold detection iteration to ensure a
        # "Goldilocks Zone" for the fraction of detected pixels.
        if not background:
            background = BackgroundList()
            median_background = 0.0
        else:
            median_background = np.median(background.getImage().array)
        self.log.warning("Original median_background = %.2f", median_background)
        # TODO: apply backgroundToPhotometricRatio here!
        exposure.image.array += background.getImage().array

        with self._restore_mask_when_done(exposure) as original_mask:
            self._dilate_original_mask(exposure, original_mask)
            self._set_adaptive_detection_mask(exposure, median_background)
            # Do not pass the original background in, since we want to wholly
            # replace it.
            return super().run(exposure=exposure, stats=stats, statsKeys=statsKeys,
                               backgroundToPhotometricRatio=backgroundToPhotometricRatio)

    def _dilate_original_mask(self, exposure, original_mask):
        nPixToDilate = 10
        detected_fraction_orig = self._compute_mask_fraction(exposure.mask)
        # Dilate the current detected mask planes and don't clear
        # them in the detection step.
        inDilating = True
        while inDilating:
            dilatedMask = original_mask.clone()
            for maskName in self._DETECTED_MASK_PLANES:
                # Compute the grown detection mask plane using SpanSet
                detectedMaskBit = dilatedMask.getPlaneBitMask(maskName)
                detectedMaskSpanSet = SpanSet.fromMask(dilatedMask, detectedMaskBit)
                detectedMaskSpanSet = detectedMaskSpanSet.dilated(nPixToDilate)
                detectedMaskSpanSet = detectedMaskSpanSet.clippedTo(dilatedMask.getBBox())
                # Clear the detected mask plane
                detectedMask = dilatedMask.getMaskPlane(maskName)
                dilatedMask.clearMaskPlane(detectedMask)
                # Set the mask plane to the dilated one
                detectedMaskSpanSet.setMask(dilatedMask, detectedMaskBit)

            detected_fraction_dilated = self._compute_mask_fraction(dilatedMask)
            if detected_fraction_dilated < self.config.maxDetFracForFinalBg or nPixToDilate == 1:
                inDilating = False
            else:
                nPixToDilate -= 1
        exposure.mask = dilatedMask
        self.log.warning("detected_fraction_orig = %.3f  detected_fraction_dilated = %.3f",
                            detected_fraction_orig, detected_fraction_dilated)
        n_above_max_per_amp = -99
        highest_detected_fraction_per_amp = float("nan")
        doCheckPerAmpDetFraction = True
        if doCheckPerAmpDetFraction:  # detected_fraction < maxDetFracForFinalBg:
            n_above_max_per_amp, highest_detected_fraction_per_amp, no_zero_det_amps = \
                self._compute_per_amp_fraction(exposure, detected_fraction_dilated)
            self.log.warning("Dilated mask: n_above_max_per_amp = %d, "
                                "highest_detected_fraction_per_amp = %.3f",
                                n_above_max_per_amp, highest_detected_fraction_per_amp)

    def _set_adaptive_detection_mask(self, exposure, median_background):
        inBackgroundDet = True
        detected_fraction = 1.0
        maxIter = 40
        nIter = 0
        nFootprintTemp = 1e12
        starBackgroundDetectionConfig = SourceDetectionConfig()
        starBackgroundDetectionConfig.doTempLocalBackground = False
        starBackgroundDetectionConfig.nSigmaToGrow = 70.0
        starBackgroundDetectionConfig.reEstimateBackground = False
        starBackgroundDetectionConfig.includeThresholdMultiplier = 1.0
        starBackgroundDetectionConfig.thresholdValue = max(2.0, 0.2*median_background)
        starBackgroundDetectionConfig.thresholdType = "pixel_stdev"  # "stdev"

        n_above_max_per_amp = -99
        highest_detected_fraction_per_amp = float("nan")
        doCheckPerAmpDetFraction = True

        while inBackgroundDet:
            currentThresh = starBackgroundDetectionConfig.thresholdValue
            if detected_fraction > self.config.maxDetFracForFinalBg:
                starBackgroundDetectionConfig.thresholdValue = 1.07*currentThresh
                if nFootprintTemp < 3 and detected_fraction > 0.9*self.config.maxDetFracForFinalBg:
                    starBackgroundDetectionConfig.thresholdValue = 1.2*currentThresh
            if n_above_max_per_amp > 1:
                starBackgroundDetectionConfig.thresholdValue = 1.1*currentThresh
            if detected_fraction < self.config.minDetFracForFinalBg:
                starBackgroundDetectionConfig.thresholdValue = 0.8*currentThresh
            starBackgroundDetectionTask = SourceDetectionTask(
                config=starBackgroundDetectionConfig)
            tempDetections = starBackgroundDetectionTask.detectFootprints(
                exposure=exposure, clearMask=True)
            exposure.mask |= dilatedMask
            nFootprintTemp = (
                (len(tempDetections.positive.getFootprints()) if tempDetections is not None else 0)
                + (len(tempDetections.negative.getFootprints()) if tempDetections.negative is not None else 0)
            )
            detected_fraction = self._compute_mask_fraction(exposure.mask)
            self.log.info("nIter = %d, thresh = %.2f: Fraction of pixels marked as DETECTED or "
                          "DETECTED_NEGATIVE in star_background_detection = %.3f "
                          "(max is %.3f; min is %.3f)",
                          nIter, starBackgroundDetectionConfig.thresholdValue,
                          detected_fraction, self.config.maxDetFracForFinalBg, self.config.minDetFracForFinalBg)

            n_amp = len(exposure.detector.getAmplifiers())
            if doCheckPerAmpDetFraction:  # detected_fraction < maxDetFracForFinalBg:
                n_above_max_per_amp, highest_detected_fraction_per_amp, no_zero_det_amps = \
                    self._compute_per_amp_fraction(exposure, detected_fraction)

            if not no_zero_det_amps:
                starBackgroundDetectionConfig.thresholdValue = 0.95*currentThresh
            nIter += 1
            if nIter > maxIter:
                inBackgroundDet = False

            if (detected_fraction < self.config.maxDetFracForFinalBg and detected_fraction > self.config.minDetFracForFinalBg
                    and n_above_max_per_amp < int(0.75*n_amp)
                    and no_zero_det_amps):
                if (n_above_max_per_amp < max(1, int(0.15*n_amp))
                        or detected_fraction < 0.85*self.config.maxDetFracForFinalBg):
                    inBackgroundDet = False
                else:
                    self.log.warning("Making small tweak....")
                    starBackgroundDetectionConfig.thresholdValue = 1.05*currentThresh
            self.log.warning("n_above_max_per_amp = %d (abs max is %d)", n_above_max_per_amp, int(0.75*n_amp))

        self.log.info("Fraction of pixels marked as DETECTED or DETECTED_NEGATIVE is now %.5f "
                      "(highest per amp section = %.5f)",
                      detected_fraction, highest_detected_fraction_per_amp)

        if detected_fraction > self.config.maxDetFracForFinalBg:
            exposure.mask = dilatedMask
            self.log.warning("Final fraction of pixels marked as DETECTED or DETECTED_NEGATIVE "
                             "was too large in star_background_detection = %.3f (max = %.3f). "
                             "Reverting to dilated mask from PSF detection...",
                             detected_fraction, self.config.maxDetFracForFinalBg)

    def _compute_mask_fraction(self, mask):
        """Evaluate the fraction of masked pixels in a (set of) mask plane(s).

        Parameters
        ----------
        mask : `lsst.afw.image.Mask`
            The mask on which to evaluate the fraction.

        Returns
        -------
        detected_fraction : `float`
            The calculated fraction of masked pixels
        """
        bad_pixel_mask = Mask.getPlaneBitMask(self.config.detectedFractionBadMaskPlanes)
        n_good_pix = np.sum(mask.array & bad_pixel_mask == 0)
        if n_good_pix == 0:
            detected_fraction = float("nan")
            return detected_fraction
        detected_pixel_mask = Mask.getPlaneBitMask(self._DETECTED_MASK_PLANES)
        n_detected_pix = np.sum((mask.array & detected_pixel_mask != 0)
                                & (mask.array & bad_pixel_mask == 0))
        detected_fraction = n_detected_pix/n_good_pix
        return detected_fraction

    def _compute_per_amp_fraction(self, exposure, detected_fraction):
        """Evaluate the maximum per-amplifier fraction of masked pixels.

        Parameters
        ----------
        exposure : `lsst.afw.image.ExposureF`
            The exposure on which to compute the per-amp masked fraction.
        detected_fraction : `float`
            The current detected_fraction of the detected mask planes for the
            full detector.

        Returns
        -------
        n_above_max_per_amp : `int`
            The number of amplifiers with masked fractions above a maximum
            value (set by the current full-detector ``detected_fraction``).
        highest_detected_fraction_per_amp : `float`
            The highest value of the per-amplifier fraction of masked pixels.
        no_zero_det_amps : `bool`
            A boolean representing whether any of the amplifiers has zero
            masked pixels.
        """
        highest_detected_fraction_per_amp = -9.99
        n_above_max_per_amp = 0
        n_no_zero_det_amps = 0
        no_zero_det_amps = True
        amps = exposure.detector.getAmplifiers()
        if amps is not None:
            for ia, amp in enumerate(amps):
                amp_bbox = amp.getBBox()
                exp_bbox = exposure.getBBox()
                if not exp_bbox.contains(amp_bbox):
                    self.log.info("Bounding box of amplifier (%s) does not fit in exposure's "
                                  "bounding box (%s).  Skipping...", amp_bbox, exp_bbox)
                    continue
                sub_image = exposure.subset(amp.getBBox())
                detected_fraction_amp = self._compute_mask_fraction(sub_image.mask)
                self.log.debug("Current detected fraction for amplifier %s = %.3f",
                               amp.getName(), detected_fraction_amp)
                if detected_fraction_amp < 0.002:
                    n_no_zero_det_amps += 1
                    if n_no_zero_det_amps > 2:
                        no_zero_det_amps = False
                        break
                highest_detected_fraction_per_amp = max(detected_fraction_amp,
                                                        highest_detected_fraction_per_amp)
                if highest_detected_fraction_per_amp > min(0.998, max(0.8, 3.0*detected_fraction)):
                    n_above_max_per_amp += 1
                    if n_above_max_per_amp > 2:
                        break
        else:
            self.log.info("No amplifier object for detector %d, so skipping per-amp "
                          "detection fraction checks.", exposure.detector.getId())
        return n_above_max_per_amp, highest_detected_fraction_per_amp, no_zero_det_amps

    @contextmanager
    def _restore_mask_when_done(self, exposure):
        original_mask = exposure.mask.clone()
        try:
            yield original_mask
        finally:
            exposure.mask = original_mask
