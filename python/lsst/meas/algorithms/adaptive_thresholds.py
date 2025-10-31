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
]

import numpy as np

from lsst.pex.config import Field, Config, DictField, FieldValidationError
from lsst.pipe.base import Struct, Task

from .detection import SourceDetectionConfig, SourceDetectionTask


class AdaptiveThresholdDetectionConfig(Config):
    """Configuration for AdaptiveThresholdDetectionTask
    """
    maxAdaptiveDetIter = Field(dtype=int, default=20,
                               doc="Maximum number of adaptive threshold detection iterations.")
    maxDynDetIter = Field(dtype=int, default=20,
                          doc="DEPRECATED...just including for config reading of prev runs...")
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

    def validate(self):
        super().validate()
        if "fallback" not in self.maxNumPeakPerBand:
            msg = ("Must include a \"fallback\" key in the config.maxNumPeakPerBand config dict. "
                   f"It is currenly: {self.maxNumPeakPerBand}.")
            raise FieldValidationError(self.__class__.maxNumPeakPerBand, self, msg)
        if self.minFootprint < self.minIsolated:
            msg = (f"The config.minFootprint (= {self.minFootprint}) must be >= that of "
                   f"config.minIsolated (= {self.minIsolated}).")
            raise FieldValidationError(self.__class__.minFootprint, self, msg)
        if self.sufficientIsolated < self.minIsolated:
            msg = (f"The config.sufficientIsolated (= {self.sufficientIsolated}) must be >= that of "
                   f"config.minIsolated (= {self.minIsolated}).")
            raise FieldValidationError(self.__class__.sufficientIsolated, self, msg)


class AdaptiveThresholdDetectionTask(Task):
    """Detection of sources on an image using an adaptive scheme for
    the detection threshold.
    """
    ConfigClass = AdaptiveThresholdDetectionConfig
    _DefaultName = "adaptiveThresholdDetection"

    def __init__(self, *args, **kwargs):
        Task.__init__(self, *args, **kwargs)

    def run(self, table, exposure, initialThreshold=None, initialThresholdMultiplier=2.0,
            doReEstimateBackground=True, backgroundToPhotometricRatio=None):
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
        initialThreshold : `float`, optional
            Initial threshold for detection of PSF sources.
        initialThresholdMultiplier : `float`, optional
            Initial threshold for detection of PSF sources.
        doReEstimateBackground: `bool`, optional
            Re-estimate the background during the final detection stage?
        backgroundToPhotometricRatio : `lsst.afw.image.Image`, optional
            Image to convert photometric-flattened image to
            background-flattened image if ``reEstimateBackground`` is `True`
            and exposure has been photometric-flattened, and
            ``config.doApplyFlatBackgroundRatio`` is `True`.

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            The adaptive threshold detection results as a struct with
            attributes:

            ``detections``
                Results of the final round of detection as a struch with
                attributes:

                ``sources``
                    Detected sources on the exposure
                    (`lsst.afw.table.SourceCatalog`).
                ``positive``
                    Positive polarity footprints
                    (`lsst.afw.detection.FootprintSet` or `None`).
                ``negative``
                    Negative polarity footprints
                    (`lsst.afw.detection.FootprintSet` or `None`).
                ``numPos``
                    Number of footprints in positive or 0 if detection polarity was
                    negative (`int`).
                ``numNeg``
                    Number of footprints in negative or 0 if detection polarity was
                    positive (`int`).
                ``background``
                    Re-estimated background. `None` if
                    ``reEstimateBackground==False``
                    (`lsst.afw.math.BackgroundList`).
                ``factor``
                    Multiplication factor applied to the configured detection
                    threshold. (`float`).
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
        inAdaptiveDetection = True
        nAdaptiveDetIter = 0
        thresholdFactor = 1.0
        if nAdaptiveDetIter == 0:
            if initialThreshold is None:
                maxSn = float(np.nanmax(exposure.image.array/np.sqrt(exposure.variance.array)))
                adaptiveDetThreshold = min(maxSn, 5.0)
            else:
                adaptiveDetThreshold = initialThreshold
            adaptiveDetectionConfig = SourceDetectionConfig()
            adaptiveDetectionConfig.thresholdValue = adaptiveDetThreshold
            adaptiveDetectionConfig.includeThresholdMultiplier = initialThresholdMultiplier
            adaptiveDetectionConfig.reEstimateBackground = False
            adaptiveDetectionConfig.doTempWideBackground = True
            adaptiveDetectionConfig.tempWideBackground.binSize = 512
            adaptiveDetectionConfig.thresholdPolarity = "both"
            self.log.info("Using adaptive detection with thresholdValue = %.2f and multiplier = %.1f",
                          adaptiveDetectionConfig.thresholdValue,
                          adaptiveDetectionConfig.includeThresholdMultiplier)
            adaptiveDetectionTask = SourceDetectionTask(config=adaptiveDetectionConfig)

        maxNumNegFactor = 1.0
        while inAdaptiveDetection:
            inAdaptiveDetection = False
            nAdaptiveDetIter += 1
            detRes = adaptiveDetectionTask.run(table=table, exposure=exposure, doSmooth=True)
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
                    inAdaptiveDetection = False
                    continue

            if nFootprint == 0 or nPosPeak == 0:
                thresholdFactor = 0.25
                maxNumNegFactor *= 10
                self.log.warning("Adaptive threshold increase went too far (nFootprint = 0). "
                                 "Decreasing threshold to %.2f and rerunning.",
                                 thresholdFactor*adaptiveDetectionConfig.thresholdValue)
                adaptiveDetectionConfig.thresholdValue = (
                    thresholdFactor*adaptiveDetectionConfig.thresholdValue)
                adaptiveDetectionTask = SourceDetectionTask(config=adaptiveDetectionConfig)
                inAdaptiveDetection = False if nAdaptiveDetIter >= self.config.maxAdaptiveDetIter else True
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
                inAdaptiveDetection = False if nAdaptiveDetIter >= self.config.maxAdaptiveDetIter else True
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
                    inAdaptiveDetection = (
                        False if nAdaptiveDetIter >= self.config.maxAdaptiveDetIter else True
                    )
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
                inAdaptiveDetection = True
            else:
                inAdaptiveDetection = False
            if nAdaptiveDetIter >= self.config.maxAdaptiveDetIter:
                inAdaptiveDetection = False
        # Final round of detection with positive polarity
        adaptiveDetectionConfig.thresholdPolarity = "positive"
        if doReEstimateBackground:
            adaptiveDetectionConfig.reEstimateBackground = True
        adaptiveDetectionTask = SourceDetectionTask(config=adaptiveDetectionConfig)
        self.log.info("Perfomring final round of detection with threshold %.2f and multiplier %.1f",
                      adaptiveDetectionConfig.thresholdValue,
                      adaptiveDetectionConfig.includeThresholdMultiplier)
        detRes = adaptiveDetectionTask.run(table=table, exposure=exposure, doSmooth=True,
                                           backgroundToPhotometricRatio=backgroundToPhotometricRatio)
        return Struct(
            detections=detRes,
            thresholdValue=adaptiveDetectionConfig.thresholdValue,
            includeThresholdMultiplier=adaptiveDetectionConfig.includeThresholdMultiplier,
        )
