
__all__ = ["DynamicDetectionConfig", "DynamicDetectionTask"]

import numpy as np

from lsst.pex.config import Field, ConfigurableField
from lsst.pipe.base import Struct, NoWorkFound

from .detection import SourceDetectionConfig, SourceDetectionTask
from .skyObjects import SkyObjectsTask

from lsst.afw.detection import FootprintSet
from lsst.afw.geom import makeCdMatrix, makeSkyWcs
from lsst.afw.table import SourceCatalog, SourceTable
from lsst.meas.base import ForcedMeasurementTask

import lsst.afw.image
import lsst.afw.math
import lsst.geom as geom


class DynamicDetectionConfig(SourceDetectionConfig):
    """Configuration for DynamicDetectionTask
    """
    prelimThresholdFactor = Field(dtype=float, default=0.5,
                                  doc="Factor by which to multiply the main detection threshold "
                                  "(thresholdValue) to use for first pass (to find sky objects).")
    prelimNegMultiplier = Field(dtype=float, default=2.5,
                                doc="Multiplier for the negative (relative to positive) polarity "
                                "detections threshold to use for first pass (to find sky objects).")
    skyObjects = ConfigurableField(target=SkyObjectsTask, doc="Generate sky objects.")
    doBackgroundTweak = Field(dtype=bool, default=True,
                              doc="Tweak background level so median PSF flux of sky objects is zero?")
    minFractionSources = Field(dtype=float, default=0.02,
                               doc="Minimum fraction of the requested number of sky sources for dynamic "
                               "detection to be considered a success. If the number of good sky sources "
                               "identified falls below this threshold, a NoWorkFound error is raised so "
                               "that this dataId is no longer considered in downstream processing.")
    doBrightPrelimDetection = Field(dtype=bool, default=True,
                                    doc="Do initial bright detection pass where footprints are grown "
                                    "by brightGrowFactor?")
    brightMultiplier = Field(dtype=float, default=2000.0,
                             doc="Multiplier to apply to the prelimThresholdFactor for the "
                             "\"bright\" detections stage (want this to be large to only "
                             "detect the brightest sources).")
    brightNegFactor = Field(dtype=float, default=2.2,
                            doc="Factor by which to multiply the threshold for the negative polatiry "
                            "detections for the \"bright\" detections stage (this needs to be fairly "
                            "low given the nature of the negative polarity detections in the very "
                            "large positive polarity threshold).")
    brightGrowFactor = Field(dtype=int, default=40,
                             doc="Factor by which to grow the footprints of sources detected in the "
                             "\"bright\" detections stage (want this to be large to mask wings of "
                             "bright sources).")
    brightMaskFractionMax = Field(dtype=float, default=0.95,
                                  doc="Maximum allowed fraction of masked pixes from the \"bright\" "
                                  "detection stage (to mask regions unsuitable for sky sourcess). "
                                  "If this fraction is exeeded, the detection threshold for this stage "
                                  "will be increased by bisectFactor until the fraction of masked "
                                  "pixels drops below this threshold.")
    bisectFactor = Field(dtype=float, default=1.2,
                         doc="Factor by which to increase thresholds in brightMaskFractionMax loop.")

    def setDefaults(self):
        SourceDetectionConfig.setDefaults(self)
        self.skyObjects.nSources = 1000  # For good statistics
        for maskStr in ["SAT"]:
            if maskStr not in self.skyObjects.avoidMask:
                self.skyObjects.avoidMask.append(maskStr)


class DynamicDetectionTask(SourceDetectionTask):
    """Detection of sources on an image with a dynamic threshold

    We first detect sources using a lower threshold than normal (see config
    parameter ``prelimThresholdFactor``) in order to identify good sky regions
    (configurable ``skyObjects``). Then we perform forced PSF photometry on
    those sky regions. Using those PSF flux measurements and estimated errors,
    we set the threshold so that the stdev of the measurements matches the
    median estimated error.

    Besides the usual initialisation of configurables, we also set up
    the forced measurement which is deliberately not represented in
    this Task's configuration parameters because we're using it as
    part of the algorithm and we don't want to allow it to be modified.
    """
    ConfigClass = DynamicDetectionConfig
    _DefaultName = "dynamicDetection"

    def __init__(self, *args, **kwargs):

        SourceDetectionTask.__init__(self, *args, **kwargs)
        self.makeSubtask("skyObjects")

        # Set up forced measurement.
        config = ForcedMeasurementTask.ConfigClass()
        config.plugins.names = ['base_TransformedCentroid', 'base_PsfFlux', 'base_LocalBackground']
        # We'll need the "centroid" and "psfFlux" slots
        for slot in ("shape", "psfShape", "apFlux", "modelFlux", "gaussianFlux", "calibFlux"):
            setattr(config.slots, slot, None)
        config.copyColumns = {}
        self.skySchema = SourceTable.makeMinimalSchema()
        self.skyMeasurement = ForcedMeasurementTask(config=config, name="skyMeasurement", parentTask=self,
                                                    refSchema=self.skySchema)

    def calculateThreshold(self, exposure, seed, sigma=None, minFractionSourcesFactor=1.0, isBgTweak=False):
        """Calculate new threshold

        This is the main functional addition to the vanilla
        `SourceDetectionTask`.

        We identify sky objects and perform forced PSF photometry on
        them. Using those PSF flux measurements and estimated errors,
        we set the threshold so that the stdev of the measurements
        matches the median estimated error.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure on which we're detecting sources.
        seed : `int`
            RNG seed to use for finding sky objects.
        sigma : `float`, optional
            Gaussian sigma of smoothing kernel; if not provided,
            will be deduced from the exposure's PSF.
        minFractionSourcesFactor : `float`
            Change the fraction of required sky sources from that set in
            ``self.config.minFractionSources`` by this factor.  NOTE: this
            is intended for use in the background tweak pass (the detection
            threshold is much lower there, so many more pixels end up marked
            as DETECTED or DETECTED_NEGATIVE, leaving less room for sky
            object placement).
        isBgTweak : `bool`
           Set to ``True`` for the background tweak pass (for more helpful
           log messages).

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with components:

            ``multiplicative``
                Multiplicative factor to be applied to the
                configured detection threshold (`float`).
            ``additive``
                Additive factor to be applied to the background
                level (`float`).

        Raises
        ------
        NoWorkFound
            Raised if the number of good sky sources found is less than the
            minimum fraction
            (``self.config.minFractionSources``*``minFractionSourcesFactor``)
            of the number requested (``self.skyObjects.config.nSources``).
        """
        wcsIsNone = exposure.getWcs() is None
        if wcsIsNone:  # create a dummy WCS as needed by ForcedMeasurementTask
            self.log.info("WCS for exposure is None.  Setting a dummy WCS for dynamic detection.")
            exposure.setWcs(makeSkyWcs(crpix=geom.Point2D(0, 0),
                                       crval=geom.SpherePoint(0, 0, geom.degrees),
                                       cdMatrix=makeCdMatrix(scale=1e-5*geom.degrees)))
        fp = self.skyObjects.run(exposure.maskedImage.mask, seed)
        skyFootprints = FootprintSet(exposure.getBBox())
        skyFootprints.setFootprints(fp)
        table = SourceTable.make(self.skyMeasurement.schema)
        catalog = SourceCatalog(table)
        catalog.reserve(len(skyFootprints.getFootprints()))
        skyFootprints.makeSources(catalog)
        key = catalog.getCentroidSlot().getMeasKey()
        for source in catalog:
            peaks = source.getFootprint().getPeaks()
            assert len(peaks) == 1
            source.set(key, peaks[0].getF())
            source.updateCoord(exposure.getWcs())

        # Forced photometry on sky objects
        self.skyMeasurement.run(catalog, exposure, catalog, exposure.getWcs())

        # Calculate new threshold
        fluxes = catalog["base_PsfFlux_instFlux"]
        area = catalog["base_PsfFlux_area"]
        bg = catalog["base_LocalBackground_instFlux"]

        good = (~catalog["base_PsfFlux_flag"] & ~catalog["base_LocalBackground_flag"]
                & np.isfinite(fluxes) & np.isfinite(area) & np.isfinite(bg))

        minNumSources = int(self.config.minFractionSources*self.skyObjects.config.nSources)
        # Reduce the number of sky sources required if requested, but ensure
        # a minumum of 3.
        if minFractionSourcesFactor != 1.0:
            minNumSources = max(3, int(minNumSources*minFractionSourcesFactor))
        if good.sum() < minNumSources:
            if not isBgTweak:
                msg = (f"Insufficient good sky source flux measurements ({good.sum()} < "
                       f"{minNumSources}) for dynamic threshold calculation.")
            else:
                msg = (f"Insufficient good sky source flux measurements ({good.sum()} < "
                       f"{minNumSources}) for background tweak calculation.")

            nPix = exposure.mask.array.size
            badPixelMask = lsst.afw.image.Mask.getPlaneBitMask(["NO_DATA", "BAD"])
            nGoodPix = np.sum(exposure.mask.array & badPixelMask == 0)
            if nGoodPix/nPix > 0.2:
                detectedPixelMask = lsst.afw.image.Mask.getPlaneBitMask(["DETECTED", "DETECTED_NEGATIVE"])
                nDetectedPix = np.sum(exposure.mask.array & detectedPixelMask != 0)
                msg += (f" However, {nGoodPix}/{nPix} pixels are not marked NO_DATA or BAD, "
                        "so there should be sufficient area to locate suitable sky sources. "
                        f"Note that {nDetectedPix} of {nGoodPix} \"good\" pixels were marked "
                        "as DETECTED or DETECTED_NEGATIVE.")
                raise RuntimeError(msg)
            raise NoWorkFound(msg)

        if not isBgTweak:
            self.log.info("Number of good sky sources used for dynamic detection: %d (of %d requested).",
                          good.sum(), self.skyObjects.config.nSources)
        else:
            self.log.info("Number of good sky sources used for dynamic detection background tweak:"
                          " %d (of %d requested).", good.sum(), self.skyObjects.config.nSources)
        bgMedian = np.median((fluxes/area)[good])

        lq, uq = np.percentile((fluxes - bg*area)[good], [25.0, 75.0])
        stdevMeas = 0.741*(uq - lq)
        medianError = np.median(catalog["base_PsfFlux_instFluxErr"][good])
        if wcsIsNone:
            exposure.setWcs(None)
        return Struct(multiplicative=medianError/stdevMeas, additive=bgMedian)

    def detectFootprints(self, exposure, doSmooth=True, sigma=None, clearMask=True, expId=None,
                         background=None):
        """Detect footprints with a dynamic threshold

        This varies from the vanilla ``detectFootprints`` method because we
        do detection three times: first with a high threshold to detect
        "bright" (both positive and negative, the latter to identify very
        over-subtracted regions) sources for which we grow the DETECTED and
        DETECTED_NEGATIVE masks significantly to account for wings.  Second,
        with a low threshold to mask all non-empty regions of the image. These
        two masks are combined and used to identify regions of sky
        uncontaminated by objects.  A final round of detection is then done
        with the new calculated threshold.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to process; DETECTED{,_NEGATIVE} mask plane will be
            set in-place.
        doSmooth : `bool`, optional
            If True, smooth the image before detection using a Gaussian
            of width ``sigma``.
        sigma : `float`, optional
            Gaussian Sigma of PSF (pixels); used for smoothing and to grow
            detections; if `None` then measure the sigma of the PSF of the
            ``exposure``.
        clearMask : `bool`, optional
            Clear both DETECTED and DETECTED_NEGATIVE planes before running
            detection.
        expId : `int`, optional
            Exposure identifier, used as a seed for the random number
            generator. If absent, the seed will be the sum of the image.
        background : `lsst.afw.math.BackgroundList`, optional
            Background that was already subtracted from the exposure; will be
            modified in-place if ``reEstimateBackground=True``.

        Returns
        -------
        resutls : `lsst.pipe.base.Struct`
            The results `~lsst.pipe.base.Struct` contains:

            ``positive``
                Positive polarity footprints.
                (`lsst.afw.detection.FootprintSet` or `None`)
            ``negative``
                Negative polarity footprints.
                (`lsst.afw.detection.FootprintSet` or `None`)
            ``numPos``
                Number of footprints in positive or 0 if detection polarity was
                negative. (`int`)
            ``numNeg``
                Number of footprints in negative or 0 if detection polarity was
                positive. (`int`)
            ``background``
                Re-estimated background.  `None` or the input ``background``
                if ``reEstimateBackground==False``.
                (`lsst.afw.math.BackgroundList`)
            ``factor``
                Multiplication factor applied to the configured detection
                threshold. (`float`)
            ``prelim``
                Results from preliminary detection pass.
                (`lsst.pipe.base.Struct`)
        """
        maskedImage = exposure.maskedImage

        if clearMask:
            self.clearMask(maskedImage.mask)
        else:
            oldDetected = maskedImage.mask.array & maskedImage.mask.getPlaneBitMask(["DETECTED",
                                                                                     "DETECTED_NEGATIVE"])
        nPix = maskedImage.mask.array.size
        badPixelMask = lsst.afw.image.Mask.getPlaneBitMask(["NO_DATA", "BAD"])
        nGoodPix = np.sum(maskedImage.mask.array & badPixelMask == 0)
        self.log.info("Number of good data pixels (i.e. not NO_DATA or BAD): {} ({:.1f}% of total)".
                      format(nGoodPix, 100*nGoodPix/nPix))

        with self.tempWideBackgroundContext(exposure):
            # Could potentially smooth with a wider kernel than the PSF in
            # order to better pick up the wings of stars and galaxies, but for
            # now sticking with the PSF as that's more simple.
            psf = self.getPsf(exposure, sigma=sigma)
            convolveResults = self.convolveImage(maskedImage, psf, doSmooth=doSmooth)

            if self.config.doBrightPrelimDetection:
                brightDetectedMask = self._computeBrightDetectionMask(maskedImage, convolveResults)

            middle = convolveResults.middle
            sigma = convolveResults.sigma
            prelim = self.applyThreshold(
                middle, maskedImage.getBBox(), factor=self.config.prelimThresholdFactor,
                factorNeg=self.config.prelimNegMultiplier*self.config.prelimThresholdFactor
            )
            self.finalizeFootprints(
                maskedImage.mask, prelim, sigma, factor=self.config.prelimThresholdFactor,
                factorNeg=self.config.prelimNegMultiplier*self.config.prelimThresholdFactor
            )
            if self.config.doBrightPrelimDetection:
                # Combine prelim and bright detection masks for multiplier
                # determination.
                maskedImage.mask.array |= brightDetectedMask

            # Calculate the proper threshold
            # seed needs to fit in a C++ 'int' so pybind doesn't choke on it
            seed = (expId if expId is not None else int(maskedImage.image.array.sum())) % (2**31 - 1)
            threshResults = self.calculateThreshold(exposure, seed, sigma=sigma)
            factor = threshResults.multiplicative
            self.log.info("Modifying configured detection threshold by factor %f to %f",
                          factor, factor*self.config.thresholdValue)

            # Blow away preliminary (low threshold) detection mask
            self.clearMask(maskedImage.mask)
            if not clearMask:
                maskedImage.mask.array |= oldDetected

            # Rinse and repeat thresholding with new calculated threshold
            results = self.applyThreshold(middle, maskedImage.getBBox(), factor)
            results.prelim = prelim
            results.background = background if background is not None else lsst.afw.math.BackgroundList()
            if self.config.doTempLocalBackground:
                self.applyTempLocalBackground(exposure, middle, results)
            self.finalizeFootprints(maskedImage.mask, results, sigma, factor=factor)

            self.clearUnwantedResults(maskedImage.mask, results)

        if self.config.reEstimateBackground:
            self.reEstimateBackground(maskedImage, results.background)

        self.display(exposure, results, middle)

        if self.config.doBackgroundTweak:
            # Re-do the background tweak after any temporary backgrounds have
            # been restored.
            #
            # But we want to keep any large-scale background (e.g., scattered
            # light from bright stars) from being selected for sky objects in
            # the calculation, so do another detection pass without either the
            # local or wide temporary background subtraction; the DETECTED
            # pixels will mark the area to ignore.
            originalMask = maskedImage.mask.array.copy()
            try:
                self.clearMask(exposure.mask)
                convolveResults = self.convolveImage(maskedImage, psf, doSmooth=doSmooth)
                tweakDetResults = self.applyThreshold(convolveResults.middle, maskedImage.getBBox(), factor)
                self.finalizeFootprints(maskedImage.mask, tweakDetResults, sigma, factor=factor)
                bgLevel = self.calculateThreshold(exposure, seed, sigma=sigma, minFractionSourcesFactor=0.5,
                                                  isBgTweak=True).additive
            finally:
                maskedImage.mask.array[:] = originalMask
            self.tweakBackground(exposure, bgLevel, results.background)

        return results

    def tweakBackground(self, exposure, bgLevel, bgList=None):
        """Modify the background by a constant value

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure for which to tweak background.
        bgLevel : `float`
            Background level to remove
        bgList : `lsst.afw.math.BackgroundList`, optional
            List of backgrounds to append to.

        Returns
        -------
        bg : `lsst.afw.math.BackgroundMI`
            Constant background model.
        """
        self.log.info("Tweaking background by %f to match sky photometry", bgLevel)
        exposure.image -= bgLevel
        bgStats = lsst.afw.image.MaskedImageF(1, 1)
        bgStats.set(bgLevel, 0, bgLevel)
        bg = lsst.afw.math.BackgroundMI(exposure.getBBox(), bgStats)
        bgData = (bg, lsst.afw.math.Interpolate.LINEAR, lsst.afw.math.REDUCE_INTERP_ORDER,
                  lsst.afw.math.ApproximateControl.UNKNOWN, 0, 0, False)
        if bgList is not None:
            bgList.append(bgData)
        return bg

    def _computeBrightDetectionMask(self, maskedImage, convolveResults):
        """Perform an initial bright source detection pass.

        Perform an initial bright object detection pass using a high detection
        threshold. The footprints in this pass are grown significantly more
        than is typical to account for wings around bright sources.  The
        negative polarity detections in this pass help in masking severely
        over-subtracted regions.

        A maximum fraction of masked pixel from this pass is ensured via
        the config ``brightMaskFractionMax``.  If the masked pixel fraction is
        above this value, the detection thresholds here are increased by
        ``bisectFactor`` in a while loop until the detected masked fraction
        falls below this value.

        Parameters
        ----------
        maskedImage : `lsst.afw.image.MaskedImage`
            Masked image on which to run the detection.
        convolveResults :  `lsst.pipe.base.Struct`
            The results of the self.convolveImage function with attributes:

            ``middle``
                Convolved image, without the edges
               (`lsst.afw.image.MaskedImage`).
            ``sigma``
                Gaussian sigma used for the convolution (`float`).

        Returns
        -------
        brightDetectedMask : `numpy.ndarray`
            Boolean array representing the union of the bright detection pass
            DETECTED and DETECTED_NEGATIVE masks.
        """
        # Initialize some parameters.
        brightPosFactor = (
            self.config.prelimThresholdFactor*self.config.brightMultiplier/self.config.bisectFactor
        )
        brightNegFactor = self.config.brightNegFactor/self.config.bisectFactor
        nPix = 1
        nPixDet = 1
        nPixDetNeg = 1
        brightMaskFractionMax = self.config.brightMaskFractionMax

        # Loop until masked fraction is smaller than
        # brightMaskFractionMax, increasing the thresholds by
        # config.bisectFactor on each iteration (rarely necessary
        # for current defaults).
        while nPixDetNeg/nPix > brightMaskFractionMax or nPixDet/nPix > brightMaskFractionMax:
            self.clearMask(maskedImage.mask)
            brightPosFactor *= self.config.bisectFactor
            brightNegFactor *= self.config.bisectFactor
            prelimBright = self.applyThreshold(convolveResults.middle, maskedImage.getBBox(),
                                               factor=brightPosFactor, factorNeg=brightNegFactor)
            self.finalizeFootprints(
                maskedImage.mask, prelimBright, convolveResults.sigma*self.config.brightGrowFactor,
                factor=brightPosFactor, factorNeg=brightNegFactor
            )
            # Check that not too many pixels got masked.
            nPix = maskedImage.mask.array.size
            nPixDet = countMaskedPixels(maskedImage, "DETECTED")
            self.log.info("Number (%) of bright DETECTED pix: {} ({:.1f}%)".
                          format(nPixDet, 100*nPixDet/nPix))
            nPixDetNeg = countMaskedPixels(maskedImage, "DETECTED_NEGATIVE")
            self.log.info("Number (%) of bright DETECTED_NEGATIVE pix: {} ({:.1f}%)".
                          format(nPixDetNeg, 100*nPixDetNeg/nPix))
            if nPixDetNeg/nPix > brightMaskFractionMax or nPixDet/nPix > brightMaskFractionMax:
                self.log.warn("Too high a fraction (%.1f > %.1f) of pixels were masked with current "
                              "\"bright\" detection round thresholds.  Increasing by a factor of %f "
                              "and trying again.", max(nPixDetNeg, nPixDet)/nPix,
                              brightMaskFractionMax, self.config.bisectFactor)

        # Save the mask planes from the "bright" detection round, then
        # clear them before moving on to the "prelim" detection phase.
        brightDetectedMask = (maskedImage.mask.array
                              & maskedImage.mask.getPlaneBitMask(["DETECTED", "DETECTED_NEGATIVE"]))
        self.clearMask(maskedImage.mask)
        return brightDetectedMask


def countMaskedPixels(maskedIm, maskPlane):
    """Count the number of pixels in a given mask plane.

    Parameters
    ----------
    maskedIm : `lsst.afw.image.MaskedImage`
        Masked image to examine.
    maskPlane : `str`
        Name of the mask plane to examine.

    Returns
    -------
    nPixMasked : `int`
        Number of pixels with ``maskPlane`` bit set.
    """
    maskBit = maskedIm.mask.getPlaneBitMask(maskPlane)
    nPixMasked = np.sum(np.bitwise_and(maskedIm.mask.array, maskBit))/maskBit
    return nPixMasked
