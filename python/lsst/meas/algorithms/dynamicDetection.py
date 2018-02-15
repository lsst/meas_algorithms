from __future__ import absolute_import, division, print_function

__all__ = ["DynamicDetectionConfig", "DynamicDetectionTask"]

import numpy as np

from lsst.pex.config import Field, ConfigurableField
from lsst.pipe.base import Task

from .detection import SourceDetectionConfig, SourceDetectionTask
from .skyObjects import SkyObjectsTask

from lsst.afw.detection import FootprintSet
from lsst.afw.table import SourceCatalog, SourceTable, IdFactory
from lsst.meas.base import ForcedMeasurementTask


class DynamicDetectionConfig(SourceDetectionConfig):
    """Configuration for DynamicDetectionTask"""
    prelimThresholdFactor = Field(dtype=float, default=0.5,
                                  doc="Fraction of the threshold to use for first pass (to find sky objects)")
    skyObjects = ConfigurableField(target=SkyObjectsTask, doc="Generate sky objects")


class DynamicDetectionTask(SourceDetectionTask):
    """Detection of sources on an image with a dynamic threshold

    We first detect sources using a lower threshold than normal (see config
    parameter ``prelimThresholdFactor``) in order to identify good sky regions
    (configurable ``skyObjects``). Then we perform forced PSF photometry on
    those sky regions. Using those PSF flux measurements and estimated errors,
    we set the threshold so that the stdev of the measurements matches the
    median estimated error.
    """
    ConfigClass = DynamicDetectionConfig
    _DefaultName = "dynamicDetection"

    def __init__(self, *args, **kwargs):
        """Constructor

        Besides the usual initialisation of configurables, we also set up
        the forced measurement which is deliberately not represented in
        this Task's configuration parameters because we're using it as part
        of the algorithm and we don't want to allow it to be modified.
        """
        SourceDetectionTask.__init__(self, *args, **kwargs)
        self.makeSubtask("skyObjects")

        # Set up forced measurement.
        config = ForcedMeasurementTask.ConfigClass()
        config.plugins.names = ['base_TransformedCentroid', 'base_PsfFlux']
        # We'll need the "centroid" and "psfFlux" slots
        for slot in ("shape", "psfShape", "apFlux", "modelFlux", "instFlux", "calibFlux"):
            setattr(config.slots, slot, None)
        config.copyColumns = {}
        self.skySchema = SourceTable.makeMinimalSchema()
        self.skyMeasurement = ForcedMeasurementTask(config=config, name="skyMeasurement", parentTask=self,
                                                    refSchema=self.skySchema)

    def calculateThreshold(self, exposure, seed, sigma=None):
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

        Returns
        -------
        factor : `float`
            Multiplication factor to be applied to the configured detection
            threshold.
        """
        # Make a catalog of sky objects
        fp = self.skyObjects.run(exposure.maskedImage.mask, seed)
        skyFootprints = FootprintSet(exposure.getBBox())
        skyFootprints.setFootprints(fp)
        table = SourceTable.make(self.skyMeasurement.schema)
        catalog = SourceCatalog(table)
        table.preallocate(len(skyFootprints.getFootprints()))
        skyFootprints.makeSources(catalog)
        key = catalog.getCentroidKey()
        for source in catalog:
            peaks = source.getFootprint().getPeaks()
            assert len(peaks) == 1
            source.set(key, peaks[0].getF())
            source.updateCoord(exposure.getWcs())

        # Forced photometry on sky objects
        self.skyMeasurement.run(catalog, exposure, catalog, exposure.getWcs())

        # Calculate new threshold
        fluxes = catalog["base_PsfFlux_flux"]
        lq, uq = np.percentile(fluxes, [25.0, 75.0])
        stdev = 0.741*(uq - lq)
        errors = catalog["base_PsfFlux_fluxSigma"]
        median = np.median(errors)
        return median/stdev

    def detectFootprints(self, exposure, doSmooth=True, sigma=None, clearMask=True, expId=None):
        """Detect footprints with a dynamic threshold

        This varies from the vanilla ``detectFootprints`` method because we
        do detection twice: one with a low threshold so that we can find
        sky uncontaminated by objects, then one more with the new calculated
        threshold.

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

        Return Struct contents
        ----------------------
        positive : `lsst.afw.detection.FootprintSet`
            Positive polarity footprints (may be `None`)
        negative : `lsst.afw.detection.FootprintSet`
            Negative polarity footprints (may be `None`)
        numPos : `int`
            Number of footprints in positive or 0 if detection polarity was
            negative.
        numNeg : `int`
            Number of footprints in negative or 0 if detection polarity was
            positive.
        background : `lsst.afw.math.BackgroundMI`
            Re-estimated background.  `None` if
            ``reEstimateBackground==False``.
        factor : `float`
            Multiplication factor applied to the configured detection
            threshold.
        """
        maskedImage = exposure.maskedImage

        if clearMask:
            self.clearMask(maskedImage.mask)
        else:
            oldDetected = maskedImage.mask.array & maskedImage.mask.getPlaneBitMask(["DETECTED",
                                                                                     "DETECTED_NEGATIVE"])

        with self.tempWideBackgroundContext(exposure):
            # Could potentially smooth with a wider kernel than the PSF in order to better pick up the
            # wings of stars and galaxies, but for now sticking with the PSF as that's more simple.
            psf = self.getPsf(exposure, sigma=sigma)
            convolveResults = self.convolveImage(maskedImage, psf, doSmooth=doSmooth)
            middle = convolveResults.middle
            sigma = convolveResults.sigma

            prelim = self.applyThreshold(middle, maskedImage.getBBox(), self.config.prelimThresholdFactor)
            self.finalizeFootprints(maskedImage.mask, prelim, sigma, self.config.prelimThresholdFactor)

            # Calculate the proper threshold
            # seed needs to fit in a C++ 'int' so pybind doesn't choke on it
            seed = (expId if expId is not None else int(maskedImage.image.array.sum())) % (2**31 - 1)
            factor = self.calculateThreshold(exposure, seed, sigma=sigma)
            self.log.info("Modifying configured detection threshold by factor %f to %f",
                          factor, factor*self.config.thresholdValue)

            # Blow away preliminary (low threshold) detection mask
            self.clearMask(maskedImage.mask)
            if not clearMask:
                maskedImage.mask.array |= oldDetected

            # Rinse and repeat thresholding with new calculated threshold
            results = self.applyThreshold(middle, maskedImage.getBBox(), factor)
            results.prelim = prelim
            if self.config.doTempLocalBackground:
                self.applyTempLocalBackground(exposure, middle, results)
            self.finalizeFootprints(maskedImage.mask, results, sigma, factor)

            if self.config.reEstimateBackground:
                self.reEstimateBackground(maskedImage, prelim)

            self.clearUnwantedResults(maskedImage.mask, results)
            self.display(exposure, results, middle)

        return results
