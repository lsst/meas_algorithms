#
# LSST Data Management System
#
# Copyright 2008-2017  AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#

__all__ = ("SourceDetectionConfig", "SourceDetectionTask", "addExposures")

from contextlib import contextmanager

import numpy as np

import lsst.geom
import lsst.afw.display as afwDisplay
import lsst.afw.detection as afwDet
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.utils.timer import timeMethod
from .subtractBackground import SubtractBackgroundTask, backgroundFlatContext


class SourceDetectionConfig(pexConfig.Config):
    """Configuration parameters for the SourceDetectionTask
    """
    minPixels = pexConfig.RangeField(
        doc="detected sources with fewer than the specified number of pixels will be ignored",
        dtype=int, optional=False, default=1, min=0,
    )
    isotropicGrow = pexConfig.Field(
        doc="Grow pixels as isotropically as possible? If False, use a Manhattan metric instead.",
        dtype=bool, default=True,
    )
    combinedGrow = pexConfig.Field(
        doc="Grow all footprints at the same time? This allows disconnected footprints to merge.",
        dtype=bool, default=True,
    )
    nSigmaToGrow = pexConfig.Field(
        doc="Grow detections by nSigmaToGrow * [PSF RMS width]; if 0 then do not grow",
        dtype=float, default=2.4,  # 2.4 pixels/sigma is roughly one pixel/FWHM
    )
    returnOriginalFootprints = pexConfig.Field(
        doc="Grow detections to set the image mask bits, but return the original (not-grown) footprints",
        dtype=bool, optional=False, default=False,
    )
    thresholdValue = pexConfig.RangeField(
        doc="Threshold for detecting footprints; exact meaning and units depend on thresholdType.",
        dtype=float, optional=False, default=5.0, min=0.0,
    )
    includeThresholdMultiplier = pexConfig.RangeField(
        doc="Multiplier on thresholdValue for whether a source is included in the output catalog."
            " For example, thresholdValue=5, includeThresholdMultiplier=10, thresholdType='pixel_stdev' "
            "results in a catalog of sources at >50 sigma with the detection mask and footprints "
            "including pixels >5 sigma.",
        dtype=float, default=1.0, min=0.0,
    )
    thresholdType = pexConfig.ChoiceField(
        doc="Specifies the meaning of thresholdValue.",
        dtype=str, optional=False, default="pixel_stdev",
        allowed={
            "variance": "threshold applied to image variance",
            "stdev": "threshold applied to image std deviation",
            "value": "threshold applied to image value",
            "pixel_stdev": "threshold applied to per-pixel std deviation",
        },
    )
    thresholdPolarity = pexConfig.ChoiceField(
        doc="Specifies whether to detect positive, or negative sources, or both.",
        dtype=str, optional=False, default="positive",
        allowed={
            "positive": "detect only positive sources",
            "negative": "detect only negative sources",
            "both": "detect both positive and negative sources",
        },
    )
    adjustBackground = pexConfig.Field(
        dtype=float,
        doc="Fiddle factor to add to the background; debugging only",
        default=0.0,
    )
    reEstimateBackground = pexConfig.Field(
        dtype=bool,
        doc="Estimate the background again after final source detection?",
        default=True, optional=False,
    )
    doApplyFlatBackgroundRatio = pexConfig.Field(
        doc="Convert from a photometrically flat image to one suitable for background subtraction? "
            "Only used if reEstimateBackground is True."
            "If True, then a backgroundToPhotometricRatio must be supplied to the task run method.",
        dtype=bool,
        default=False,
    )
    background = pexConfig.ConfigurableField(
        doc="Background re-estimation; ignored if reEstimateBackground false",
        target=SubtractBackgroundTask,
    )
    tempLocalBackground = pexConfig.ConfigurableField(
        doc=("A local (small-scale), temporary background estimation step run between "
             "detecting above-threshold regions and detecting the peaks within "
             "them; used to avoid detecting spuerious peaks in the wings."),
        target=SubtractBackgroundTask,
    )
    doTempLocalBackground = pexConfig.Field(
        dtype=bool,
        doc="Enable temporary local background subtraction? (see tempLocalBackground)",
        default=True,
    )
    tempWideBackground = pexConfig.ConfigurableField(
        doc=("A wide (large-scale) background estimation and removal before footprint and peak detection. "
             "It is added back into the image after detection. The purpose is to suppress very large "
             "footprints (e.g., from large artifacts) that the deblender may choke on."),
        target=SubtractBackgroundTask,
    )
    doTempWideBackground = pexConfig.Field(
        dtype=bool,
        doc="Do temporary wide (large-scale) background subtraction before footprint detection?",
        default=False,
    )
    nPeaksMaxSimple = pexConfig.Field(
        dtype=int,
        doc=("The maximum number of peaks in a Footprint before trying to "
             "replace its peaks using the temporary local background"),
        default=1,
    )
    nSigmaForKernel = pexConfig.Field(
        dtype=float,
        doc=("Multiple of PSF RMS size to use for convolution kernel bounding box size; "
             "note that this is not a half-size. The size will be rounded up to the nearest odd integer"),
        default=7.0,
    )
    statsMask = pexConfig.ListField(
        dtype=str,
        doc="Mask planes to ignore when calculating statistics of image (for thresholdType=stdev)",
        default=['BAD', 'SAT', 'EDGE', 'NO_DATA'],
    )
    excludeMaskPlanes = lsst.pex.config.ListField(
        dtype=str,
        default=[],
        doc="Mask planes to exclude when detecting sources."
    )

    def setDefaults(self):
        self.tempLocalBackground.binSize = 64
        self.tempLocalBackground.algorithm = "AKIMA_SPLINE"
        self.tempLocalBackground.useApprox = False
        # Background subtraction to remove a large-scale background (e.g., scattered light); restored later.
        # Want to keep it from exceeding the deblender size limit of 1 Mpix, so half that is reasonable.
        self.tempWideBackground.binSize = 512
        self.tempWideBackground.algorithm = "AKIMA_SPLINE"
        self.tempWideBackground.useApprox = False
        # Ensure we can remove even bright scattered light that is DETECTED
        for maskPlane in ("DETECTED", "DETECTED_NEGATIVE"):
            if maskPlane in self.tempWideBackground.ignoredPixelMask:
                self.tempWideBackground.ignoredPixelMask.remove(maskPlane)


class SourceDetectionTask(pipeBase.Task):
    """Detect peaks and footprints of sources in an image.

    This task expects the image to have been background subtracted first.
    Running detection on images with a non-zero-centered background may result
    in a single source detected on the entire image containing thousands of
    peaks, or other pathological outputs.

    Parameters
    ----------
    schema : `lsst.afw.table.Schema`
        Schema object used to create the output `lsst.afw.table.SourceCatalog`
    **kwds
        Keyword arguments passed to `lsst.pipe.base.Task.__init__`

    If schema is not None and configured for 'both' detections,
    a 'flags.negative' field will be added to label detections made with a
    negative threshold.

    Notes
    -----
    This task convolves the image with a Gaussian approximation to the PSF,
    matched to the sigma of the input exposure, because this is separable and
    fast. The PSF would have to be very non-Gaussian or non-circular for this
    approximation to have a significant impact on the signal-to-noise of the
    detected sources.
    """
    ConfigClass = SourceDetectionConfig
    _DefaultName = "sourceDetection"

    def __init__(self, schema=None, **kwds):
        pipeBase.Task.__init__(self, **kwds)
        if schema is not None and self.config.thresholdPolarity == "both":
            self.negativeFlagKey = schema.addField(
                "is_negative", type="Flag",
                doc="set if source was detected as significantly negative"
            )
        else:
            if self.config.thresholdPolarity == "both":
                self.log.warning("Detection polarity set to 'both', but no flag will be "
                                 "set to distinguish between positive and negative detections")
            self.negativeFlagKey = None
        if self.config.reEstimateBackground:
            self.makeSubtask("background")
        if self.config.doTempLocalBackground:
            self.makeSubtask("tempLocalBackground")
        if self.config.doTempWideBackground:
            self.makeSubtask("tempWideBackground")

    @timeMethod
    def run(self, table, exposure, doSmooth=True, sigma=None, clearMask=True, expId=None,
            background=None, backgroundToPhotometricRatio=None):
        r"""Detect sources and return catalog(s) of detections.

        Parameters
        ----------
        table : `lsst.afw.table.SourceTable`
            Table object that will be used to create the SourceCatalog.
        exposure : `lsst.afw.image.Exposure`
            Exposure to process; DETECTED mask plane will be set in-place.
        doSmooth : `bool`, optional
            If True, smooth the image before detection using a Gaussian of width
            ``sigma``, or the measured PSF width. Set to False when running on
            e.g. a pre-convolved image, or a mask plane.
        sigma : `float`, optional
            Sigma of PSF (pixels); used for smoothing and to grow detections;
            if None then measure the sigma of the PSF of the exposure
        clearMask : `bool`, optional
            Clear DETECTED{,_NEGATIVE} planes before running detection.
        expId : `int`, optional
            Exposure identifier; unused by this implementation, but used for
            RNG seed by subclasses.
        background : `lsst.afw.math.BackgroundList`, optional
            Background that was already subtracted from the exposure; will be
            modified in-place if ``reEstimateBackground=True``.
        backgroundToPhotometricRatio : `lsst.afw.image.Image`, optional
            Image to convert photometric-flattened image to
            background-flattened image if ``reEstimateBackground=True`` and
            exposure has been photometric-flattened.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            The `~lsst.pipe.base.Struct` contains:

            ``sources``
                Detected sources on the exposure.
                (`lsst.afw.table.SourceCatalog`)
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
                Re-estimated background. `None` if
                ``reEstimateBackground==False``.
                (`lsst.afw.math.BackgroundList`)
            ``factor``
                Multiplication factor applied to the configured detection
                threshold. (`float`)

        Raises
        ------
        ValueError
            Raised if flags.negative is needed, but isn't in table's schema.
        lsst.pipe.base.TaskError
            Raised if sigma=None, doSmooth=True and the exposure has no PSF.

        Notes
        -----
        If you want to avoid dealing with Sources and Tables, you can use
        `detectFootprints()` to just get the
        `~lsst.afw.detection.FootprintSet`\s.
        """
        if self.negativeFlagKey is not None and self.negativeFlagKey not in table.getSchema():
            raise ValueError("Table has incorrect Schema")
        results = self.detectFootprints(exposure=exposure, doSmooth=doSmooth, sigma=sigma,
                                        clearMask=clearMask, expId=expId, background=background,
                                        backgroundToPhotometricRatio=backgroundToPhotometricRatio)
        sources = afwTable.SourceCatalog(table)
        sources.reserve(results.numPos + results.numNeg)
        if results.negative:
            results.negative.makeSources(sources)
            if self.negativeFlagKey:
                for record in sources:
                    record.set(self.negativeFlagKey, True)
        if results.positive:
            results.positive.makeSources(sources)
        results.sources = sources
        return results

    def display(self, exposure, results, convolvedImage=None):
        """Display detections if so configured

        Displays the ``exposure`` in frame 0, overlays the detection peaks.

        Requires that ``lsstDebug`` has been set up correctly, so that
        ``lsstDebug.Info("lsst.meas.algorithms.detection")`` evaluates `True`.

        If the ``convolvedImage`` is non-`None` and
        ``lsstDebug.Info("lsst.meas.algorithms.detection") > 1``, the
        ``convolvedImage`` will be displayed in frame 1.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to display, on which will be plotted the detections.
        results : `lsst.pipe.base.Struct`
            Results of the 'detectFootprints' method, containing positive and
            negative footprints (which contain the peak positions that we will
            plot). This is a `Struct` with ``positive`` and ``negative``
            elements that are of type `lsst.afw.detection.FootprintSet`.
        convolvedImage : `lsst.afw.image.Image`, optional
            Convolved image used for thresholding.
        """
        try:
            import lsstDebug
            display = lsstDebug.Info(__name__).display
        except ImportError:
            try:
                display
            except NameError:
                display = False
        if not display:
            return

        afwDisplay.setDefaultMaskTransparency(75)

        disp0 = afwDisplay.Display(frame=0)
        disp0.mtv(exposure, title="detection")

        def plotPeaks(fps, ctype):
            if fps is None:
                return
            with disp0.Buffering():
                for fp in fps.getFootprints():
                    for pp in fp.getPeaks():
                        disp0.dot("+", pp.getFx(), pp.getFy(), ctype=ctype)
        plotPeaks(results.positive, "yellow")
        plotPeaks(results.negative, "red")

        if convolvedImage and display > 1:
            disp1 = afwDisplay.Display(frame=1)
            disp1.mtv(convolvedImage, title="PSF smoothed")

        disp2 = afwDisplay.Display(frame=2)
        disp2.mtv(afwImage.ImageF(np.sqrt(exposure.variance.array)), title="stddev")

    def applyTempLocalBackground(self, exposure, middle, results):
        """Apply a temporary local background subtraction

        This temporary local background serves to suppress noise fluctuations
        in the wings of bright objects.

        Peaks in the footprints will be updated.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure for which to fit local background.
        middle : `lsst.afw.image.MaskedImage`
            Convolved image on which detection will be performed
            (typically smaller than ``exposure`` because the
            half-kernel has been removed around the edges).
        results : `lsst.pipe.base.Struct`
            Results of the 'detectFootprints' method, containing positive and
            negative footprints (which contain the peak positions that we will
            plot). This is a `Struct` with ``positive`` and ``negative``
            elements that are of type `lsst.afw.detection.FootprintSet`.
        """
        # Subtract the local background from the smoothed image. Since we
        # never use the smoothed again we don't need to worry about adding
        # it back in.
        bg = self.tempLocalBackground.fitBackground(exposure.getMaskedImage())
        bgImage = bg.getImageF(self.tempLocalBackground.config.algorithm,
                               self.tempLocalBackground.config.undersampleStyle)
        middle -= bgImage.Factory(bgImage, middle.getBBox())
        if self.config.thresholdPolarity != "negative":
            results.positiveThreshold = self.makeThreshold(middle, "positive")
            self.updatePeaks(results.positive, middle, results.positiveThreshold)
        if self.config.thresholdPolarity != "positive":
            results.negativeThreshold = self.makeThreshold(middle, "negative")
            self.updatePeaks(results.negative, middle, results.negativeThreshold)

    def clearMask(self, mask):
        """Clear the DETECTED and DETECTED_NEGATIVE mask planes.

        Removes any previous detection mask in preparation for a new
        detection pass.

        Parameters
        ----------
        mask : `lsst.afw.image.Mask`
            Mask to be cleared.
        """
        mask &= ~(mask.getPlaneBitMask("DETECTED") | mask.getPlaneBitMask("DETECTED_NEGATIVE"))

    def calculateKernelSize(self, sigma):
        """Calculate the size of the smoothing kernel.

        Uses the ``nSigmaForKernel`` configuration parameter. Note
        that that is the full width of the kernel bounding box
        (so a value of 7 means 3.5 sigma on either side of center).
        The value will be rounded up to the nearest odd integer.

        Parameters
        ----------
        sigma : `float`
            Gaussian sigma of smoothing kernel.

        Returns
        -------
        size : `int`
            Size of the smoothing kernel.
        """
        return (int(sigma * self.config.nSigmaForKernel + 0.5)//2)*2 + 1  # make sure it is odd

    def getPsf(self, exposure, sigma=None):
        """Create a single Gaussian PSF for an exposure.

        If ``sigma`` is provided, we make a `~lsst.afw.detection.GaussianPsf`
        with that, otherwise use the sigma from the psf of the ``exposure`` to
        make the `~lsst.afw.detection.GaussianPsf`.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure from which to retrieve the PSF.
        sigma : `float`, optional
            Gaussian sigma to use if provided.

        Returns
        -------
        psf : `lsst.afw.detection.GaussianPsf`
            PSF to use for detection.

        Raises
        ------
        RuntimeError
            Raised if ``sigma`` is not provided and ``exposure`` does not
            contain a ``Psf`` object.
        """
        if sigma is None:
            psf = exposure.getPsf()
            if psf is None:
                raise RuntimeError("Unable to determine PSF to use for detection: no sigma provided")
            sigma = psf.computeShape(psf.getAveragePosition()).getDeterminantRadius()
        size = self.calculateKernelSize(sigma)
        psf = afwDet.GaussianPsf(size, size, sigma)
        return psf

    def convolveImage(self, maskedImage, psf, doSmooth=True):
        """Convolve the image with the PSF.

        We convolve the image with a Gaussian approximation to the PSF,
        because this is separable and therefore fast. It's technically a
        correlation rather than a convolution, but since we use a symmetric
        Gaussian there's no difference.

        The convolution can be disabled with ``doSmooth=False``. If we do
        convolve, we mask the edges as ``EDGE`` and return the convolved image
        with the edges removed. This is because we can't convolve the edges
        because the kernel would extend off the image.

        Parameters
        ----------
        maskedImage : `lsst.afw.image.MaskedImage`
            Image to convolve.
        psf : `lsst.afw.detection.Psf`
            PSF to convolve with (actually with a Gaussian approximation
            to it).
        doSmooth : `bool`
            Actually do the convolution? Set to False when running on
            e.g. a pre-convolved image, or a mask plane.

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            The `~lsst.pipe.base.Struct` contains:

            ``middle``
                Convolved image, without the edges. (`lsst.afw.image.MaskedImage`)
            ``sigma``
                Gaussian sigma used for the convolution. (`float`)
        """
        self.metadata["doSmooth"] = doSmooth
        sigma = psf.computeShape(psf.getAveragePosition()).getDeterminantRadius()
        self.metadata["sigma"] = sigma

        if not doSmooth:
            middle = maskedImage.Factory(maskedImage, deep=True)
            return pipeBase.Struct(middle=middle, sigma=sigma)

        # Smooth using a Gaussian (which is separable, hence fast) of width sigma
        # Make a SingleGaussian (separable) kernel with the 'sigma'
        kWidth = self.calculateKernelSize(sigma)
        self.metadata["smoothingKernelWidth"] = kWidth
        gaussFunc = afwMath.GaussianFunction1D(sigma)
        gaussKernel = afwMath.SeparableKernel(kWidth, kWidth, gaussFunc, gaussFunc)

        convolvedImage = maskedImage.Factory(maskedImage.getBBox())

        afwMath.convolve(convolvedImage, maskedImage, gaussKernel, afwMath.ConvolutionControl())

        # Only search psf-smoothed part of frame
        goodBBox = gaussKernel.shrinkBBox(convolvedImage.getBBox())
        middle = convolvedImage.Factory(convolvedImage, goodBBox, afwImage.PARENT, False)

        # Mark the parts of the image outside goodBBox as EDGE
        self.setEdgeBits(maskedImage, goodBBox, maskedImage.getMask().getPlaneBitMask("EDGE"))

        return pipeBase.Struct(middle=middle, sigma=sigma)

    def applyThreshold(self, middle, bbox, factor=1.0, factorNeg=None):
        r"""Apply thresholds to the convolved image

        Identifies `~lsst.afw.detection.Footprint`\s, both positive and negative.
        The threshold can be modified by the provided multiplication
        ``factor``.

        Parameters
        ----------
        middle : `lsst.afw.image.MaskedImage`
            Convolved image to threshold.
        bbox : `lsst.geom.Box2I`
            Bounding box of unconvolved image.
        factor : `float`
            Multiplier for the configured threshold.
        factorNeg : `float` or `None`
            Multiplier for the configured threshold for negative detection polarity.
            If `None`, will be set equal to ``factor`` (i.e. equal to the factor used
            for positive detection polarity).

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            The `~lsst.pipe.base.Struct` contains:

            ``positive``
                Positive detection footprints, if configured.
                (`lsst.afw.detection.FootprintSet` or `None`)
            ``negative``
                Negative detection footprints, if configured.
                (`lsst.afw.detection.FootprintSet` or `None`)
            ``factor``
                Multiplier for the configured threshold.
                (`float`)
            ``factorNeg``
                Multiplier for the configured threshold for negative detection polarity.
                (`float`)
        """
        if factorNeg is None:
            factorNeg = factor
            self.log.info("Setting factor for negative detections equal to that for positive "
                          "detections: %f", factor)
        results = pipeBase.Struct(positive=None, negative=None, factor=factor, factorNeg=factorNeg,
                                  positiveThreshold=None, negativeThreshold=None)
        # Detect the Footprints (peaks may be replaced if doTempLocalBackground)
        if self.config.reEstimateBackground or self.config.thresholdPolarity != "negative":
            results.positiveThreshold = self.makeThreshold(middle, "positive", factor=factor)
            results.positive = afwDet.FootprintSet(
                middle,
                results.positiveThreshold,
                "DETECTED",
                self.config.minPixels
            )
            results.positive.setRegion(bbox)
        if self.config.reEstimateBackground or self.config.thresholdPolarity != "positive":
            results.negativeThreshold = self.makeThreshold(middle, "negative", factor=factorNeg)
            results.negative = afwDet.FootprintSet(
                middle,
                results.negativeThreshold,
                "DETECTED_NEGATIVE",
                self.config.minPixels
            )
            results.negative.setRegion(bbox)

        return results

    def finalizeFootprints(self, mask, results, sigma, factor=1.0, factorNeg=None):
        """Finalize the detected footprints.

        Grow the footprints, set the ``DETECTED`` and ``DETECTED_NEGATIVE``
        mask planes, and log the results.

        ``numPos`` (number of positive footprints), ``numPosPeaks`` (number
        of positive peaks), ``numNeg`` (number of negative footprints),
        ``numNegPeaks`` (number of negative peaks) entries are added to the
        ``results`` struct.

        Parameters
        ----------
        mask : `lsst.afw.image.Mask`
            Mask image on which to flag detected pixels.
        results : `lsst.pipe.base.Struct`
            Struct of detection results, including ``positive`` and
            ``negative`` entries; modified.
        sigma : `float`
            Gaussian sigma of PSF.
        factor : `float`
            Multiplier for the configured threshold. Note that this is only
            used here for logging purposes.
        factorNeg : `float` or `None`
            Multiplier used for the negative detection polarity threshold.
            If `None`, a factor equal to ``factor`` (i.e. equal to the one used
            for positive detection polarity) is assumed. Note that this is only
            used here for logging purposes.
        """
        factorNeg = factor if factorNeg is None else factorNeg
        for polarity, maskName in (("positive", "DETECTED"), ("negative", "DETECTED_NEGATIVE")):
            fpSet = getattr(results, polarity)
            if fpSet is None:
                continue
            if self.config.nSigmaToGrow > 0:
                nGrow = int((self.config.nSigmaToGrow * sigma) + 0.5)
                self.metadata["nGrow"] = nGrow
                if self.config.combinedGrow:
                    fpSet = afwDet.FootprintSet(fpSet, nGrow, self.config.isotropicGrow)
                else:
                    stencil = (afwGeom.Stencil.CIRCLE if self.config.isotropicGrow else
                               afwGeom.Stencil.MANHATTAN)
                    for fp in fpSet:
                        fp.dilate(nGrow, stencil)
            fpSet.setMask(mask, maskName)
            if not self.config.returnOriginalFootprints:
                setattr(results, polarity, fpSet)

        results.numPos = 0
        results.numPosPeaks = 0
        results.numNeg = 0
        results.numNegPeaks = 0
        positive = ""
        negative = ""

        if results.positive is not None:
            results.numPos = len(results.positive.getFootprints())
            results.numPosPeaks = sum(len(fp.getPeaks()) for fp in results.positive.getFootprints())
            positive = " %d positive peaks in %d footprints" % (results.numPosPeaks, results.numPos)
        if results.negative is not None:
            results.numNeg = len(results.negative.getFootprints())
            results.numNegPeaks = sum(len(fp.getPeaks()) for fp in results.negative.getFootprints())
            negative = " %d negative peaks in %d footprints" % (results.numNegPeaks, results.numNeg)

        self.log.info("Detected%s%s%s to %g +ve and %g -ve %s",
                      positive, " and" if positive and negative else "", negative,
                      self.config.thresholdValue*self.config.includeThresholdMultiplier*factor,
                      self.config.thresholdValue*self.config.includeThresholdMultiplier*factorNeg,
                      "DN" if self.config.thresholdType == "value" else "sigma")

    def reEstimateBackground(self, maskedImage, backgrounds, backgroundToPhotometricRatio=None):
        """Estimate the background after detection

        Parameters
        ----------
        maskedImage : `lsst.afw.image.MaskedImage`
            Image on which to estimate the background.
        backgrounds : `lsst.afw.math.BackgroundList`
            List of backgrounds; modified.
        backgroundToPhotometricRatio : `lsst.afw.image.Image`, optional
            Image to convert photometric-flattened image to
            background-flattened image.

        Returns
        -------
        bg : `lsst.afw.math.backgroundMI`
            Empirical background model.
        """
        with backgroundFlatContext(
            maskedImage,
            self.config.doApplyFlatBackgroundRatio,
            backgroundToPhotometricRatio=backgroundToPhotometricRatio,
        ):
            bg = self.background.fitBackground(maskedImage)
            if self.config.adjustBackground:
                self.log.warning("Fiddling the background by %g", self.config.adjustBackground)
                bg += self.config.adjustBackground
            self.log.info("Resubtracting the background after object detection")
            maskedImage -= bg.getImageF(self.background.config.algorithm,
                                        self.background.config.undersampleStyle)

            actrl = bg.getBackgroundControl().getApproximateControl()
            backgrounds.append((bg, getattr(afwMath.Interpolate, self.background.config.algorithm),
                                bg.getAsUsedUndersampleStyle(), actrl.getStyle(), actrl.getOrderX(),
                                actrl.getOrderY(), actrl.getWeighting()))
        return bg

    def clearUnwantedResults(self, mask, results):
        """Clear unwanted results from the Struct of results

        If we specifically want only positive or only negative detections,
        drop the ones we don't want, and its associated mask plane.

        Parameters
        ----------
        mask : `lsst.afw.image.Mask`
            Mask image.
        results : `lsst.pipe.base.Struct`
            Detection results, with ``positive`` and ``negative`` elements;
            modified.
        """
        if self.config.thresholdPolarity == "positive":
            if self.config.reEstimateBackground:
                mask &= ~mask.getPlaneBitMask("DETECTED_NEGATIVE")
            results.negative = None
        elif self.config.thresholdPolarity == "negative":
            if self.config.reEstimateBackground:
                mask &= ~mask.getPlaneBitMask("DETECTED")
            results.positive = None

    @timeMethod
    def detectFootprints(self, exposure, doSmooth=True, sigma=None, clearMask=True, expId=None,
                         background=None, backgroundToPhotometricRatio=None):
        """Detect footprints on an exposure.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to process; DETECTED{,_NEGATIVE} mask plane will be
            set in-place.
        doSmooth : `bool`, optional
            If True, smooth the image before detection using a Gaussian
            of width ``sigma``, or the measured PSF width of ``exposure``.
            Set to False when running on e.g. a pre-convolved image, or a mask
            plane.
        sigma : `float`, optional
            Gaussian Sigma of PSF (pixels); used for smoothing and to grow
            detections; if `None` then measure the sigma of the PSF of the
            ``exposure``.
        clearMask : `bool`, optional
            Clear both DETECTED and DETECTED_NEGATIVE planes before running
            detection.
        expId : `dict`, optional
            Exposure identifier; unused by this implementation, but used for
            RNG seed by subclasses.
        background : `lsst.afw.math.BackgroundList`, optional
            Background that was already subtracted from the exposure; will be
            modified in-place if ``reEstimateBackground=True``.
        backgroundToPhotometricRatio : `lsst.afw.image.Image`, optional
            Image to convert photometric-flattened image to
            background-flattened image if ``reEstimateBackground=True`` and
            exposure has been photometric-flattened.

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
        """
        maskedImage = exposure.maskedImage

        if clearMask:
            self.clearMask(maskedImage.getMask())

        psf = self.getPsf(exposure, sigma=sigma)
        with self.tempWideBackgroundContext(exposure):
            convolveResults = self.convolveImage(maskedImage, psf, doSmooth=doSmooth)
            middle = convolveResults.middle
            sigma = convolveResults.sigma
            self.removeBadPixels(middle)

            results = self.applyThreshold(middle, maskedImage.getBBox())
            results.background = background if background is not None else afwMath.BackgroundList()

            if self.config.doTempLocalBackground:
                self.applyTempLocalBackground(exposure, middle, results)
            self.finalizeFootprints(maskedImage.mask, results, sigma)

            # Compute the significance of peaks after the peaks have been
            # finalized and after local background correction/updatePeaks, so
            # that the significance represents the "final" detection S/N.
            results.positive = self.setPeakSignificance(middle, results.positive, results.positiveThreshold)
            results.negative = self.setPeakSignificance(middle, results.negative, results.negativeThreshold,
                                                        negative=True)

            if self.config.reEstimateBackground:
                self.reEstimateBackground(
                    maskedImage,
                    results.background,
                    backgroundToPhotometricRatio=backgroundToPhotometricRatio,
                )

            self.clearUnwantedResults(maskedImage.getMask(), results)

            self.display(exposure, results, middle)

        return results

    def removeBadPixels(self, middle):
        """Set the significance of flagged pixels to zero.

        Parameters
        ----------
        middle : `lsst.afw.image.Exposure`
            Score or maximum likelihood difference image.
            The image plane will be modified in place.
        """
        badPixelMask = middle.mask.getPlaneBitMask(self.config.excludeMaskPlanes)
        badPixels = middle.mask.array & badPixelMask > 0
        middle.image.array[badPixels] = 0

    def setPeakSignificance(self, exposure, footprints, threshold, negative=False):
        """Set the significance of each detected peak to the pixel value divided
        by the appropriate standard-deviation for ``config.thresholdType``.

        Only sets significance for "stdev" and "pixel_stdev" thresholdTypes;
        we leave it undefined for "value" and "variance" as it does not have a
        well-defined meaning in those cases.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure that footprints were detected on, likely the convolved,
            local background-subtracted image.
        footprints : `lsst.afw.detection.FootprintSet`
            Footprints detected on the image.
        threshold : `lsst.afw.detection.Threshold`
            Threshold used to find footprints.
        negative : `bool`, optional
            Are we calculating for negative sources?
        """
        if footprints is None or footprints.getFootprints() == []:
            return footprints
        polarity = -1 if negative else 1

        # All incoming footprints have the same schema.
        mapper = afwTable.SchemaMapper(footprints.getFootprints()[0].peaks.schema)
        mapper.addMinimalSchema(footprints.getFootprints()[0].peaks.schema)
        mapper.addOutputField("significance", type=float,
                              doc="Ratio of peak value to configured standard deviation.")

        # Copy the old peaks to the new ones with a significance field.
        # Do this independent of the threshold type, so we always have a
        # significance field.
        newFootprints = afwDet.FootprintSet(footprints)
        for old, new in zip(footprints.getFootprints(), newFootprints.getFootprints()):
            newPeaks = afwDet.PeakCatalog(mapper.getOutputSchema())
            newPeaks.extend(old.peaks, mapper=mapper)
            new.getPeaks().clear()
            new.setPeakCatalog(newPeaks)

        # Compute the significance values.
        if self.config.thresholdType == "pixel_stdev":
            for footprint in newFootprints.getFootprints():
                footprint.updatePeakSignificance(exposure.variance, polarity)
        elif self.config.thresholdType == "stdev":
            sigma = threshold.getValue() / self.config.thresholdValue
            for footprint in newFootprints.getFootprints():
                footprint.updatePeakSignificance(polarity*sigma)
        else:
            for footprint in newFootprints.getFootprints():
                for peak in footprint.peaks:
                    peak["significance"] = 0

        return newFootprints

    def makeThreshold(self, image, thresholdParity, factor=1.0):
        """Make an afw.detection.Threshold object corresponding to the task's
        configuration and the statistics of the given image.

        Parameters
        ----------
        image : `afw.image.MaskedImage`
            Image to measure noise statistics from if needed.
        thresholdParity: `str`
            One of "positive" or "negative", to set the kind of fluctuations
            the Threshold will detect.
        factor : `float`
            Factor by which to multiply the configured detection threshold.
            This is useful for tweaking the detection threshold slightly.

        Returns
        -------
        threshold : `lsst.afw.detection.Threshold`
            Detection threshold.
        """
        parity = False if thresholdParity == "negative" else True
        thresholdValue = self.config.thresholdValue
        thresholdType = self.config.thresholdType
        if self.config.thresholdType == 'stdev':
            bad = image.getMask().getPlaneBitMask(self.config.statsMask)
            sctrl = afwMath.StatisticsControl()
            sctrl.setAndMask(bad)
            stats = afwMath.makeStatistics(image, afwMath.STDEVCLIP, sctrl)
            thresholdValue *= stats.getValue(afwMath.STDEVCLIP)
            thresholdType = 'value'

        threshold = afwDet.createThreshold(thresholdValue*factor, thresholdType, parity)
        threshold.setIncludeMultiplier(self.config.includeThresholdMultiplier)
        self.log.debug("Detection threshold: %s", threshold)
        return threshold

    def updatePeaks(self, fpSet, image, threshold):
        """Update the Peaks in a FootprintSet by detecting new Footprints and
        Peaks in an image and using the new Peaks instead of the old ones.

        Parameters
        ----------
        fpSet : `afw.detection.FootprintSet`
            Set of Footprints whose Peaks should be updated.
        image : `afw.image.MaskedImage`
            Image to detect new Footprints and Peak in.
        threshold : `afw.detection.Threshold`
            Threshold object for detection.

        Input Footprints with fewer Peaks than self.config.nPeaksMaxSimple
        are not modified, and if no new Peaks are detected in an input
        Footprint, the brightest original Peak in that Footprint is kept.
        """
        for footprint in fpSet.getFootprints():
            oldPeaks = footprint.getPeaks()
            if len(oldPeaks) <= self.config.nPeaksMaxSimple:
                continue
            # We detect a new FootprintSet within each non-simple Footprint's
            # bbox to avoid a big O(N^2) comparison between the two sets of
            # Footprints.
            sub = image.Factory(image, footprint.getBBox())
            fpSetForPeaks = afwDet.FootprintSet(
                sub,
                threshold,
                "",  # don't set a mask plane
                self.config.minPixels
            )
            newPeaks = afwDet.PeakCatalog(oldPeaks.getTable())
            for fpForPeaks in fpSetForPeaks.getFootprints():
                for peak in fpForPeaks.getPeaks():
                    if footprint.contains(peak.getI()):
                        newPeaks.append(peak)
            if len(newPeaks) > 0:
                del oldPeaks[:]
                oldPeaks.extend(newPeaks)
            else:
                del oldPeaks[1:]

    @staticmethod
    def setEdgeBits(maskedImage, goodBBox, edgeBitmask):
        """Set the edgeBitmask bits for all of maskedImage outside goodBBox

        Parameters
        ----------
        maskedImage : `lsst.afw.image.MaskedImage`
            Image on which to set edge bits in the mask.
        goodBBox : `lsst.geom.Box2I`
            Bounding box of good pixels, in ``LOCAL`` coordinates.
        edgeBitmask : `lsst.afw.image.MaskPixel`
            Bit mask to OR with the existing mask bits in the region
            outside ``goodBBox``.
        """
        msk = maskedImage.getMask()

        mx0, my0 = maskedImage.getXY0()
        for x0, y0, w, h in ([0, 0,
                              msk.getWidth(), goodBBox.getBeginY() - my0],
                             [0, goodBBox.getEndY() - my0, msk.getWidth(),
                              maskedImage.getHeight() - (goodBBox.getEndY() - my0)],
                             [0, 0,
                              goodBBox.getBeginX() - mx0, msk.getHeight()],
                             [goodBBox.getEndX() - mx0, 0,
                              maskedImage.getWidth() - (goodBBox.getEndX() - mx0), msk.getHeight()],
                             ):
            edgeMask = msk.Factory(msk, lsst.geom.BoxI(lsst.geom.PointI(x0, y0),
                                                       lsst.geom.ExtentI(w, h)), afwImage.LOCAL)
            edgeMask |= edgeBitmask

    @contextmanager
    def tempWideBackgroundContext(self, exposure):
        """Context manager for removing wide (large-scale) background

        Removing a wide (large-scale) background helps to suppress the
        detection of large footprints that may overwhelm the deblender.
        It does, however, set a limit on the maximum scale of objects.

        The background that we remove will be restored upon exit from
        the context manager.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure on which to remove large-scale background.

        Returns
        -------
        context : context manager
            Context manager that will ensure the temporary wide background
            is restored.
        """
        doTempWideBackground = self.config.doTempWideBackground
        if doTempWideBackground:
            self.log.info("Applying temporary wide background subtraction")
            original = exposure.maskedImage.image.array[:].copy()
            self.tempWideBackground.run(exposure).background
            # Remove NO_DATA regions (e.g., edge of the field-of-view); these can cause detections after
            # subtraction because of extrapolation of the background model into areas with no constraints.
            image = exposure.maskedImage.image
            mask = exposure.maskedImage.mask
            noData = mask.array & mask.getPlaneBitMask("NO_DATA") > 0
            isGood = mask.array & mask.getPlaneBitMask(self.config.statsMask) == 0
            image.array[noData] = np.median(image.array[~noData & isGood])
        try:
            yield
        finally:
            if doTempWideBackground:
                exposure.maskedImage.image.array[:] = original


def addExposures(exposureList):
    """Add a set of exposures together.

    Parameters
    ----------
    exposureList : `list` of `lsst.afw.image.Exposure`
        Sequence of exposures to add.

    Returns
    -------
    addedExposure : `lsst.afw.image.Exposure`
        An exposure of the same size as each exposure in ``exposureList``,
        with the metadata from ``exposureList[0]`` and a masked image equal
        to the sum of all the exposure's masked images.
    """
    exposure0 = exposureList[0]
    image0 = exposure0.getMaskedImage()

    addedImage = image0.Factory(image0, True)
    addedImage.setXY0(image0.getXY0())

    for exposure in exposureList[1:]:
        image = exposure.getMaskedImage()
        addedImage += image

    addedExposure = exposure0.Factory(addedImage, exposure0.getWcs())
    return addedExposure
