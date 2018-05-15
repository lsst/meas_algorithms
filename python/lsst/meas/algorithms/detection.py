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

import lsst.afw.display
import lsst.afw.detection as afwDet
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from .subtractBackground import SubtractBackgroundTask


class SourceDetectionConfig(pexConfig.Config):
    """!Configuration parameters for the SourceDetectionTask
    """
    minPixels = pexConfig.RangeField(
        doc="detected sources with fewer than the specified number of pixels will be ignored",
        dtype=int, optional=False, default=1, min=0,
    )
    isotropicGrow = pexConfig.Field(
        doc="Pixels should be grown as isotropically as possible (slower)",
        dtype=bool, optional=False, default=False,
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
        doc="Threshold for footprints",
        dtype=float, optional=False, default=5.0, min=0.0,
    )
    includeThresholdMultiplier = pexConfig.RangeField(
        doc="Include threshold relative to thresholdValue",
        dtype=float, default=1.0, min=0.0,
    )
    thresholdType = pexConfig.ChoiceField(
        doc="specifies the desired flavor of Threshold",
        dtype=str, optional=False, default="stdev",
        allowed={
            "variance": "threshold applied to image variance",
            "stdev": "threshold applied to image std deviation",
            "value": "threshold applied to image value",
            "pixel_stdev": "threshold applied to per-pixel std deviation",
        },
    )
    thresholdPolarity = pexConfig.ChoiceField(
        doc="specifies whether to detect positive, or negative sources, or both",
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

## @addtogroup LSST_task_documentation
## @{
## @page sourceDetectionTask
## @ref SourceDetectionTask_ "SourceDetectionTask"
## @copybrief SourceDetectionTask
## @}


class SourceDetectionTask(pipeBase.Task):
    """!
@anchor SourceDetectionTask_

@brief Detect positive and negative sources on an exposure and return a new @link table.SourceCatalog@endlink.

@section meas_algorithms_detection_Contents Contents

 - @ref meas_algorithms_detection_Purpose
 - @ref meas_algorithms_detection_Initialize
 - @ref meas_algorithms_detection_Invoke
 - @ref meas_algorithms_detection_Config
 - @ref meas_algorithms_detection_Debug
 - @ref meas_algorithms_detection_Example

@section meas_algorithms_detection_Purpose      Description

@copybrief SourceDetectionTask

@section meas_algorithms_detection_Initialize   Task initialisation

@copydoc \_\_init\_\_

@section meas_algorithms_detection_Invoke       Invoking the Task

@copydoc run

@section meas_algorithms_detection_Config       Configuration parameters

See @ref SourceDetectionConfig

@section meas_algorithms_detection_Debug                Debug variables

The @link lsst.pipe.base.cmdLineTask.CmdLineTask command line task@endlink interface supports a
flag @c -d to import @b debug.py from your @c PYTHONPATH; see @ref baseDebug for more about @b debug.py files.

The available variables in SourceDetectionTask are:
<DL>
  <DT> @c display
  <DD>
  - If True, display the exposure on ds9's frame 0.  +ve detections in blue, -ve detections in cyan
  - If display > 1, display the convolved exposure on frame 1
</DL>

@section meas_algorithms_detection_Example      A complete example of using SourceDetectionTask

This code is in @link measAlgTasks.py@endlink in the examples directory, and can be run as @em e.g.
@code
examples/measAlgTasks.py --ds9
@endcode
@dontinclude measAlgTasks.py
The example also runs the SourceMeasurementTask; see @ref meas_algorithms_measurement_Example for more
explanation.

Import the task (there are some other standard imports; read the file if you're confused)
@skipline SourceDetectionTask

We need to create our task before processing any data as the task constructor
can add an extra column to the schema, but first we need an almost-empty Schema
@skipline makeMinimalSchema
after which we can call the constructor:
@skip SourceDetectionTask.ConfigClass
@until detectionTask

We're now ready to process the data (we could loop over multiple exposures/catalogues using the same
task objects).  First create the output table:
@skipline afwTable

And process the image
@skipline result
(You may not be happy that the threshold was set in the config before creating the Task rather than being set
separately for each exposure.  You @em can reset it just before calling the run method if you must, but we
should really implement a better solution).

We can then unpack and use the results:
@skip sources
@until print

<HR>
To investigate the @ref meas_algorithms_detection_Debug, put something like
@code{.py}
    import lsstDebug
    def DebugInfo(name):
        di = lsstDebug.getInfo(name)        # N.b. lsstDebug.Info(name) would call us recursively
        if name == "lsst.meas.algorithms.detection":
            di.display = 1

        return di

    lsstDebug.Info = DebugInfo
@endcode
into your debug.py file and run measAlgTasks.py with the @c --debug flag.
    """
    ConfigClass = SourceDetectionConfig
    _DefaultName = "sourceDetection"

    def __init__(self, schema=None, **kwds):
        """!Create the detection task.  Most arguments are simply passed onto pipe.base.Task.

        @param schema An lsst::afw::table::Schema used to create the output lsst.afw.table.SourceCatalog
        @param **kwds Keyword arguments passed to lsst.pipe.base.task.Task.__init__.

        If schema is not None and configured for 'both' detections,
        a 'flags.negative' field will be added to label detections made with a
        negative threshold.

        @note This task can add fields to the schema, so any code calling this task must ensure that
        these columns are indeed present in the input match list; see @ref Example
        """
        pipeBase.Task.__init__(self, **kwds)
        if schema is not None and self.config.thresholdPolarity == "both":
            self.negativeFlagKey = schema.addField(
                "flags_negative", type="Flag",
                doc="set if source was detected as significantly negative"
            )
        else:
            if self.config.thresholdPolarity == "both":
                self.log.warn("Detection polarity set to 'both', but no flag will be "
                              "set to distinguish between positive and negative detections")
            self.negativeFlagKey = None
        if self.config.reEstimateBackground:
            self.makeSubtask("background")
        if self.config.doTempLocalBackground:
            self.makeSubtask("tempLocalBackground")
        if self.config.doTempWideBackground:
            self.makeSubtask("tempWideBackground")

    @pipeBase.timeMethod
    def run(self, table, exposure, doSmooth=True, sigma=None, clearMask=True, expId=None):
        """!Run source detection and create a SourceCatalog.

        @param table    lsst.afw.table.SourceTable object that will be used to create the SourceCatalog.
        @param exposure Exposure to process; DETECTED mask plane will be set in-place.
        @param doSmooth if True, smooth the image before detection using a Gaussian of width sigma
                        (default: True)
        @param sigma    sigma of PSF (pixels); used for smoothing and to grow detections;
            if None then measure the sigma of the PSF of the exposure (default: None)
        @param clearMask Clear DETECTED{,_NEGATIVE} planes before running detection (default: True)
        @param expId    Exposure identifier (integer); unused by this implementation, but used for
            RNG seed by subclasses.

        @return a lsst.pipe.base.Struct with:
          - sources -- an lsst.afw.table.SourceCatalog object
          - fpSets --- lsst.pipe.base.Struct returned by @link detectFootprints @endlink

        @throws ValueError if flags.negative is needed, but isn't in table's schema
        @throws lsst.pipe.base.TaskError if sigma=None, doSmooth=True and the exposure has no PSF

        @note
        If you want to avoid dealing with Sources and Tables, you can use detectFootprints()
        to just get the afw::detection::FootprintSet%s.
        """
        if self.negativeFlagKey is not None and self.negativeFlagKey not in table.getSchema():
            raise ValueError("Table has incorrect Schema")
        results = self.detectFootprints(exposure=exposure, doSmooth=doSmooth, sigma=sigma,
                                        clearMask=clearMask, expId=expId)
        sources = afwTable.SourceCatalog(table)
        sources.reserve(results.numPos + results.numNeg)
        if results.negative:
            results.negative.makeSources(sources)
            if self.negativeFlagKey:
                for record in sources:
                    record.set(self.negativeFlagKey, True)
        if results.positive:
            results.positive.makeSources(sources)
        results.fpSets = results.copy()  # Backward compatibility
        results.sources = sources
        return results

    ## An alias for run             @deprecated Remove this alias after checking for where it's used
    makeSourceCatalog = run

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

        disp0 = lsst.afw.display.Display(frame=0)
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
            disp1 = lsst.afw.display.Display(frame=1)
            disp1.mtv(convolvedImage, title="PSF smoothed")

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
        bgImage = bg.getImageF()
        middle -= bgImage.Factory(bgImage, middle.getBBox())
        thresholdPos = self.makeThreshold(middle, "positive")
        thresholdNeg = self.makeThreshold(middle, "negative")
        if self.config.thresholdPolarity != "negative":
            self.updatePeaks(results.positive, middle, thresholdPos)
        if self.config.thresholdPolarity != "positive":
            self.updatePeaks(results.negative, middle, thresholdNeg)

    def clearMask(self, mask):
        """Clear the DETECTED and DETECTED_NEGATIVE mask planes

        Removes any previous detection mask in preparation for a new
        detection pass.

        Parameters
        ----------
        mask : `lsst.afw.image.Mask`
            Mask to be cleared.
        """
        mask &= ~(mask.getPlaneBitMask("DETECTED") | mask.getPlaneBitMask("DETECTED_NEGATIVE"))

    def calculateKernelSize(self, sigma):
        """Calculate size of smoothing kernel

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
        """Retrieve the PSF for an exposure

        If ``sigma`` is provided, we make a ``GaussianPsf`` with that,
        otherwise use the one from the ``exposure``.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure from which to retrieve the PSF.
        sigma : `float`, optional
            Gaussian sigma to use if provided.

        Returns
        -------
        psf : `lsst.afw.detection.Psf`
            PSF to use for detection.
        """
        if sigma is None:
            psf = exposure.getPsf()
            if psf is None:
                raise RuntimeError("Unable to determine PSF to use for detection: no sigma provided")
            sigma = psf.computeShape().getDeterminantRadius()
        size = self.calculateKernelSize(sigma)
        psf = afwDet.GaussianPsf(size, size, sigma)
        return psf

    def convolveImage(self, maskedImage, psf, doSmooth=True):
        """Convolve the image with the PSF

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
            Actually do the convolution?

        Return Struct contents
        ----------------------
        middle : `lsst.afw.image.MaskedImage`
            Convolved image, without the edges.
        sigma : `float`
            Gaussian sigma used for the convolution.
        """
        self.metadata.set("doSmooth", doSmooth)
        sigma = psf.computeShape().getDeterminantRadius()
        self.metadata.set("sigma", sigma)

        if not doSmooth:
            middle = maskedImage.Factory(maskedImage)
            return pipeBase.Struct(middle=middle, sigma=sigma)

        # Smooth using a Gaussian (which is separable, hence fast) of width sigma
        # Make a SingleGaussian (separable) kernel with the 'sigma'
        kWidth = self.calculateKernelSize(sigma)
        self.metadata.set("smoothingKernelWidth", kWidth)
        gaussFunc = afwMath.GaussianFunction1D(sigma)
        gaussKernel = afwMath.SeparableKernel(kWidth, kWidth, gaussFunc, gaussFunc)

        convolvedImage = maskedImage.Factory(maskedImage.getBBox())

        afwMath.convolve(convolvedImage, maskedImage, gaussKernel, afwMath.ConvolutionControl())
        #
        # Only search psf-smoothed part of frame
        #
        goodBBox = gaussKernel.shrinkBBox(convolvedImage.getBBox())
        middle = convolvedImage.Factory(convolvedImage, goodBBox, afwImage.PARENT, False)
        #
        # Mark the parts of the image outside goodBBox as EDGE
        #
        self.setEdgeBits(maskedImage, goodBBox, maskedImage.getMask().getPlaneBitMask("EDGE"))

        return pipeBase.Struct(middle=middle, sigma=sigma)

    def applyThreshold(self, middle, bbox, factor=1.0):
        """Apply thresholds to the convolved image

        Identifies ``Footprint``s, both positive and negative.

        The threshold can be modified by the provided multiplication
        ``factor``.

        Parameters
        ----------
        middle : `lsst.afw.image.MaskedImage`
            Convolved image to threshold.
        bbox : `lsst.afw.geom.Box2I`
            Bounding box of unconvolved image.
        factor : `float`
            Multiplier for the configured threshold.

        Return Struct contents
        ----------------------
        positive : `lsst.afw.detection.FootprintSet` or `None`
            Positive detection footprints, if configured.
        negative : `lsst.afw.detection.FootprintSet` or `None`
            Negative detection footprints, if configured.
        factor : `float`
            Multiplier for the configured threshold.
        """
        results = pipeBase.Struct(positive=None, negative=None, factor=factor)
        # Detect the Footprints (peaks may be replaced if doTempLocalBackground)
        if self.config.reEstimateBackground or self.config.thresholdPolarity != "negative":
            threshold = self.makeThreshold(middle, "positive", factor=factor)
            results.positive = afwDet.FootprintSet(
                middle,
                threshold,
                "DETECTED",
                self.config.minPixels
            )
            results.positive.setRegion(bbox)
        if self.config.reEstimateBackground or self.config.thresholdPolarity != "positive":
            threshold = self.makeThreshold(middle, "negative", factor=factor)
            results.negative = afwDet.FootprintSet(
                middle,
                threshold,
                "DETECTED_NEGATIVE",
                self.config.minPixels
            )
            results.negative.setRegion(bbox)

        return results

    def finalizeFootprints(self, mask, results, sigma, factor=1.0):
        """Finalize the detected footprints

        Grows the footprints, sets the ``DETECTED`` and ``DETECTED_NEGATIVE``
        mask planes, and logs the results.

        ``numPos`` (number of positive footprints), ``numPosPeaks`` (number
        of positive peaks), ``numNeg`` (number of negative footprints),
        ``numNegPeaks`` (number of negative peaks) entries are added to the
        detection results.

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
            Multiplier for the configured threshold.
        """
        for polarity, maskName in (("positive", "DETECTED"), ("negative", "DETECTED_NEGATIVE")):
            fpSet = getattr(results, polarity)
            if fpSet is None:
                continue
            if self.config.nSigmaToGrow > 0:
                nGrow = int((self.config.nSigmaToGrow * sigma) + 0.5)
                self.metadata.set("nGrow", nGrow)
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

        self.log.info("Detected%s%s%s to %g %s" %
                      (positive, " and" if positive and negative else "", negative,
                       self.config.thresholdValue*self.config.includeThresholdMultiplier*factor,
                       "DN" if self.config.thresholdType == "value" else "sigma"))

    def reEstimateBackground(self, maskedImage, backgrounds):
        """Estimate the background after detection

        Parameters
        ----------
        maskedImage : `lsst.afw.image.MaskedImage`
            Image on which to estimate the background.
        backgrounds : `lsst.afw.math.BackgroundList`
            List of backgrounds; modified.

        Returns
        -------
        bg : `lsst.afw.math.backgroundMI`
            Empirical background model.
        """
        bg = self.background.fitBackground(maskedImage)
        if self.config.adjustBackground:
            self.log.warn("Fiddling the background by %g", self.config.adjustBackground)
            bg += self.config.adjustBackground
        self.log.info("Resubtracting the background after object detection")
        maskedImage -= bg.getImageF()
        backgrounds.append(bg)
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

    @pipeBase.timeMethod
    def detectFootprints(self, exposure, doSmooth=True, sigma=None, clearMask=True, expId=None):
        """Detect footprints.

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
        expId : `dict`, optional
            Exposure identifier; unused by this implementation, but used for
            RNG seed by subclasses.

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
        background : `lsst.afw.math.BackgroundList`
            Re-estimated background.  `None` if
            ``reEstimateBackground==False``.
        factor : `float`
            Multiplication factor applied to the configured detection
            threshold.
        """
        maskedImage = exposure.maskedImage

        if clearMask:
            self.clearMask(maskedImage.getMask())

        psf = self.getPsf(exposure, sigma=sigma)
        with self.tempWideBackgroundContext(exposure):
            convolveResults = self.convolveImage(maskedImage, psf, doSmooth=doSmooth)
            middle = convolveResults.middle
            sigma = convolveResults.sigma

            results = self.applyThreshold(middle, maskedImage.getBBox())
            results.background = afwMath.BackgroundList()
            if self.config.doTempLocalBackground:
                self.applyTempLocalBackground(exposure, middle, results)
            self.finalizeFootprints(maskedImage.mask, results, sigma)

            if self.config.reEstimateBackground:
                self.reEstimateBackground(maskedImage, results.background)

            self.clearUnwantedResults(maskedImage.getMask(), results)
            self.display(exposure, results, middle)

        return results

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
        goodBBox : `lsst.afw.geom.Box2I`
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
            edgeMask = msk.Factory(msk, afwGeom.BoxI(afwGeom.PointI(x0, y0),
                                                     afwGeom.ExtentI(w, h)), afwImage.LOCAL)
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
