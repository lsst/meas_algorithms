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

import lsst.afw.detection as afwDet
import lsst.afw.display.ds9 as ds9
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
    nSigmaToGrow = pexConfig.Field(
        doc="Grow detections by nSigmaToGrow * sigma; if 0 then do not grow",
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
        doc=("A seperate background estimation and removal before footprint and peak detection. "
             "It is added back into the image after detection."),
        target=SubtractBackgroundTask,
    )
    doTempLocalBackground = pexConfig.Field(
        dtype=bool,
        doc="Do temporary interpolated background subtraction before footprint detection?",
        default=True,
    )
    nPeaksMaxSimple = pexConfig.Field(
        dtype=int,
        doc=("The maximum number of peaks in a Footprint before trying to "
             "replace its peaks using the temporary local background"),
        default=1,
    )

    def setDefaults(self):
        self.tempLocalBackground.binSize = 64
        self.tempLocalBackground.algorithm = "AKIMA_SPLINE"
        self.tempLocalBackground.useApprox = False

## \addtogroup LSST_task_documentation
## \{
## \page sourceDetectionTask
## \ref SourceDetectionTask_ "SourceDetectionTask"
## \copybrief SourceDetectionTask
## \}


class SourceDetectionTask(pipeBase.Task):
    """!
\anchor SourceDetectionTask_

\brief Detect positive and negative sources on an exposure and return a new \link table.SourceCatalog\endlink.

\section meas_algorithms_detection_Contents Contents

 - \ref meas_algorithms_detection_Purpose
 - \ref meas_algorithms_detection_Initialize
 - \ref meas_algorithms_detection_Invoke
 - \ref meas_algorithms_detection_Config
 - \ref meas_algorithms_detection_Debug
 - \ref meas_algorithms_detection_Example

\section meas_algorithms_detection_Purpose      Description

\copybrief SourceDetectionTask

\section meas_algorithms_detection_Initialize   Task initialisation

\copydoc \_\_init\_\_

\section meas_algorithms_detection_Invoke       Invoking the Task

\copydoc run

\section meas_algorithms_detection_Config       Configuration parameters

See \ref SourceDetectionConfig

\section meas_algorithms_detection_Debug                Debug variables

The \link lsst.pipe.base.cmdLineTask.CmdLineTask command line task\endlink interface supports a
flag \c -d to import \b debug.py from your \c PYTHONPATH; see \ref baseDebug for more about \b debug.py files.

The available variables in SourceDetectionTask are:
<DL>
  <DT> \c display
  <DD>
  - If True, display the exposure on ds9's frame 0.  +ve detections in blue, -ve detections in cyan
  - If display > 1, display the convolved exposure on frame 1
</DL>

\section meas_algorithms_detection_Example      A complete example of using SourceDetectionTask

This code is in \link measAlgTasks.py\endlink in the examples directory, and can be run as \em e.g.
\code
examples/measAlgTasks.py --ds9
\endcode
\dontinclude measAlgTasks.py
The example also runs the SourceMeasurementTask; see \ref meas_algorithms_measurement_Example for more
explanation.

Import the task (there are some other standard imports; read the file if you're confused)
\skipline SourceDetectionTask

We need to create our task before processing any data as the task constructor
can add an extra column to the schema, but first we need an almost-empty Schema
\skipline makeMinimalSchema
after which we can call the constructor:
\skip SourceDetectionTask.ConfigClass
@until detectionTask

We're now ready to process the data (we could loop over multiple exposures/catalogues using the same
task objects).  First create the output table:
\skipline afwTable

And process the image
\skipline result
(You may not be happy that the threshold was set in the config before creating the Task rather than being set
separately for each exposure.  You \em can reset it just before calling the run method if you must, but we
should really implement a better solution).

We can then unpack and use the results:
\skip sources
@until print

<HR>
To investigate the \ref meas_algorithms_detection_Debug, put something like
\code{.py}
    import lsstDebug
    def DebugInfo(name):
        di = lsstDebug.getInfo(name)        # N.b. lsstDebug.Info(name) would call us recursively
        if name == "lsst.meas.algorithms.detection":
            di.display = 1

        return di

    lsstDebug.Info = DebugInfo
\endcode
into your debug.py file and run measAlgTasks.py with the \c --debug flag.
    """
    ConfigClass = SourceDetectionConfig
    _DefaultName = "sourceDetection"

    def __init__(self, schema=None, **kwds):
        """!Create the detection task.  Most arguments are simply passed onto pipe.base.Task.

        \param schema An lsst::afw::table::Schema used to create the output lsst.afw.table.SourceCatalog
        \param **kwds Keyword arguments passed to lsst.pipe.base.task.Task.__init__.

        If schema is not None and configured for 'both' detections,
        a 'flags.negative' field will be added to label detections made with a
        negative threshold.

        \note This task can add fields to the schema, so any code calling this task must ensure that
        these columns are indeed present in the input match list; see \ref Example
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

    @pipeBase.timeMethod
    def run(self, table, exposure, doSmooth=True, sigma=None, clearMask=True):
        """!Run source detection and create a SourceCatalog.

        \param table    lsst.afw.table.SourceTable object that will be used to create the SourceCatalog.
        \param exposure Exposure to process; DETECTED mask plane will be set in-place.
        \param doSmooth if True, smooth the image before detection using a Gaussian of width sigma
                        (default: True)
        \param sigma    sigma of PSF (pixels); used for smoothing and to grow detections;
            if None then measure the sigma of the PSF of the exposure (default: None)
        \param clearMask Clear DETECTED{,_NEGATIVE} planes before running detection (default: True)

        \return a lsst.pipe.base.Struct with:
          - sources -- an lsst.afw.table.SourceCatalog object
          - fpSets --- lsst.pipe.base.Struct returned by \link detectFootprints \endlink

        \throws ValueError if flags.negative is needed, but isn't in table's schema
        \throws lsst.pipe.base.TaskError if sigma=None, doSmooth=True and the exposure has no PSF

        \note
        If you want to avoid dealing with Sources and Tables, you can use detectFootprints()
        to just get the afw::detection::FootprintSet%s.
        """
        if self.negativeFlagKey is not None and self.negativeFlagKey not in table.getSchema():
            raise ValueError("Table has incorrect Schema")
        fpSets = self.detectFootprints(exposure=exposure, doSmooth=doSmooth, sigma=sigma,
                                       clearMask=clearMask)
        sources = afwTable.SourceCatalog(table)
        table.preallocate(fpSets.numPos + fpSets.numNeg)  # not required, but nice
        if fpSets.negative:
            fpSets.negative.makeSources(sources)
            if self.negativeFlagKey:
                for record in sources:
                    record.set(self.negativeFlagKey, True)
        if fpSets.positive:
            fpSets.positive.makeSources(sources)
        return pipeBase.Struct(
            sources=sources,
            fpSets=fpSets
        )

    ## An alias for run             \deprecated Remove this alias after checking for where it's used
    makeSourceCatalog = run

    @pipeBase.timeMethod
    def detectFootprints(self, exposure, doSmooth=True, sigma=None, clearMask=True):
        """!Detect footprints.

        \param exposure Exposure to process; DETECTED{,_NEGATIVE} mask plane will be set in-place.
        \param doSmooth if True, smooth the image before detection using a Gaussian of width sigma
        \param sigma    sigma of PSF (pixels); used for smoothing and to grow detections;
            if None then measure the sigma of the PSF of the exposure
        \param clearMask Clear both DETECTED and DETECTED_NEGATIVE planes before running detection

        \return a lsst.pipe.base.Struct with fields:
        - positive: lsst.afw.detection.FootprintSet with positive polarity footprints (may be None)
        - negative: lsst.afw.detection.FootprintSet with negative polarity footprints (may be None)
        - numPos: number of footprints in positive or 0 if detection polarity was negative
        - numNeg: number of footprints in negative or 0 if detection polarity was positive
        - background: re-estimated background.  None if reEstimateBackground==False

        \throws lsst.pipe.base.TaskError if sigma=None and the exposure has no PSF
        """
        try:
            import lsstDebug
            display = lsstDebug.Info(__name__).display
        except ImportError:
            try:
                display
            except NameError:
                display = False

        if exposure is None:
            raise RuntimeError("No exposure for detection")

        maskedImage = exposure.getMaskedImage()
        region = maskedImage.getBBox()

        if clearMask:
            mask = maskedImage.getMask()
            mask &= ~(mask.getPlaneBitMask("DETECTED") | mask.getPlaneBitMask("DETECTED_NEGATIVE"))
            del mask

        if self.config.doTempLocalBackground:
            # Estimate the background, but add it back in instead of leaving
            # it subtracted (for now); we'll want to smooth before we
            # subtract it.
            tempBg = self.tempLocalBackground.fitBackground(
                exposure.getMaskedImage()
            )
            tempLocalBkgdImage = tempBg.getImageF()

        if sigma is None:
            psf = exposure.getPsf()
            if psf is None:
                raise pipeBase.TaskError("exposure has no PSF; must specify sigma")
            shape = psf.computeShape()
            sigma = shape.getDeterminantRadius()

        self.metadata.set("sigma", sigma)
        self.metadata.set("doSmooth", doSmooth)

        if not doSmooth:
            convolvedImage = maskedImage.Factory(maskedImage)
            middle = convolvedImage
        else:
            # smooth using a Gaussian (which is separate, hence fast) of width sigma
            # make a SingleGaussian (separable) kernel with the 'sigma'
            psf = exposure.getPsf()
            kWidth = (int(sigma * 7 + 0.5) // 2) * 2 + 1  # make sure it is odd
            self.metadata.set("smoothingKernelWidth", kWidth)
            gaussFunc = afwMath.GaussianFunction1D(sigma)
            gaussKernel = afwMath.SeparableKernel(kWidth, kWidth, gaussFunc, gaussFunc)

            convolvedImage = maskedImage.Factory(maskedImage.getBBox())

            afwMath.convolve(convolvedImage, maskedImage, gaussKernel, afwMath.ConvolutionControl())
            #
            # Only search psf-smooth part of frame
            #
            goodBBox = gaussKernel.shrinkBBox(convolvedImage.getBBox())
            middle = convolvedImage.Factory(convolvedImage, goodBBox, afwImage.PARENT, False)
            #
            # Mark the parts of the image outside goodBBox as EDGE
            #
            self.setEdgeBits(maskedImage, goodBBox, maskedImage.getMask().getPlaneBitMask("EDGE"))

        fpSets = pipeBase.Struct(positive=None, negative=None)

        # Detect the Footprints (peaks may be replaced if doTempLocalBackground)
        if self.config.thresholdPolarity != "negative":
            threshold = self.makeThreshold(middle, "positive")
            fpSets.positive = afwDet.FootprintSet(
                middle,
                threshold,
                "DETECTED",
                self.config.minPixels
            )
        if self.config.reEstimateBackground or self.config.thresholdPolarity != "positive":
            threshold = self.makeThreshold(middle, "negative")
            fpSets.negative = afwDet.FootprintSet(
                middle,
                threshold,
                "DETECTED_NEGATIVE",
                self.config.minPixels
            )

        if self.config.doTempLocalBackground:
            # Subtract the local background from the smoothed image. Since we
            # never use the smoothed again we don't need to worry about adding
            # it back in.
            tempLocalBkgdImage = tempLocalBkgdImage.Factory(tempLocalBkgdImage,
                                                            middle.getBBox())
            middle -= tempLocalBkgdImage
            thresholdPos = self.makeThreshold(middle, "positive")
            thresholdNeg = self.makeThreshold(middle, "negative")
            if self.config.thresholdPolarity != "negative":
                self.updatePeaks(fpSets.positive, middle, thresholdPos)
            if self.config.thresholdPolarity != "positive":
                self.updatePeaks(fpSets.negative, middle, thresholdNeg)

        for polarity, maskName in (("positive", "DETECTED"), ("negative", "DETECTED_NEGATIVE")):
            fpSet = getattr(fpSets, polarity)
            if fpSet is None:
                continue
            fpSet.setRegion(region)
            if self.config.nSigmaToGrow > 0:
                nGrow = int((self.config.nSigmaToGrow * sigma) + 0.5)
                self.metadata.set("nGrow", nGrow)
                fpSet = afwDet.FootprintSet(fpSet, nGrow, self.config.isotropicGrow)
            fpSet.setMask(maskedImage.getMask(), maskName)
            if not self.config.returnOriginalFootprints:
                setattr(fpSets, polarity, fpSet)

        fpSets.numPos = len(fpSets.positive.getFootprints()) if fpSets.positive is not None else 0
        fpSets.numNeg = len(fpSets.negative.getFootprints()) if fpSets.negative is not None else 0

        if self.config.thresholdPolarity != "negative":
            self.log.info("Detected %d positive sources to %g sigma.",
                          fpSets.numPos, self.config.thresholdValue*self.config.includeThresholdMultiplier)

        fpSets.background = None
        if self.config.reEstimateBackground:
            mi = exposure.getMaskedImage()
            bkgd = self.background.fitBackground(mi)

            if self.config.adjustBackground:
                self.log.warn("Fiddling the background by %g", self.config.adjustBackground)

                bkgd += self.config.adjustBackground
            fpSets.background = bkgd
            self.log.info("Resubtracting the background after object detection")

            mi -= bkgd.getImageF()
            del mi

        if self.config.thresholdPolarity == "positive":
            if self.config.reEstimateBackground:
                mask = maskedImage.getMask()
                mask &= ~mask.getPlaneBitMask("DETECTED_NEGATIVE")
                del mask
            fpSets.negative = None
        else:
            self.log.info("Detected %d negative sources to %g %s",
                          fpSets.numNeg, self.config.thresholdValue,
                          ("DN" if self.config.thresholdType == "value" else "sigma"))

        if display:
            ds9.mtv(exposure, frame=0, title="detection")
            x0, y0 = exposure.getXY0()

            def plotPeaks(fps, ctype):
                if fps is None:
                    return
                with ds9.Buffering():
                    for fp in fps.getFootprints():
                        for pp in fp.getPeaks():
                            ds9.dot("+", pp.getFx() - x0, pp.getFy() - y0, ctype=ctype)
            plotPeaks(fpSets.positive, "yellow")
            plotPeaks(fpSets.negative, "red")

            if convolvedImage and display and display > 1:
                ds9.mtv(convolvedImage, frame=1, title="PSF smoothed")

        return fpSets

    def makeThreshold(self, image, thresholdParity):
        """Make an afw.detection.Threshold object corresponding to the task's
        configuration and the statistics of the given image.

        Parameters
        ----------
        image : `afw.image.MaskedImage`
            Image to measure noise statistics from if needed.
        thresholdParity: `str`
            One of "positive" or "negative", to set the kind of fluctuations
            the Threshold will detect.
        """
        parity = False if thresholdParity == "negative" else True
        threshold = afwDet.createThreshold(self.config.thresholdValue,
                                           self.config.thresholdType, parity)
        threshold.setIncludeMultiplier(self.config.includeThresholdMultiplier)

        if self.config.thresholdType == 'stdev':
            bad = image.getMask().getPlaneBitMask(['BAD', 'SAT', 'EDGE',
                                                   'NO_DATA', ])
            sctrl = afwMath.StatisticsControl()
            sctrl.setAndMask(bad)
            stats = afwMath.makeStatistics(image, afwMath.STDEVCLIP, sctrl)
            thres = (stats.getValue(afwMath.STDEVCLIP) *
                     self.config.thresholdValue)
            threshold = afwDet.createThreshold(thres, 'value', parity)
            threshold.setIncludeMultiplier(
                self.config.includeThresholdMultiplier
            )

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
        """!Set the edgeBitmask bits for all of maskedImage outside goodBBox

        \param[in,out] maskedImage  image on which to set edge bits in the mask
        \param[in] goodBBox  bounding box of good pixels, in LOCAL coordinates
        \param[in] edgeBitmask  bit mask to OR with the existing mask bits in the region outside goodBBox
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


def addExposures(exposureList):
    """!Add a set of exposures together.

    \param[in] exposureList  sequence of exposures to add

    \return an exposure of the same size as each exposure in exposureList,
    with the metadata from exposureList[0] and a masked image equal to the
    sum of all the exposure's masked images.

    \throw LsstException if the exposures do not all have the same dimensions (but does not check xy0)
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
