# 
# LSST Data Management System
# Copyright 2008-2015 AURA/LSST.
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
import lsstDebug
import lsst.pex.logging as pexLogging

import lsst.pex.config as pexConfig
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.detection as afwDet
import lsst.pipe.base as pipeBase

__all__ = ("SourceDetectionConfig", "SourceDetectionTask", "getBackground",
           "estimateBackground", "BackgroundConfig", "addExposures")

import lsst.afw.display.ds9 as ds9

class BackgroundConfig(pexConfig.Config):
    """!Config for background estimation
    """
    statisticsProperty = pexConfig.ChoiceField(
        doc="type of statistic to use for grid points",
        dtype=str, default="MEANCLIP",
        allowed={
            "MEANCLIP": "clipped mean",
            "MEAN": "unclipped mean",
            "MEDIAN": "median",
            }
        )
    undersampleStyle = pexConfig.ChoiceField(
        doc="behaviour if there are too few points in grid for requested interpolation style",
        dtype=str, default="REDUCE_INTERP_ORDER",
        allowed={
            "THROW_EXCEPTION": "throw an exception if there are too few points",
            "REDUCE_INTERP_ORDER": "use an interpolation style with a lower order.",
            "INCREASE_NXNYSAMPLE": "Increase the number of samples used to make the interpolation grid.",
            }
        )
    binSize = pexConfig.RangeField(
        doc="how large a region of the sky should be used for each background point",
        dtype=int, default=256, min=10
        )
    algorithm = pexConfig.ChoiceField(
        doc="how to interpolate the background values. This maps to an enum; see afw::math::Background",
        dtype=str, default="NATURAL_SPLINE", optional=True,
        allowed={
            "CONSTANT" : "Use a single constant value",
            "LINEAR" : "Use linear interpolation",
            "NATURAL_SPLINE" : "cubic spline with zero second derivative at endpoints",
            "AKIMA_SPLINE": "higher-level nonlinear spline that is more robust to outliers",
            "NONE": "No background estimation is to be attempted",
            }
        )
    ignoredPixelMask = pexConfig.ListField(
        doc="Names of mask planes to ignore while estimating the background",
        dtype=str, default = ["BAD", "EDGE", "DETECTED", "DETECTED_NEGATIVE", "NO_DATA",],
        itemCheck = lambda x: x in afwImage.MaskU().getMaskPlaneDict().keys(),
        )
    isNanSafe = pexConfig.Field(
        doc="Ignore NaNs when estimating the background",
        dtype=bool, default=False,
    )

    useApprox = pexConfig.Field(
        doc="Use Approximate (Chebyshev) to model background.",
        dtype=bool, default=False,
    )
    approxOrderX = pexConfig.Field(
        doc="Approximation order in X for background Chebyshev (valid only with useApprox=True)",
        dtype=int, default=6,
    )
    # Note: Currently X- and Y-orders must be equal due to a limitation in math::Chebyshev1Function2
    # The following is being added so that the weighting attribute can also be configurable for the
    # call to afwMath.ApproximateControl
    approxOrderY = pexConfig.Field(
        doc="Approximation order in Y for background Chebyshev (valid only with useApprox=True)",
        dtype=int, default=-1,
    )
    weighting = pexConfig.Field(
        doc="Use inverse variance weighting in calculation (valid only with useApprox=True)",
        dtype=bool, default=True,
    )

    def validate(self):
        pexConfig.Config.validate(self)
        # Allow None to be used as an equivalent for "NONE", even though C++ expects the latter.
        if self.algorithm is None:
            self.algorithm = "NONE"

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
        dtype=float, default=2.4, # 2.4 pixels/sigma is roughly one pixel/FWHM
    )
    returnOriginalFootprints = pexConfig.Field(
        doc="Grow detections to set the image mask bits, but return the original (not-grown) footprints",
        dtype=bool, optional=False, default=False
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
        }
    )
    thresholdPolarity = pexConfig.ChoiceField(
        doc="specifies whether to detect positive, or negative sources, or both",
        dtype=str, optional=False, default="positive",
        allowed={
            "positive": "detect only positive sources",
            "negative": "detect only negative sources",
            "both": "detect both positive and negative sources",
        }
    )
    adjustBackground = pexConfig.Field(
        dtype = float,
        doc = "Fiddle factor to add to the background; debugging only",
        default=0.0,
    )
    reEstimateBackground = pexConfig.Field(
        dtype = bool,
        doc = "Estimate the background again after final source detection?",
        default=True, optional=False,
    )
    background = pexConfig.ConfigField(
        dtype=BackgroundConfig,
        doc="Background re-estimation configuration"
    )

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

\copydoc init

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
The example also runs the SourceMeasurementTask; see \ref meas_algorithms_measurement_Example for more explanation.

Import the task (there are some other standard imports; read the file if you're confused)
\skipline SourceDetectionTask

We need to create our task before processing any data as the task constructor
can add an extra column to the schema, but first we need an almost-empty Schema
\skipline makeMinimalSchema
after which we can call the constructor:
\skip SourceDetectionTask.ConfigClass
\until detectionTask

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
\until print

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

    # Need init as well as __init__ because "\copydoc __init__" fails (doxygen bug 732264)
    def init(self, schema=None, **kwds):
        """!Create the detection task.  Most arguments are simply passed onto pipe.base.Task.

        \param schema An lsst::afw::table::Schema used to create the output lsst.afw.table.SourceCatalog
        \param **kwds Keyword arguments passed to lsst.pipe.base.task.Task.__init__.

        If schema is not None, a 'flags.negative' field will be added to label detections
        made with a negative threshold.

        \note This task can add fields to the schema, so any code calling this task must ensure that
        these columns are indeed present in the input match list; see \ref Example
        """
        self.__init__(schema, **kwds)

    def __init__(self, schema=None, **kwds):
        """!Create the detection task.  See SourceDetectionTask.init for documentation
        """
        pipeBase.Task.__init__(self, **kwds)
        if schema is not None:
            self.negativeFlagKey = schema.addField(
                "flags_negative", type="Flag",
                doc="set if source was detected as significantly negative"
            )
        else:
            if self.config.thresholdPolarity == "both":
                self.log.log(self.log.WARN, "Detection polarity set to 'both', but no flag will be "\
                             "set to distinguish between positive and negative detections")
            self.negativeFlagKey = None

    @pipeBase.timeMethod
    def run(self, table, exposure, doSmooth=True, sigma=None, clearMask=True):
        """!Run source detection and create a SourceCatalog.

        \param table    lsst.afw.table.SourceTable object that will be used to create the SourceCatalog.
        \param exposure Exposure to process; DETECTED mask plane will be set in-place.
        \param doSmooth if True, smooth the image before detection using a Gaussian of width sigma (default: True)
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
        table.preallocate(fpSets.numPos + fpSets.numNeg) # not required, but nice
        if fpSets.negative:
            fpSets.negative.makeSources(sources)
            if self.negativeFlagKey:
                for record in sources:
                    record.set(self.negativeFlagKey, True)
        if fpSets.positive:
            fpSets.positive.makeSources(sources)
        return pipeBase.Struct(
            sources = sources,
            fpSets = fpSets
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
            kWidth = (int(sigma * 7 + 0.5) // 2) * 2 + 1 # make sure it is odd
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

        if self.config.thresholdPolarity != "negative":
            fpSets.positive = self.thresholdImage(middle, "positive")
        if self.config.reEstimateBackground or self.config.thresholdPolarity != "positive":
            fpSets.negative = self.thresholdImage(middle, "negative")

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
            self.log.log(self.log.INFO, "Detected %d positive sources to %g sigma." %
                         (fpSets.numPos, self.config.thresholdValue))

        fpSets.background = None
        if self.config.reEstimateBackground:
            mi = exposure.getMaskedImage()
            bkgd = getBackground(mi, self.config.background)

            if self.config.adjustBackground:
                self.log.log(self.log.WARN, "Fiddling the background by %g" % self.config.adjustBackground)

                bkgd += self.config.adjustBackground
            fpSets.background = bkgd
            self.log.log(self.log.INFO, "Resubtracting the background after object detection")

            mi -= bkgd.getImageF()
            del mi

        if self.config.thresholdPolarity == "positive":
            if self.config.reEstimateBackground:
                mask = maskedImage.getMask()
                mask &= ~mask.getPlaneBitMask("DETECTED_NEGATIVE")
                del mask
            fpSets.negative = None
        else:
            self.log.log(self.log.INFO, "Detected %d negative sources to %g %s" %
                         (fpSets.numNeg, self.config.thresholdValue,
                          ("DN" if self.config.thresholdType == "value" else "sigma")))

        if display:
            ds9.mtv(exposure, frame=0, title="detection")

            if convolvedImage and display and display > 1:
                ds9.mtv(convolvedImage, frame=1, title="PSF smoothed")

        return fpSets

    def thresholdImage(self, image, thresholdParity, maskName="DETECTED"):
        """!Threshold the convolved image, returning a FootprintSet.
        Helper function for detect().

        \param image The (optionally convolved) MaskedImage to threshold
        \param thresholdParity Parity of threshold
        \param maskName Name of mask to set

        \return FootprintSet
        """
        parity = False if thresholdParity == "negative" else True
        threshold = afwDet.createThreshold(self.config.thresholdValue, self.config.thresholdType, parity)
        threshold.setIncludeMultiplier(self.config.includeThresholdMultiplier)

        if self.config.thresholdType == 'stdev':
            bad = image.getMask().getPlaneBitMask(['BAD', 'SAT', 'EDGE', 'NO_DATA',])
            sctrl = afwMath.StatisticsControl()
            sctrl.setAndMask(bad)
            stats = afwMath.makeStatistics(image, afwMath.STDEVCLIP, sctrl)
            thres = stats.getValue(afwMath.STDEVCLIP) * self.config.thresholdValue
            threshold = afwDet.createThreshold(thres, 'value', parity)
            threshold.setIncludeMultiplier(self.config.includeThresholdMultiplier)

        fpSet = afwDet.FootprintSet(image, threshold, maskName, self.config.minPixels)
        return fpSet

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

def getBackground(image, backgroundConfig, nx=0, ny=0, algorithm=None):
    """!Estimate the background of an image (a thin layer on lsst.afw.math.makeBackground)

    \param[in] image  image whose background is to be computed
    \param[in] backgroundConfig  configuration (a BackgroundConfig)
    \param[in] nx  number of x bands; 0 for default
    \param[in] ny  number of y bands; 0 for default
    \param[in] algorithm  name of interpolation algorithm; see lsst.afw.math.BackgroundControl for details
    """
    backgroundConfig.validate();

    logger = pexLogging.getDefaultLog()
    logger = pexLogging.Log(logger,"lsst.meas.algorithms.detection.getBackground")

    if not nx:
        nx = image.getWidth()//backgroundConfig.binSize + 1
    if not ny:
        ny = image.getHeight()//backgroundConfig.binSize + 1

    sctrl = afwMath.StatisticsControl()
    sctrl.setAndMask(reduce(lambda x, y: x | image.getMask().getPlaneBitMask(y),
                            backgroundConfig.ignoredPixelMask, 0x0))
    sctrl.setNanSafe(backgroundConfig.isNanSafe)

    pl = pexLogging.Debug("lsst.meas.algorithms.detection.getBackground")
    pl.debug(3, "Ignoring mask planes: %s" % ", ".join(backgroundConfig.ignoredPixelMask))

    if not algorithm:
        algorithm = backgroundConfig.algorithm

    bctrl = afwMath.BackgroundControl(algorithm, nx, ny,
                                      backgroundConfig.undersampleStyle, sctrl,
                                      backgroundConfig.statisticsProperty)

    # TODO: The following check should really be done within afw/math.  With the
    #       current code structure, it would need to be accounted for in the
    #       doGetImage() funtion in BackgroundMI.cc (which currently only checks
    #       against the interpoation settings which is not appropriate when
    #       useApprox=True) and/or the makeApproximate() function in
    #       afw/Approximate.cc.
    #       See ticket DM-2920: "Clean up code in afw for Approximate background
    #       estimation" (which includes a note to remove the following and the
    #       similar checks in pipe_tasks/matchBackgrounds.py once implemented)
    #
    # Check that config setting of approxOrder/binSize make sense
    # (i.e. ngrid (= shortDimension/binSize) > approxOrderX) and perform
    # appropriate undersampleStlye behavior.
    if backgroundConfig.useApprox:
        if not backgroundConfig.approxOrderY in (backgroundConfig.approxOrderX,-1):
            raise ValueError("Error: approxOrderY not in (approxOrderX, -1)")
        order = backgroundConfig.approxOrderX
        minNumberGridPoints = backgroundConfig.approxOrderX + 1
        if min(nx,ny) <= backgroundConfig.approxOrderX:
            logger.warn("Too few points in grid to constrain fit: min(nx, ny) < approxOrder) "+
                        "[min(%d, %d) < %d]" % (nx, ny, backgroundConfig.approxOrderX))
            if backgroundConfig.undersampleStyle == "THROW_EXCEPTION":
                raise ValueError("Too few points in grid (%d, %d) for order (%d) and binsize (%d)" % (
                        nx, ny, backgroundConfig.approxOrderX, backgroundConfig.binSize))
            elif backgroundConfig.undersampleStyle == "REDUCE_INTERP_ORDER":
                if order < 1:
                    raise ValueError("Cannot reduce approxOrder below 0.  " +
                                     "Try using undersampleStyle = \"INCREASE_NXNYSAMPLE\" instead?")
                order = min(nx, ny) - 1
                logger.warn("Reducing approxOrder to %d" % order)
            elif backgroundConfig.undersampleStyle == "INCREASE_NXNYSAMPLE":
                newBinSize = min(image.getWidth(),image.getHeight())//(minNumberGridPoints-1)
                if newBinSize < 1:
                    raise ValueError("Binsize must be greater than 0")
                newNx = image.getWidth()//newBinSize + 1
                newNy = image.getHeight()//newBinSize + 1
                bctrl.setNxSample(newNx)
                bctrl.setNySample(newNy)
                logger.warn("Decreasing binSize from %d to %d for a grid of (%d, %d)" %
                            (backgroundConfig.binSize, newBinSize, newNx, newNy))

        actrl = afwMath.ApproximateControl(afwMath.ApproximateControl.CHEBYSHEV, order, order,
                                           backgroundConfig.weighting)
        bctrl.setApproximateControl(actrl)

    return afwMath.makeBackground(image, bctrl)

getBackground.ConfigClass = BackgroundConfig

def estimateBackground(exposure, backgroundConfig, subtract=True, stats=True,
                       statsKeys=None):
    """!Estimate exposure's background using parameters in backgroundConfig.

    If subtract is true, make a copy of the exposure and subtract the background.
    If `stats` is True, measure the mean and variance of the background and
    add them to the background-subtracted exposure's metadata with keys
    "BGMEAN" and "BGVAR", or the keys given in `statsKeys` (2-tuple of strings).

    Return background, backgroundSubtractedExposure
    """

    displayBackground = lsstDebug.Info(__name__).displayBackground

    maskedImage = exposure.getMaskedImage()

    background = getBackground(maskedImage, backgroundConfig)

    if not background:
        raise RuntimeError, "Unable to estimate background for exposure"

    bgimg = None

    if displayBackground > 1:
        bgimg = background.getImageF()
        ds9.mtv(bgimg, title="background", frame=1)

    if not subtract:
        return background, None

    bbox = maskedImage.getBBox()
    backgroundSubtractedExposure = exposure.Factory(exposure, bbox, afwImage.PARENT, True)
    copyImage = backgroundSubtractedExposure.getMaskedImage().getImage()
    if bgimg is None:
        bgimg = background.getImageF()
    copyImage -= bgimg

    # Record statistics of the background in the bgsub exposure metadata.
    # (lsst.daf.base.PropertySet)
    if stats:
        if statsKeys is None:
            mnkey  = 'BGMEAN'
            varkey = 'BGVAR'
        else:
            mnkey, varkey = statsKeys
        meta = backgroundSubtractedExposure.getMetadata()
        s = afwMath.makeStatistics(bgimg, afwMath.MEAN | afwMath.VARIANCE)
        bgmean = s.getValue(afwMath.MEAN)
        bgvar = s.getValue(afwMath.VARIANCE)
        meta.addDouble(mnkey, bgmean)
        meta.addDouble(varkey, bgvar)

    if displayBackground:
        ds9.mtv(backgroundSubtractedExposure, title="subtracted")

    return background, backgroundSubtractedExposure
estimateBackground.ConfigClass = BackgroundConfig
