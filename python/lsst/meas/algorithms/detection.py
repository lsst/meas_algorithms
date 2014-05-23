 
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2011 LSST Corporation.
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
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import numpy

import lsstDebug
import lsst.pex.logging as pexLogging 

import lsst.pex.config as pexConfig
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.detection as afwDet
import lsst.pipe.base as pipeBase

from . import algorithmsLib

__all__ = ("SourceDetectionConfig", "SourceDetectionTask", "getBackground",
           "estimateBackground", "BackgroundConfig", "addExposures")

import lsst.daf.persistence as dafPersist
import lsst.pex.config as pexConfig
import lsst.afw.detection as afwDet
import lsst.afw.display.ds9 as ds9
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath

class BackgroundConfig(pexConfig.Config):
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

    def validate(self):
        pexConfig.Config.validate(self)
        # Allow None to be used as an equivalent for "NONE", even though C++ expects the latter.
        if self.algorithm is None:
            self.algorithm = "NONE"

class SourceDetectionConfig(pexConfig.Config):
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
        dtype=bool, optional=False, default=True    # TODO: set default to False once we have a deblender; ticket #2138
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
        default = 0.0,
    )
    reEstimateBackground = pexConfig.Field(
        dtype = bool,
        doc = "Estimate the background again after final source detection?",
        default = True, optional=False,
    )
    background = pexConfig.ConfigField(
        dtype=BackgroundConfig,
        doc="Background re-estimation configuration"
        )

class SourceDetectionTask(pipeBase.Task):
    """
    Detect positive and negative sources on an exposure and return a new SourceCatalog.
    """
    ConfigClass = SourceDetectionConfig
    _DefaultName = "sourceDetection"

    def __init__(self, schema=None, **kwds):
        """Create the detection task.  Most arguments are simply passed onto pipe_base.Task.

        If schema is not None, it will be used to register a 'flags.negative' flag field
        that will be set for negative detections.
        """
        pipeBase.Task.__init__(self, **kwds)
        if schema is not None:
            self.negativeFlagKey = schema.addField(
                "flags.negative", type="Flag",
                doc="set if source was detected as significantly negative"
                )
        else:
            if self.config.thresholdPolarity == "both":
                self.log.log(self.log.WARN, "Detection polarity set to 'both', but no flag will be "\
                             "set to distinguish between positive and negative detections")
            self.negativeFlagKey = None

    @pipeBase.timeMethod
    def makeSourceCatalog(self, table, exposure, doSmooth=True, sigma=None, clearMask=True):
        """Run source detection and create a SourceCatalog.

        To avoid dealing with sources and tables, use detectFootprints() to just get the FootprintSets.

        @param table    lsst.afw.table.SourceTable object that will be used to created the SourceCatalog.
        @param exposure Exposure to process; DETECTED mask plane will be set in-place.
        @param doSmooth if True, smooth the image before detection using a Gaussian of width sigma
        @param sigma    sigma of PSF (pixels); used for smoothing and to grow detections;
            if None then measure the sigma of the PSF of the exposure
        @param clearMask Clear DETECTED{,_NEGATIVE} planes before running detection
        
        @return a Struct with:
          sources -- an lsst.afw.table.SourceCatalog object
          fpSets --- Struct returned by detectFootprints
        
        @raise pipe_base TaskError if sigma=None, doSmooth=True and the exposure has no PSF
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

    @pipeBase.timeMethod
    def detectFootprints(self, exposure, doSmooth=True, sigma=None, clearMask=True):
        """Detect footprints.

        @param exposure Exposure to process; DETECTED{,_NEGATIVE} mask plane will be set in-place.
        @param doSmooth if True, smooth the image before detection using a Gaussian of width sigma
        @param sigma    sigma of PSF (pixels); used for smoothing and to grow detections;
            if None then measure the sigma of the PSF of the exposure
        @param clearMask Clear both DETECTED and DETECTED_NEGATIVE planes before running detection

        @return a lsst.pipe.base.Struct with fields:
        - positive: lsst.afw.detection.FootprintSet with positive polarity footprints (may be None)
        - negative: lsst.afw.detection.FootprintSet with negative polarity footprints (may be None)
        - numPos: number of footprints in positive or 0 if detection polarity was negative
        - numNeg: number of footprints in negative or 0 if detection polarity was positive
        - background: re-estimated background.  None if reEstimateBackground==False
        
        @raise pipe_base TaskError if sigma=None and the exposure has no PSF
        """
        try:
            import lsstDebug
            display = lsstDebug.Info(__name__).display
        except ImportError, e:
            try:
                display
            except NameError:
                display = False

        if exposure is None:
            raise RuntimeException("No exposure for detection")

        maskedImage = exposure.getMaskedImage()
        region = maskedImage.getBBox(afwImage.PARENT)

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
            if not numpy.isfinite(sigma):
                raise RuntimeError("Non-finite PSF width: %f" % sigma)

        self.metadata.set("sigma", sigma)
        self.metadata.set("doSmooth", doSmooth)
        
        if not doSmooth:
            convolvedImage = maskedImage.Factory(maskedImage)
            middle = convolvedImage
        else:
            # smooth using a Gaussian (which is separate, hence fast) of width sigma
            # make a SingleGaussian (separable) kernel with the 'sigma'
            psf = exposure.getPsf()
            kWidth = (int(sigma * 7 + 0.5) / 2) * 2 + 1 # make sure it is odd
            self.metadata.set("smoothingKernelWidth", kWidth)
            gaussFunc = afwMath.GaussianFunction1D(sigma)
            gaussKernel = afwMath.SeparableKernel(kWidth, kWidth, gaussFunc, gaussFunc)

            convolvedImage = maskedImage.Factory(maskedImage.getBBox(afwImage.PARENT))

            afwMath.convolve(convolvedImage, maskedImage, gaussKernel, afwMath.ConvolutionControl())
            #
            # Only search psf-smooth part of frame
            #
            goodBBox = gaussKernel.shrinkBBox(convolvedImage.getBBox(afwImage.PARENT))
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
                fpSet = afwDet.FootprintSet(fpSet, nGrow, False)
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

            if middle and display and display > 1:
                ds9.mtv(middle, frame=2, title="middle")

        return fpSets

    def thresholdImage(self, image, thresholdParity, maskName="DETECTED"):
        """Threshold the convolved image, returning a FootprintSet.
        Helper function for detect().

        @param image The (optionally convolved) MaskedImage to threshold
        @param thresholdParity Parity of threshold
        @param maskName Name of mask to set

        @return FootprintSet
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
            threshold = afwDet.createThreshold(thres, 'value')
            threshold.setIncludeMultiplier(self.config.includeThresholdMultiplier)

        fpSet = afwDet.FootprintSet(image, threshold, maskName, self.config.minPixels)
        return fpSet

    @staticmethod
    def setEdgeBits(maskedImage, goodBBox, edgeBitmask):
        """Set the edgeBitmask bits for all of maskedImage outside goodBBox"""
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
    """
    Add a set of exposures together. 
    Assumes that all exposures in set have the same dimensions
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
    """
    Make a new Exposure which is exposure - background
    """
    backgroundConfig.validate();

    if not nx:
        nx = image.getWidth()//backgroundConfig.binSize + 1
    if not ny:
        ny = image.getHeight()//backgroundConfig.binSize + 1

    displayBackground = lsstDebug.Info(__name__).displayBackground
    if displayBackground:
        import itertools
        ds9.mtv(image, frame=1)
        xPosts = numpy.rint(numpy.linspace(0, image.getWidth() + 1, num=nx, endpoint=True))
        yPosts = numpy.rint(numpy.linspace(0, image.getHeight() + 1, num=ny, endpoint=True))
        with ds9.Buffering():
            for (xMin, xMax), (yMin, yMax) in itertools.product(zip(xPosts[:-1], xPosts[1:]),
                                                                zip(yPosts[:-1], yPosts[1:])):
                ds9.line([(xMin, yMin), (xMin, yMax), (xMax, yMax), (xMax, yMin), (xMin, yMin)], frame=1)


    sctrl = afwMath.StatisticsControl()
    sctrl.setAndMask(reduce(lambda x, y: x | image.getMask().getPlaneBitMask(y),
                            backgroundConfig.ignoredPixelMask, 0x0))
    sctrl.setNanSafe(backgroundConfig.isNanSafe)

    pl = pexLogging.Debug("meas.utils.sourceDetection.getBackground")
    pl.debug(3, "Ignoring mask planes: %s" % ", ".join(backgroundConfig.ignoredPixelMask))

    if not algorithm:
        algorithm = backgroundConfig.algorithm
        
    bctrl = afwMath.BackgroundControl(algorithm, nx, ny,
                                      backgroundConfig.undersampleStyle, sctrl,
                                      backgroundConfig.statisticsProperty)

    return afwMath.makeBackground(image, bctrl)

getBackground.ConfigClass = BackgroundConfig
    
def estimateBackground(exposure, backgroundConfig, subtract=True, stats=True,
                       statsKeys=None):
    """
    Estimate exposure's background using parameters in backgroundConfig.  
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
        ds9.mtv(bgimg, title="background", frame=3)

    if not subtract:
        return background, None

    bbox = maskedImage.getBBox(afwImage.PARENT)
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
            mnkey,varkey = statsKeys
        meta = backgroundSubtractedExposure.getMetadata()
        s = afwMath.makeStatistics(bgimg, afwMath.MEAN | afwMath.VARIANCE)
        bgmean = s.getValue(afwMath.MEAN)
        bgvar  = s.getValue(afwMath.VARIANCE)
        meta.addDouble(mnkey,  bgmean)
        meta.addDouble(varkey,  bgvar)
    
    if displayBackground:
        ds9.mtv(backgroundSubtractedExposure, title="subtracted")

    return background, backgroundSubtractedExposure
estimateBackground.ConfigClass = BackgroundConfig
