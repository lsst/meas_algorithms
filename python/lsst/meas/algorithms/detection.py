# 
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

import lsst.pex.config as pexConfig
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.detection as afwDet
import lsst.pipe.base as pipeBase

from . import algorithmsLib

__all__ = "SourceDetectionConfig", "SourceDetectionTask"

class SourceDetectionConfig(pexConfig.Config):
    minPixels = pexConfig.RangeField(
        doc="detected sources with fewer than the specified number of pixels will be ignored",
        dtype=int, optional=True, default=1, min=0,
    )
    nGrow = pexConfig.RangeField(
        doc="How many pixels to to grow detections",
        dtype=int, optional=True, default=1, min=0,
    )
    thresholdValue = pexConfig.RangeField(
        doc="Threshold for footprints",
        dtype=float, optional=True, default=5.0, min=0.0,
    )
    includeThresholdMultiplier = pexConfig.RangeField(
        doc="Include threshold relative to thresholdValue",
        dtype=float, default=1.0, min=0.0,
        )        
    thresholdType = pexConfig.ChoiceField(
        doc="specifies the desired flavor of Threshold",
        dtype=str, optional=True, default="stdev",
        allowed={
            "variance": "threshold applied to image variance",
            "stdev": "threshold applied to image std deviation",
            "value": "threshold applied to image value"
        }
    )
    thresholdPolarity = pexConfig.ChoiceField(
        doc="specifies whether to detect positive, or negative sources, or both",
        dtype=str, optional=True, default="positive",
        allowed={
            "positive": "detect only positive sources",
            "negative": "detect only negative sources",
            "both": "detect both positive and negative sources",
        }
    )

class SourceDetectionTask(pipeBase.Task):
    """
    Detect positive and negative sources on an exposure and return a new SourceCatalog.
    """
    ConfigClass = SourceDetectionConfig

    def __init__(self, config=None, schema=None, **kwds):
        """Create the detection task.  Most arguments are simply passed onto pipe_base.Task.

        If schema is not None, it will be used to register a 'flags.negative' flag field
        that will be set for negative detections.
        """
        pipeBase.Task.__init__(self, config, **kwds)
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
    def makeSourceCatalog(self, table, exposure):
        """Run source detection and create a SourceCatalog.

        To avoid dealing with sources and tables, use detectFootprints() to just get the FootprintSets.

        @param table    lsst.afw.table.SourceTable object that will be used to created the SourceCatalog.
        @param exposure Exposure to process; DETECTED mask plane will be set in-place.
        
        @return an lsst.afw.table.SourceCatalog object
        """
        assert exposure, "No exposure provided"
        assert self.negativeFlagKey is None or self.negativeFlagKey in table.getSchema(), \
            "Table has incorrect Schema"
        fpSets = self.detectFootprints(exposure)
        sources = afwTable.SourceCatalog(table)
        table.preallocate(fpSets.numPos + fpSets.numNeg) # not required, but nice
        if fpSets.negative:
            fpSets.positive.makeSources(sources)
            if self.negativeFlagKey:
                for record in sources:
                    record.set(self.negativeFlagKey, True)
        if fpSets.positive:
            fpSets.positive.makeSources(sources)
        return sources

    @pipeBase.timeMethod
    def detectFootprints(self, exposure):
        """Detect footprints.

        @param exposure Exposure to process; DETECTED mask plane will be set in-place.

        @return a lsst.pipe.base.Struct with fields:
        - positive: lsst.afw.detection.FootprintSet with positive polarity footprints (may be None)
        - negative: lsst.afw.detection.FootprintSet with negative polarity footprints (may be None)
        - numPos: number of footprints in positive or 0 if detection polarity was negative
        - numNeg: number of footprints in negative or 0 if detection polarity was positive
        """
        assert exposure, "No exposure provided"

        maskedImage = exposure.getMaskedImage()
        region = maskedImage.getBBox(afwImage.PARENT)

        mask = maskedImage.getMask()
        # not resetting DETECTED_NEGATIVE anymore because we never set it later (?)
        mask &= ~(mask.getPlaneBitMask("DETECTED"))

        psf = exposure.getPsf()

        if psf is None:
            convolvedImage = maskedImage.Factory(maskedImage)
            middle = convolvedImage
        else:
            # use a separable psf for convolution ... the psf width for the center of the image will do

            xCen = maskedImage.getX0() + maskedImage.getWidth()/2
            yCen = maskedImage.getY0() + maskedImage.getHeight()/2

            # measure the 'sigma' of the psf we were given
            psfAttrib = algorithmsLib.PsfAttributes(psf, xCen, yCen)
            sigma = psfAttrib.computeGaussianWidth()

            # make a SingleGaussian (separable) kernel with the 'sigma'
            gaussFunc = afwMath.GaussianFunction1D(sigma)
            gaussKernel = afwMath.SeparableKernel(psf.getKernel().getWidth(), psf.getKernel().getHeight(),
                                                  gaussFunc, gaussFunc)

            convolvedImage = maskedImage.Factory(maskedImage.getDimensions())
            convolvedImage.setXY0(maskedImage.getXY0())

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

        # Using dicts here so we can reassign them to new FootprintSets later;
        # old implementation wasn't returning the footprints that were actually
        # "growed" because it assigned the result to a (temporary) loop variable.

        fpSets = pipeBase.Struct(positive=None, negative=None)
        if self.config.thresholdPolarity != "negative":
            fpSets.positive = self.thresholdImage(middle, "positive")
        if self.config.thresholdPolarity != "positive":
            fpSets.negative = self.thresholdImage(middle, "negative")

        for polarity in ("positive", "negative"):
            fpSet = getattr(fpSets, polarity)
            if fpSet is None:
                continue
            fpSet.setRegion(region)
            if self.config.nGrow > 0:
                setattr(fpSets, polarity, afwDet.FootprintSet(fpSet, self.config.nGrow, False))
            fpSet.setMask(maskedImage.getMask(), "DETECTED")

        fpSets.numPos = len(fpSets.positive.getFootprints()) if fpSets.positive is not None else 0
        fpSets.numNeg = len(fpSets.negative.getFootprints()) if fpSets.negative is not None else 0

        self.log.log(self.log.INFO, "Detected %d positive sources to %g sigma." %
                     (fpSets.numPos, self.config.thresholdValue))
        if fpSets.numNeg > 0:
            self.log.log(self.log.INFO, "Detected %d negative sources to %g sigma" %
                         (fpSets.numNeg, self.config.thresholdValue))


        return fpSets

    def thresholdImage(self, image, thresholdParity):
        """Threshold the convolved image, returning a FootprintSet.
        Helper function for detect().

        @param image The (optionally convolved) MaskedImage to threshold
        @param thresholdParity Parity of threshold

        @return FootprintSet
        """
        parity = False if thresholdParity == "negative" else True
        threshold = afwDet.createThreshold(self.config.thresholdValue, self.config.thresholdType, parity)
        threshold.setIncludeMultiplier(self.config.includeThresholdMultiplier)
        fpSet = afwDet.FootprintSet(image, threshold, "DETECTED", self.config.minPixels)
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
