#
# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
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
import itertools

import numpy

import lsstDebug
import lsst.afw.display.ds9 as ds9
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

__all__ = ("SubtractBackgroundConfig", "SubtractBackgroundTask")

class SubtractBackgroundConfig(pexConfig.Config):
    """!Config for SubtractBackgroundTask

    @note Many of these fields match fields in lsst.afw.math.BackgroundControl,
    the control class for lsst.afw.math.makeBackground
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
            },
    )
    binSize = pexConfig.RangeField(
        doc="how large a region of the sky should be used for each background point",
        dtype=int, default=256, min=10,
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
            },
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

class SubtractBackgroundTask(pipeBase.Task):
    ConfigClass = SubtractBackgroundConfig
    _DefaultName = "subtractBackground"

    def run(self, exposure, background=None, stats=True, statsKeys=None):
        """!Fit and subtract the background of an exposure

        @param[in,out] exposure  exposure whose background is to be subtracted
        @param[in,out] background  initial background model already subtracted from exposure
            (an lsst.afw.math.BackgroundList). May be None if no background has been subtracted.
        @param[in] stats  if True then measure the mean and variance of the full background model
                        and record the results in the exposure's metadata
        @param[in] statsKeys  key names used to store the mean and variance of the background
            in the exposure's metadata (a pair of strings); if None then use ("BGMEAN", "BGVAR");
            ignored if stats is false

        @return an lsst.pipe.base.Struct containing:
        - background  full background model (initial model with changes), an lsst.afw.math.BackgroundList

        @throw RuntimeError if fitBackground returns a background of None
        """
        displayBackground = lsstDebug.Info(__name__).displayBackground

        if background is None:
            background = afwMath.BackgroundList()

        maskedImage = exposure.getMaskedImage()

        fitBg = self.fitBackground(maskedImage)
        if not fitBg:
            raise RuntimeError("Unable to estimate background for exposure")

        fitBgImg = fitBg.getImageF()
        maskedImage -= fitBgImg
        background.append(fitBg)
        del fitBgImg

        if displayBackground > 1 or stats:
            netBgImg = background.getImage()

        if displayBackground > 1:
            ds9.mtv(netBgImg, title="background", frame=3)

        # Record statistics of the background in the bgsub exposure metadata.
        # (lsst.daf.base.PropertySet)
        if stats:
            if statsKeys is None:
                statsKeys = ("BGMEAN", "BGVAR")
            mnkey, varkey = statsKeys
            meta = exposure.getMetadata()
            s = afwMath.makeStatistics(netBgImg, afwMath.MEAN | afwMath.VARIANCE)
            bgmean = s.getValue(afwMath.MEAN)
            bgvar = s.getValue(afwMath.VARIANCE)
            meta.addDouble(mnkey, bgmean)
            meta.addDouble(varkey, bgvar)

        if displayBackground:
            ds9.mtv(exposure, title="subtracted")

        return pipeBase.Struct(
            background = background,
        )

    def fitBackground(self, maskedImage, nx=0, ny=0, algorithm=None):
        """!Estimate the background of a masked image

        This is a thin layer on lsst.afw.math.makeBackground

        @param[in] maskedImage  masked image whose background is to be computed
        @param[in] nx  number of x bands; if 0 compute from width and config.binSize
        @param[in] ny  number of y bands; if 0 compute from height and config.binSize
        @param[in] algorithm  name of interpolation algorithm; if None use self.config.algorithm

        @return fit background as an lsst.afw.math.Background, or None if the estimation fails
                    (unless the code raises an exception)
        """
        if not nx:
            nx = maskedImage.getWidth()//self.config.binSize + 1
        if not ny:
            ny = maskedImage.getHeight()//self.config.binSize + 1

        displayBackground = lsstDebug.Info(__name__).displayBackground
        if displayBackground:
            ds9.mtv(maskedImage, frame=1)
            xPosts = numpy.rint(numpy.linspace(0, maskedImage.getWidth() + 1, num=nx, endpoint=True))
            yPosts = numpy.rint(numpy.linspace(0, maskedImage.getHeight() + 1, num=ny, endpoint=True))
            with ds9.Buffering():
                for (xMin, xMax), (yMin, yMax) in itertools.product(zip(xPosts[:-1], xPosts[1:]),
                                                                    zip(yPosts[:-1], yPosts[1:])):
                    ds9.line([(xMin, yMin), (xMin, yMax), (xMax, yMax), (xMax, yMin), (xMin, yMin)], frame=1)


        sctrl = afwMath.StatisticsControl()
        sctrl.setAndMask(reduce(lambda x, y: x | maskedImage.getMask().getPlaneBitMask(y),
                                self.config.ignoredPixelMask, 0x0))
        sctrl.setNanSafe(self.config.isNanSafe)

        self.log.logdebug("Ignoring mask planes: %s" % ", ".join(self.config.ignoredPixelMask))

        if algorithm is None:
            algorithm = self.config.algorithm

        bctrl = afwMath.BackgroundControl(algorithm, nx, ny,
                                          self.config.undersampleStyle, sctrl,
                                          self.config.statisticsProperty)

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
        if self.config.useApprox:
            if self.config.approxOrderY not in (self.config.approxOrderX,-1):
                raise ValueError("Error: approxOrderY not in (approxOrderX, -1)")
            order = self.config.approxOrderX
            minNumberGridPoints = self.config.approxOrderX + 1
            if min(nx,ny) <= self.config.approxOrderX:
                self.log.warn("Too few points in grid to constrain fit: min(nx, ny) < approxOrder) "+
                            "[min(%d, %d) < %d]" % (nx, ny, self.config.approxOrderX))
                if self.config.undersampleStyle == "THROW_EXCEPTION":
                    raise ValueError("Too few points in grid (%d, %d) for order (%d) and binsize (%d)" % (
                            nx, ny, self.config.approxOrderX, self.config.binSize))
                elif self.config.undersampleStyle == "REDUCE_INTERP_ORDER":
                    if order < 1:
                        raise ValueError("Cannot reduce approxOrder below 0.  " +
                                         "Try using undersampleStyle = \"INCREASE_NXNYSAMPLE\" instead?")
                    order = min(nx, ny) - 1
                    self.log.warn("Reducing approxOrder to %d" % order)
                elif self.config.undersampleStyle == "INCREASE_NXNYSAMPLE":
                    newBinSize = min(maskedImage.getWidth(),maskedImage.getHeight())//(minNumberGridPoints-1)
                    if newBinSize < 1:
                        raise ValueError("Binsize must be greater than 0")
                    newNx = maskedImage.getWidth()//newBinSize + 1
                    newNy = maskedImage.getHeight()//newBinSize + 1
                    bctrl.setNxSample(newNx)
                    bctrl.setNySample(newNy)
                    self.log.warn("Decreasing binSize from %d to %d for a grid of (%d, %d)" %
                                (self.config.binSize, newBinSize, newNx, newNy))

            actrl = afwMath.ApproximateControl(afwMath.ApproximateControl.CHEBYSHEV, order, order,
                                               self.config.weighting)
            bctrl.setApproximateControl(actrl)

        return afwMath.makeBackground(maskedImage, bctrl)
