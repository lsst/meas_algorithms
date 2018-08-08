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
__all__ = ("SubtractBackgroundConfig", "SubtractBackgroundTask")

import itertools

import numpy

from lsstDebug import getDebugFrame
import lsst.afw.display as afwDisplay
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from functools import reduce


class SubtractBackgroundConfig(pexConfig.Config):
    """Config for SubtractBackgroundTask

    Notes
    -----
    Many of these fields match fields in lsst.afw.math.BackgroundControl,
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
        dtype=int, default=128, min=1,
    )
    binSizeX = pexConfig.RangeField(
        doc=("Sky region size to be used for each background point in X direction. "
             "If 0, the binSize config is used."),
        dtype=int, default=0, min=0,
    )
    binSizeY = pexConfig.RangeField(
        doc=("Sky region size to be used for each background point in Y direction. "
             "If 0, the binSize config is used."),
        dtype=int, default=0, min=0,
    )
    algorithm = pexConfig.ChoiceField(
        doc="how to interpolate the background values. This maps to an enum; see afw::math::Background",
        dtype=str, default="AKIMA_SPLINE", optional=True,
        allowed={
            "CONSTANT": "Use a single constant value",
            "LINEAR": "Use linear interpolation",
            "NATURAL_SPLINE": "cubic spline with zero second derivative at endpoints",
            "AKIMA_SPLINE": "higher-level nonlinear spline that is more robust to outliers",
            "NONE": "No background estimation is to be attempted",
        },
    )
    ignoredPixelMask = pexConfig.ListField(
        doc="Names of mask planes to ignore while estimating the background",
        dtype=str, default=["BAD", "EDGE", "DETECTED", "DETECTED_NEGATIVE", "NO_DATA", ],
        itemCheck=lambda x: x in afwImage.Mask().getMaskPlaneDict().keys(),
    )
    isNanSafe = pexConfig.Field(
        doc="Ignore NaNs when estimating the background",
        dtype=bool, default=False,
    )

    useApprox = pexConfig.Field(
        doc="Use Approximate (Chebyshev) to model background.",
        dtype=bool, default=True,
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


## @addtogroup LSST_task_documentation
## @{
## @page SubtractBackgroundTask
## @ref SubtractBackgroundTask_ "SubtractBackgroundTask"
## @copybrief SubtractBackgroundTask
## @}

class SubtractBackgroundTask(pipeBase.Task):
    """Subtract the background from an exposure

    Notes
    -----
    The `run` method will optionally set the following items of exposure metadata;
    the names may be overridden
    
    BGMEAN <dd>mean value of background
    BGVAR  <dd>standard deviation of background

    SubtractBackgroundTask has a debug dictionary containing three integer keys:
    
    unsubtracted
    If >0: `fitBackground` displays the unsubtracted masked image overlaid with the grid of cells
                used to fit the background in the specified frame
    subtracted
    If >0: `run` displays the background-subtracted exposure is the specified frame
    background
    If >0: `run` displays the background image in the specified frame
    
    Examples
    --------
    For example, put something like
        import lsstDebug
        def DebugInfo(name):
            di = lsstDebug.getInfo(name)  # N.b. lsstDebug.Info(name) would call us recursively
            if name == "lsst.meas.algorithms.subtractBackground":
                di.display = dict(
                    unsubtracted = 1,
                    subtracted = 2,
                    background = 3,
                )

            return di

        lsstDebug.Info = DebugInfo
    into your `debug.py` file and run your task with the `--debug` flag.
    """
    ConfigClass = SubtractBackgroundConfig
    _DefaultName = "subtractBackground"

    def run(self, exposure, background=None, stats=True, statsKeys=None):
        """Fit and subtract the background of an exposure

        Parameters
        -----------
        exposure:  
        exposure whose background is to be subtracted

        background:  
        initial background model already subtracted from exposure
            (an lsst.afw.math.BackgroundList). May be None if no background has been subtracted.
        stats:
        if True then measure the mean and variance of the full background model
                        and record the results in the exposure's metadata
        statsKeys:
        key names used to store the mean and variance of the background
            in the exposure's metadata (a pair of strings); if None then use ("BGMEAN", "BGVAR");
            ignored if stats is false

        Returns
        ----------
        pipeBase.Struct():
        an lsst.pipe.base.Struct containing:
        - background  full background model (initial model with changes), an lsst.afw.math.BackgroundList
        """
        if background is None:
            background = afwMath.BackgroundList()

        maskedImage = exposure.getMaskedImage()
        fitBg = self.fitBackground(maskedImage)
        maskedImage -= fitBg.getImageF()
        background.append(fitBg)

        if stats:
            self._addStats(exposure, background, statsKeys=statsKeys)

        subFrame = getDebugFrame(self._display, "subtracted")
        if subFrame:
            subDisp = afwDisplay.getDisplay(frame=subFrame)
            subDisp.mtv(exposure, title="subtracted")

        bgFrame = getDebugFrame(self._display, "background")
        if bgFrame:
            bgDisp = afwDisplay.getDisplay(frame=bgFrame)
            bgImage = background.getImage()
            bgDisp.mtv(bgImage, title="background")

        return pipeBase.Struct(
            background=background,
        )

    def _addStats(self, exposure, background, statsKeys=None):
        """Add statistics about the background to the exposure's metadata

        Parameters
        ----------
        exposure:  
        exposure whose background was subtracted

        background:  
        background model (an lsst.afw.math.BackgroundList)

        statsKeys: 
        key names used to store the mean and variance of the background
            in the exposure's metadata (a pair of strings); if None then use ("BGMEAN", "BGVAR");
            ignored if stats is false
        """
        netBgImg = background.getImage()
        if statsKeys is None:
            statsKeys = ("BGMEAN", "BGVAR")
        mnkey, varkey = statsKeys
        meta = exposure.getMetadata()
        s = afwMath.makeStatistics(netBgImg, afwMath.MEAN | afwMath.VARIANCE)
        bgmean = s.getValue(afwMath.MEAN)
        bgvar = s.getValue(afwMath.VARIANCE)
        meta.addDouble(mnkey, bgmean)
        meta.addDouble(varkey, bgvar)

    def fitBackground(self, maskedImage, nx=0, ny=0, algorithm=None):
        """Estimate the background of a masked image

        Parameters
        ----------
        maskedImage:  
            masked image whose background is to be computed

        nx:  
            number of x bands; if 0 compute from width and config.binSizeX

        ny:  
            number of y bands; if 0 compute from height and config.binSizeY

        algorithm:  
            name of interpolation algorithm; if None use self.config.algorithm

        Returns
        -----------
        bg:
            fit background as an lsst.afw.math.Background

        Raises
        -----------
        RuntimeError 
            if lsst.afw.math.makeBackground returns None,
            which is apparently one way it indicates failure"""
        binSizeX = self.config.binSize if self.config.binSizeX == 0 else self.config.binSizeX
        binSizeY = self.config.binSize if self.config.binSizeY == 0 else self.config.binSizeY

        if not nx:
            nx = maskedImage.getWidth()//binSizeX + 1
        if not ny:
            ny = maskedImage.getHeight()//binSizeY + 1

        unsubFrame = getDebugFrame(self._display, "unsubtracted")
        if unsubFrame:
            unsubDisp = afwDisplay.getDisplay(frame=unsubFrame)
            unsubDisp.mtv(maskedImage, title="unsubtracted")
            xPosts = numpy.rint(numpy.linspace(0, maskedImage.getWidth() + 1, num=nx, endpoint=True))
            yPosts = numpy.rint(numpy.linspace(0, maskedImage.getHeight() + 1, num=ny, endpoint=True))
            with unsubDisp.Buffering():
                for (xMin, xMax), (yMin, yMax) in itertools.product(zip(xPosts[:-1], xPosts[1:]),
                                                                    zip(yPosts[:-1], yPosts[1:])):
                    unsubDisp.line([(xMin, yMin), (xMin, yMax), (xMax, yMax), (xMax, yMin), (xMin, yMin)])

        sctrl = afwMath.StatisticsControl()
        sctrl.setAndMask(reduce(lambda x, y: x | maskedImage.getMask().getPlaneBitMask(y),
                                self.config.ignoredPixelMask, 0x0))
        sctrl.setNanSafe(self.config.isNanSafe)

        self.log.debug("Ignoring mask planes: %s" % ", ".join(self.config.ignoredPixelMask))

        if algorithm is None:
            algorithm = self.config.algorithm

        bctrl = afwMath.BackgroundControl(algorithm, nx, ny,
                                          self.config.undersampleStyle, sctrl,
                                          self.config.statisticsProperty)

        # TODO: The following check should really be done within lsst.afw.math.
        #       With the current code structure, it would need to be accounted for in the doGetImage()
        #       function in BackgroundMI.cc (which currently only checks against the interpolation settings,
        #       which is not appropriate when useApprox=True)
        #       and/or the makeApproximate() function in afw/Approximate.cc.
        #       See ticket DM-2920: "Clean up code in afw for Approximate background
        #       estimation" (which includes a note to remove the following and the
        #       similar checks in pipe_tasks/matchBackgrounds.py once implemented)
        #
        # Check that config setting of approxOrder/binSize make sense
        # (i.e. ngrid (= shortDimension/binSize) > approxOrderX) and perform
        # appropriate undersampleStlye behavior.
        if self.config.useApprox:
            if self.config.approxOrderY not in (self.config.approxOrderX, -1):
                raise ValueError("Error: approxOrderY not in (approxOrderX, -1)")
            order = self.config.approxOrderX
            minNumberGridPoints = order + 1
            if min(nx, ny) <= order:
                self.log.warn("Too few points in grid to constrain fit: min(nx, ny) < approxOrder) "
                              "[min(%d, %d) < %d]" % (nx, ny, order))
                if self.config.undersampleStyle == "THROW_EXCEPTION":
                    raise ValueError("Too few points in grid (%d, %d) for order (%d) and binSize (%d, %d)" %
                                     (nx, ny, order, binSizeX, binSizeY))
                elif self.config.undersampleStyle == "REDUCE_INTERP_ORDER":
                    if order < 1:
                        raise ValueError("Cannot reduce approxOrder below 0.  " +
                                         "Try using undersampleStyle = \"INCREASE_NXNYSAMPLE\" instead?")
                    order = min(nx, ny) - 1
                    self.log.warn("Reducing approxOrder to %d" % order)
                elif self.config.undersampleStyle == "INCREASE_NXNYSAMPLE":
                    # Reduce bin size to the largest acceptable square bins
                    newBinSize = min(maskedImage.getWidth(), maskedImage.getHeight())//(minNumberGridPoints-1)
                    if newBinSize < 1:
                        raise ValueError("Binsize must be greater than 0")
                    newNx = maskedImage.getWidth()//newBinSize + 1
                    newNy = maskedImage.getHeight()//newBinSize + 1
                    bctrl.setNxSample(newNx)
                    bctrl.setNySample(newNy)
                    self.log.warn("Decreasing binSize from (%d, %d) to %d for a grid of (%d, %d)" %
                                  (binSizeX, binSizeY, newBinSize, newNx, newNy))

            actrl = afwMath.ApproximateControl(afwMath.ApproximateControl.CHEBYSHEV, order, order,
                                               self.config.weighting)
            bctrl.setApproximateControl(actrl)

        bg = afwMath.makeBackground(maskedImage, bctrl)
        if bg is None:
            raise RuntimeError("lsst.afw.math.makeBackground failed to fit a background model")
        return bg
