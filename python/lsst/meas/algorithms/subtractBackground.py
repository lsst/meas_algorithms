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
__all__ = ("SubtractBackgroundConfig", "SubtractBackgroundTask", "backgroundFlatContext")

import itertools

from contextlib import contextmanager
import numpy

from lsstDebug import getDebugFrame
from lsst.utils import suppress_deprecations

import lsst.afw.display as afwDisplay
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase


@contextmanager
def backgroundFlatContext(maskedImage, doApply, backgroundToPhotometricRatio=None):
    """Context manager to convert from photometric-flattened to background-
    flattened image.

    Parameters
    ----------
    maskedImage : `lsst.afw.image.MaskedImage`
        Masked image (image + mask + variance) to convert from a
        photometrically flat image to an image suitable for background
        subtraction.
    doApply : `bool`
        Apply the conversion? If False, this context manager will not
        do anything.
    backgroundToPhotometricRatio : `lsst.afw.image.Image`, optional
        Image to multiply a photometrically-flattened image by to obtain a
        background-flattened image.
        Only used if ``doApply`` is ``True``.

    Yields
    ------
    maskedImage : `lsst.afw.image.MaskedImage`
        Masked image converted into an image suitable for background
        subtraction.

    Raises
    ------
    RuntimeError if doApply is True and no ratio is supplied.
    ValueError if the ratio is not an `lsst.afw.image.Image`.
    """
    if doApply:
        if backgroundToPhotometricRatio is None:
            raise RuntimeError("backgroundFlatContext called with doApply=True, "
                               "but without a backgroundToPhotometricRatio")
        if not isinstance(backgroundToPhotometricRatio, afwImage.Image):
            raise ValueError("The backgroundToPhotometricRatio must be an lsst.afw.image.Image")

        maskedImage *= backgroundToPhotometricRatio

    try:
        yield maskedImage
    finally:
        if doApply:
            maskedImage /= backgroundToPhotometricRatio


class TooManyMaskedPixelsError(pipeBase.AlgorithmError):
    """Raised when all pixels in the image are masked and no background
    can be estimated.
    """
    def metadata(self) -> dict:
        """There is no metadata associated with this error.
        """
        return {}


class SubtractBackgroundConfig(pexConfig.Config):
    """Config for SubtractBackgroundTask

    Many of these fields match fields in `lsst.afw.math.BackgroundControl`,
    the control class for `lsst.afw.math.makeBackground`
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
    doApplyFlatBackgroundRatio = pexConfig.Field(
        doc="Convert from a photometrically flat image to one suitable to background subtraction? "
            "If True, then a backgroundToPhotometricRatio must be supplied to the task run method.",
        dtype=bool,
        default=False,
    )


class SubtractBackgroundTask(pipeBase.Task):
    """Subtract the background from an exposure
    """
    ConfigClass = SubtractBackgroundConfig
    _DefaultName = "subtractBackground"

    def run(self, exposure, background=None, stats=True, statsKeys=None, backgroundToPhotometricRatio=None):
        """Fit and subtract the background of an exposure.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure whose background is to be subtracted.
        background : `lsst.afw.math.BackgroundList`
            Initial background model already subtracted. May be None if no background
            has been subtracted.
        stats : `bool`
            If True then measure the mean and variance of the full background model and
            record the results in the exposure's metadata.
        statsKeys : `tuple`
            Key names used to store the mean and variance of the background in the
            exposure's metadata (another tuple); if None then use ("BGMEAN", "BGVAR");
            ignored if stats is false.
        backgroundToPhotometricRatio : `lsst.afw.image.Image`, optional
            Image to multiply a photometrically-flattened image by to obtain a
            background-flattened image.
            Only used if config.doApplyFlatBackgroundRatio = True.

        Returns
        -------
        background : `lsst.afw.math.BackgroundList`
            Full background model (initial model with changes), contained in an
            `lsst.pipe.base.Struct`.
        """
        if background is None:
            background = afwMath.BackgroundList()

        maskedImage = exposure.maskedImage

        with backgroundFlatContext(
            maskedImage,
            self.config.doApplyFlatBackgroundRatio,
            backgroundToPhotometricRatio=backgroundToPhotometricRatio,
        ):
            fitBg = self.fitBackground(maskedImage)
            maskedImage -= fitBg.getImageF(self.config.algorithm, self.config.undersampleStyle)

            actrl = fitBg.getBackgroundControl().getApproximateControl()
            background.append((fitBg, getattr(afwMath.Interpolate, self.config.algorithm),
                               fitBg.getAsUsedUndersampleStyle(), actrl.getStyle(),
                               actrl.getOrderX(), actrl.getOrderY(), actrl.getWeighting()))

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
        exposure : `lsst.afw.image.Exposure`
            Exposure whose background was subtracted.
        background : `lsst.afw.math.BackgroundList`
            Background model
        statsKeys : `tuple`
            Key names used to store the mean and variance of the background in
            the exposure's metadata (a tuple); if None then use
            ("BGMEAN", "BGVAR"); ignored if stats is false.
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
        maskedImage : `lsst.afw.image.maskedImage`
            Masked image whose background is to be computed
        nx : 'int`
            Number of x bands; if 0 compute from width and `self.config.binSizeX`
        ny : `int`
            Number of y bands; if 0 compute from height and `self.config.binSizeY`
        algorithm : `str`
            Name of interpolation algorithm; if None use `self.config.algorithm`

        Returns
        -------
        bg : `lsst.afw.math.Background`
            A fit background

        Raises
        ------
        RuntimeError
            Raised if lsst.afw.math.makeBackground returns None, an indicator
            of failure.
        """

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
        badMask = maskedImage.mask.getPlaneBitMask(self.config.ignoredPixelMask)

        sctrl.setAndMask(badMask)
        sctrl.setNanSafe(self.config.isNanSafe)

        self.log.debug("Ignoring mask planes: %s", ", ".join(self.config.ignoredPixelMask))
        if (maskedImage.mask.getArray() & badMask).all():
            raise TooManyMaskedPixelsError("All pixels masked. Cannot estimate background.")

        if algorithm is None:
            algorithm = self.config.algorithm

        # TODO: DM-22814. This call to a deprecated BackgroundControl constructor
        # is necessary to support the algorithm parameter; it # should be replaced with
        #
        #     afwMath.BackgroundControl(nx, ny, sctrl, self.config.statisticsProperty)
        #
        # when algorithm has been deprecated and removed.
        with suppress_deprecations():
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
                self.log.warning("Too few points in grid to constrain fit: min(nx, ny) < approxOrder) "
                                 "[min(%d, %d) < %d]", nx, ny, order)
                if self.config.undersampleStyle == "THROW_EXCEPTION":
                    raise ValueError("Too few points in grid (%d, %d) for order (%d) and binSize (%d, %d)" %
                                     (nx, ny, order, binSizeX, binSizeY))
                elif self.config.undersampleStyle == "REDUCE_INTERP_ORDER":
                    if order < 1:
                        raise ValueError("Cannot reduce approxOrder below 0.  "
                                         "Try using undersampleStyle = \"INCREASE_NXNYSAMPLE\" instead?")
                    order = min(nx, ny) - 1
                    self.log.warning("Reducing approxOrder to %d", order)
                elif self.config.undersampleStyle == "INCREASE_NXNYSAMPLE":
                    # Reduce bin size to the largest acceptable square bins
                    newBinSize = min(maskedImage.getWidth(), maskedImage.getHeight())//(minNumberGridPoints-1)
                    if newBinSize < 1:
                        raise ValueError("Binsize must be greater than 0")
                    newNx = maskedImage.getWidth()//newBinSize + 1
                    newNy = maskedImage.getHeight()//newBinSize + 1
                    bctrl.setNxSample(newNx)
                    bctrl.setNySample(newNy)
                    self.log.warning("Decreasing binSize from (%d, %d) to %d for a grid of (%d, %d)",
                                     binSizeX, binSizeY, newBinSize, newNx, newNy)

            actrl = afwMath.ApproximateControl(afwMath.ApproximateControl.CHEBYSHEV, order, order,
                                               self.config.weighting)
            bctrl.setApproximateControl(actrl)

        bg = afwMath.makeBackground(maskedImage, bctrl)
        if bg is None:
            raise RuntimeError("lsst.afw.math.makeBackground failed to fit a background model")
        return bg
