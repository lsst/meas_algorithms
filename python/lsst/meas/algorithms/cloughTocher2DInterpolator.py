# This file is part of meas_algorithms.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

__all__ = (
    "CloughTocher2DInterpolateConfig",
    "CloughTocher2DInterpolateTask",
)


from lsst.pex.config import Config, Field, ListField
from lsst.pipe.base import Task
from scipy.interpolate import CloughTocher2DInterpolator

from . import CloughTocher2DInterpolatorUtils as ctUtils


class CloughTocher2DInterpolateConfig(Config):
    """Config for CloughTocher2DInterpolateTask."""

    badMaskPlanes = ListField[str](
        doc="List of mask planes to interpolate over.",
        default=["BAD", "SAT", "CR"],
    )
    fillValue = Field[float](
        doc="Constant value to fill outside of the convex hull of the good "
        "pixels. A long (longer than twice the ``interpLength``) streak of "
        "bad pixels at an edge will be set to this value.",
        default=0.0,
    )
    interpLength = Field[int](
        doc="Maximum number of pixels away from a bad pixel to include in "
        "building the interpolant. Must be greater than or equal to 1.",
        default=4,
        check=lambda x: x >= 1,
    )
    flipXY = Field[bool](
        doc="Whether to flip the x and y coordinates before constructing the "
        "Delaunay triangulation. This may produce a slightly different result "
        "since the triangulation is not guaranteed to be invariant under "
        "coordinate flips.",
        default=False,
    )


class CloughTocher2DInterpolateTask(Task):
    """Interpolated over bad pixels using CloughTocher2DInterpolator.

    Pixels with mask bits set to any of those listed ``badMaskPlanes`` config
    are considered bad and are interpolated over. All good (non-bad) pixels
    within ``interpLength`` pixels of a bad pixel in either direction are used
    to construct the interpolant.  An extended streak of bad pixels at an edge,
    longer than ``interpLength``, is set to `fillValue`` specified in config.
    """

    ConfigClass = CloughTocher2DInterpolateConfig
    _DefaultName = "cloughTocher2DInterpolate"

    def run(
        self,
        maskedImage,
        badpix=None,
        goodpix=None,
    ):
        """Interpolate over bad pixels in a masked image.

        This modifies the ``image`` attribute of the ``maskedImage`` in place.
        This method returns, and accepts, the coordinates of the bad pixels
        that were interpolated over, and the coordinates and values of the
        good pixels that were used to construct the interpolant. This avoids
        having to search for the bad and the good pixels repeatedly when the
        mask plane is shared among many images, as would be the case with
        noise realizations.

        Parameters
        ----------
        maskedImage : `~lsst.afw.image.MaskedImageF`
            Image on which to perform interpolation (and modify in-place).
        badpix: `numpy.ndarray`, optional
            N x 3 numpy array, where N is the number of bad pixels.
            The coordinates of the bad pixels to interpolate over in the first
            two columns, and the pixel values (unused) in the third column.
            If None, then the coordinates of the bad pixels are determined by
            an exhaustive search over the image. If ``goodpix`` is not
            provided, then this parameter is ignored.
        goodpix: `numpy.ndarray`, optional
            M x 3 numpy array, where M is the number of good pixels.
            The first two columns are the coordinates of the good pixels around
            ``badpix`` that must be included when constructing the interpolant.
            the interpolant. The values are populated from the image plane of
            the ``maskedImage`` and returned (provided values will be ignored).
            If ``badpix`` is not provided, then this parameter is ignored.

        Returns
        -------
        badpix: `numpy.ndarray`
            N x 3 numpy array, where N is the number of bad pixels.
            The coordinates of the bad pixels that were interpolated over are
            in the first two columns, and the corresponding pixel values in the
            third. If ``badpix`` was provided, this is the same as the input.
        goodpix: `numpy.ndarray`
            M x 3 numpy array, where M is the number of bad pixels.
            The coordinates of the good pixels that were used to construct the
            interpolant arein the first two columns, and the corresponding
            pixel values in the third. If ``goodpix`` was provided, the first
            two columns are same as the input, with the third column updated
            with the pixel values from the image plane of the ``maskedImage``.
        """

        if badpix is None or goodpix is None:
            badpix, goodpix = ctUtils.findGoodPixelsAroundBadPixels(
                maskedImage,
                self.config.badMaskPlanes,
                buffer=self.config.interpLength,
            )
        else:
            # Even if badpix and goodpix is provided, make sure to update
            # the values of goodpix.
            ctUtils.updateArrayFromImage(goodpix, maskedImage.image)

        # Construct the interpolant with goodpix.
        interpolator = CloughTocher2DInterpolator(
            list(zip(goodpix[:, 0], goodpix[:, 1])),
            goodpix[:, 2],
            fill_value=self.config.fillValue,
        )

        # Compute the interpolated values at bad pixel locations.
        badpix[:, 2] = interpolator(badpix[:, :2])

        # Fill in the bad pixels.
        ctUtils.updateImageFromArray(maskedImage.image, badpix)

        return badpix, goodpix
