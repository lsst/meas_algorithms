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

"""Collection of image planes, derived from the Stamps class."""

__all__ = ["ImagePlane", "MultiImage"]

import numpy as np

from dataclasses import dataclass, field

from lsst.afw.image import ImageF
from lsst.afw.table.io import Persistable
from lsst.geom import SpherePoint, Angle

from .stamps import AbstractStamp, Stamps, _default_position, readFitsWithOptions


@dataclass
class ImagePlane(AbstractStamp):
    """A single image plane.

    Parameters
    ----------
    stamp_im : `~lsst.afw.image.ImageF`
        The actual pixel values for the postage stamp.
    archive_element : `~lsst.afw.table.io.Persistable` or `None`, optional
        An archive element (e.g. Transform or WCS) associated with
        this stamp.
    position : `~lsst.geom.SpherePoint` or `None`, optional
        Position of the center of the stamp.  Note the user must keep
        track of the coordinate system.
    """

    stamp_im: ImageF
    archive_element: Persistable | None = None
    position: SpherePoint | None = None

    @classmethod
    def factory(cls, stamp_im, metadata, index, archive_element=None):
        """This method is needed to service the FITS reader.  We need a
        standard interface to construct objects like this. Parameters
        needed to construct this object are passed in via a metadata
        dictionary and then passed to the constructor of this
        class. If lists of values are passed with the following keys,
        they will be passed to the constructor, otherwise dummy values
        will be passed: RA_DEG, DEC_DEG. They should each point to
        lists of values.

        Parameters
        ----------
        stamp_im : `~lsst.afw.image.ImageF`
            Pixel data to pass to the constructor.
        metadata : `dict`
            Dictionary containing the information needed by the
            constructor.
        index : `int`
            Index into the lists in ``metadata``
        archive_element : `~lsst.afw.table.io.Persistable`, optional
            An archive element (e.g. Transform or WCS) associated with
            this stamp.

        Returns
        -------
        imagePlane : `ImagePlane`
            An instance of this class.
        """
        return cls(
            stamp_im=stamp_im,
            archive_element=archive_element,
            position=None,
        )


class MultiImage(Stamps):
    """Collection of image planes and associated metadata.

    Parameters
    ----------
    imagePlanes : `collections.abc.Sequence` [`ImagePlane`]
        Sequence of image planes.
    """

    def _refresh_metadata(self):
        pass

    def getPositions(self):
        """Image planes do not have positions."""
        return None

    @classmethod
    def readFitsWithOptions(cls, filename, options):
        """Build an instance of this class with options.

        Parameters
        ----------
        filename : `str`
            Name of the file to read.
        options : `PropertyList` or `dict`
            Collection of metadata parameters.

        Returns
        -------
        object : `Stamps`
            An instance of this class.
        """
        stamps, metadata = readFitsWithOptions(filename, ImagePlane.factory, options)
        return cls(
            stamps,
            metadata=metadata,
            use_mask=False,
            use_variance=False,
            use_archive=metadata["HAS_ARCHIVE"],
        )
