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
#
"""Collection of small images (stamps), each centered on a bright star.
"""

__all__ = ["Stamp", "Stamps", "StampsBase", "writeFits", "readFitsWithOptions"]

from collections.abc import Sequence
import abc
from dataclasses import dataclass
import numpy

import lsst.afw.image as afwImage
import lsst.afw.fits as afwFits
from lsst.geom import Box2I, Point2I, Extent2I, Angle, degrees
from lsst.daf.base import PropertySet


def writeFits(filename, stamps, metadata):
    """Write a single FITS file containing all stamps.

    Parameters
    ----------
    filename : `str`
        A string indicating the output filename
    stamps : iterable of `lsst.afw.image.MaskedImageF`
        An iterable of maked images
    metadata : `PropertySet`
        A collection of key, value metadata pairs
    """
    # create primary HDU with global metadata
    fitsPrimary = afwFits.Fits(filename, "w")
    fitsPrimary.createEmpty()
    fitsPrimary.writeMetadata(metadata)
    fitsPrimary.closeFile()

    # add all stamps and mask planes
    for stamp in stamps:
        stamp.getImage().writeFits(filename, mode='a')
        if metadata['HAS_MASK']:
            stamp.getMask().writeFits(filename, mode='a')
        if metadata['HAS_VARIANCE']:
            stamp.getVariance().writeFits(filename, mode='a')
    return None


def readFitsWithOptions(filename, stamp_factory, options):
    """Read stamps from FITS file, allowing for only a
    subregion of the stamps to be read.

    Parameters
    ----------
    filename : `str`
        A string indicating the file to read
    stamp_factory : classmethod
        A factory function defined on a dataclass for constructing
        stamp objects a la `lsst.meas.alrogithm.Stamp`
    options : `PropertySet`
        A collection of parameters.  If certain keys are available
        (``llcX``, ``llcY``, ``width``, ``height``), a bounding box
        is constructed and passed to the ``FitsReader`` in order
        to return a sub-image.

    Returns
    -------
    stamps, metadata : `list` of `lsst.afw.image.MaskedImageF`, PropertySet
        A tuple of a list of masked images and a collection of metadata.
    """
    # extract necessary info from metadata
    metadata = afwFits.readMetadata(filename, hdu=0)
    nStamps = metadata["N_STAMPS"]
    # check if a bbox was provided
    kwargs = {}
    if options and options.exists("llcX"):
        llcX = options["llcX"]
        llcY = options["llcY"]
        width = options["width"]
        height = options["height"]
        bbox = Box2I(Point2I(llcX, llcY), Extent2I(width, height))
        kwargs["bbox"] = bbox
    # read stamps themselves
    stamps = []
    for idx in range(nStamps):
        n_parts = 1  # Always have at least an image
        for key in ['HAS_MASK', 'HAS_VARIANCE']:
            if metadata[key]:
                n_parts += 1
        parts = {}
        imReader = afwImage.ImageFitsReader(filename, hdu=n_parts*idx + 1)
        parts['image'] = imReader.read(**kwargs)
        if metadata['HAS_MASK']:
            # Alwasy stored right after image if present
            maskReader = afwImage.MaskFitsReader(filename, hdu=n_parts*idx + 2)
            parts['mask'] = maskReader.read(**kwargs)
        if metadata['HAS_VARIANCE']:
            # Maybe either right after image or two after
            varReader = afwImage.ImageFitsReader(filename, hdu=n_parts*idx + n_parts)
            parts['variance'] = varReader.read(**kwargs)
        maskedImage = afwImage.MaskedImageF(**parts)
        stamps.append(stamp_factory(maskedImage, metadata.toDict(), idx))
    return stamps, metadata


@dataclass
class Stamp:
    """Single stamp

    Parameters
    ----------
    stamp : `lsst.afw.image.MaskedImageF`
        The actual pixel values for the postage stamp
    ra : `lsst.geom.Angle`
        The Right Ascention of the center of the stamp
    dec : `lsst.geom.Angle`
        The Declination of the center of the stamp
    size : `int`
        The size of the stamp in pixels
    """
    stamp: afwImage.maskedImage.MaskedImageF
    ra: Angle
    dec: Angle
    size: int

    @classmethod
    def factory(cls, stamp, metadata, index):
        """A factory method to construct an instance of this class
        given a masked image and some metadata.

        Parameteres
        -----------
        stamp : `lsst.afw.image.MaskedImageF`
            The pixel data to associate with the instance
        metadata : `dict`
            A mapping of metadata used to populate this instance
        index : `int`
            The index value (0 based) for looking up values from ``metadata``

        Returns
        -------
        stamp : `lsst.meas.algorithms.Stamp`
            An instance of this class
        """
        if 'RA_DEG' in metadata.keys() and 'DEC_DEG' in metadata.keys() and 'SIZE' in metadata.keys():
            return cls(stamp,
                       Angle(metadata['RA_DEG'][index], degrees),
                       Angle(metadata['DEC_DEG'][index], degrees),
                       metadata['SIZE'][index])
        else:
            return cls(stamp=stamp, ra=Angle(numpy.nan), dec=Angle(numpy.nan), size=-1)


class StampsBase(abc.ABC, Sequence):
    """Collection of  stamps and associated metadata.

    Parameters
    ----------
    stamps : iterable
        This should be an iterable of dataclass objects
        a la ``lsst.meas.algorithms.Stamp``.
    metadata : `lsst.daf.base.PropertyList`, optional
        Metadata associated with the bright stars.
    has_mask : `bool`, optional
        If ``True`` read and write the mask data.  Default ``True``.
    has_variance : `bool`, optional
        If ``True`` read and write the variance data. Default ``True``.

    Notes
    -----
    A (gen2) butler can be used to read only a part of the stamps,
    specified by a bbox:

    >>> starSubregions = butler.get("brightStarStamps_sub", dataId, bbox=bbox)
    """

    def __init__(self, stamps, metadata=None, has_mask=True, has_variance=True):
        if not hasattr(stamps, '__iter__'):
            raise ValueError('The stamps parameter must be iterable.')
        self._stamps = stamps
        self._metadata = PropertySet() if metadata is None else metadata.deepCopy()
        self._metadata['HAS_MASK'] = has_mask
        self._metadata['HAS_VARIANCE'] = has_variance

    @abc.abstractmethod
    def _refresh_metadata(self):
        """Refresh the metadata.  Should be called before writing this object out.
        """
        raise NotImplementedError

    @classmethod
    @abc.abstractmethod
    def readFits(cls, filename):
        """Build an instance of this class from a file.

        Parameters
        ----------
        filename : `str`
            Name of the file to read
        """
        raise NotImplementedError

    @classmethod
    @abc.abstractmethod
    def readFitsWithOptions(cls, filename, options):
        """Build an instance of this class with options.

        Parameters
        ----------
        filename : `str`
            Name of the file to read
        options : `PropertySet`
            Collection of metadata parameters
        """
        raise NotImplementedError

    @abc.abstractmethod
    def writeFits(self, filename):
        """Write this object to a file.

        Parameters
        ----------
        filename : `str`
            Name of file to write
        """
        raise NotImplementedError

    def __len__(self):
        return len(self._stamps)

    def __getitem__(self, index):
        return self._stamps[index]

    def __iter__(self):
        return iter(self._stamps)

    def append(self, item):
        """Add an additional bright star stamp.

        Parameters
        ----------
        item : `object`
            Stamp-like object to append.
        """
        if not hasattr(item, 'stamp'):
            raise ValueError("Ojbects added must contain a stamp attribute.")
        self._stamps.append(item)
        return None

    def extend(self, s):
        """Extend Stamps instance by appending elements from another instance.

        Parameters
        ----------
        s : `list` [`object`]
            List of Stamp-like object to append.
        """
        self._stamps += s._stamps

    def getMaskedImages(self):
        """Retrieve star images.

        Returns
        -------
        maskedImages :
            `list` [`lsst.afw.image.maskedImage.maskedImage.MaskedImageF`]
        """
        return [stamp.stamp for stamp in self._stamps]

    @property
    def metadata(self):
        return self._metadata.deepCopy()


class Stamps(StampsBase):
    def _refresh_metadata(self):
        """Refresh the metadata.  Should be called before writing this object out.
        """
        # ensure metadata contains current number of objects
        self._metadata["N_STAMPS"] = len(self)

    @classmethod
    def readFits(cls, filename):
        """Build an instance of this class from a file.

        Parameters
        ----------
        filename : `str`
            Name of the file to read
        """
        return cls.readFitsWithOptions(filename, None)

    @classmethod
    def readFitsWithOptions(cls, filename, options):
        """Build an instance of this class with options.

        Parameters
        ----------
        filename : `str`
            Name of the file to read
        options : `PropertySet`
            Collection of metadata parameters
        """
        stamps, metadata = readFitsWithOptions(filename, Stamp.factory, options)
        return cls(stamps, metadata=metadata, has_mask=metadata['HAS_MASK'],
                   has_variance=metadata['HAS_VARIANCE'])

    def writeFits(self, filename):
        """Write this object to a file.

        Parameters
        ----------
        filename : `str`
            Name of file to write
        """
        self._refresh_metadata()
        stamps = self.getMaskedImages()
        writeFits(filename, stamps, self._metadata)
