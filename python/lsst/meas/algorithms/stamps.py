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
"""Collection of small images (stamps).
"""

__all__ = ["Stamp", "Stamps", "StampsBase", "writeFits", "readFitsWithOptions"]

from collections.abc import Sequence
import abc
from dataclasses import dataclass
import numpy

import lsst.afw.image as afwImage
import lsst.afw.fits as afwFits
from lsst.geom import Box2I, Point2I, Extent2I, Angle, degrees, SpherePoint
from lsst.daf.base import PropertyList


def writeFits(filename, stamp_ims, metadata, write_mask, write_variance):
    """Write a single FITS file containing all stamps.

    Parameters
    ----------
    filename : `str`
        A string indicating the output filename
    stamps_ims : iterable of `lsst.afw.image.MaskedImageF`
        An iterable of masked images
    metadata : `PropertyList`
        A collection of key, value metadata pairs to be
        written to the primary header
    write_mask : `bool`
        Write the mask data to the output file?
    write_variance : `bool`
        Write the variance data to the output file?
    """
    metadata['HAS_MASK'] = write_mask
    metadata['HAS_VARIANCE'] = write_variance
    metadata['N_STAMPS'] = len(stamp_ims)
    # create primary HDU with global metadata
    fitsPrimary = afwFits.Fits(filename, "w")
    fitsPrimary.createEmpty()
    fitsPrimary.writeMetadata(metadata)
    fitsPrimary.closeFile()

    # add all pixel data optionally writing mask and variance information
    for i, stamp in enumerate(stamp_ims):
        metadata = PropertyList()
        # EXTVER should be 1-based, the index from enumerate is 0-based
        metadata.update({'EXTVER': i+1, 'EXTNAME': 'IMAGE'})
        stamp.getImage().writeFits(filename, metadata=metadata, mode='a')
        if write_mask:
            metadata = PropertyList()
            metadata.update({'EXTVER': i+1, 'EXTNAME': 'MASK'})
            stamp.getMask().writeFits(filename, metadata=metadata, mode='a')
        if write_variance:
            metadata = PropertyList()
            metadata.update({'EXTVER': i+1, 'EXTNAME': 'VARIANCE'})
            stamp.getVariance().writeFits(filename, metadata=metadata, mode='a')
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
    options : `PropertyList`
        A collection of parameters.  If certain keys are available
        (``llcX``, ``llcY``, ``width``, ``height``), a bounding box
        is constructed and passed to the ``FitsReader`` in order
        to return a sub-image.

    Returns
    -------
    stamps : `list` of dataclass objects like `Stamp`, PropertyList
        A tuple of a list of `Stamp`-like objects
    metadata : `PropertyList`
        The metadata
    """
    # extract necessary info from metadata
    metadata = afwFits.readMetadata(filename, hdu=0)
    f = afwFits.Fits(filename, 'r')
    nExtensions = f.countHdus()
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
    stamp_parts = {}
    # We need to be careful because nExtensions includes the primary
    # header data unit
    for idx in range(nExtensions-1):
        md = afwFits.readMetadata(filename, hdu=idx+1)
        if md['EXTNAME'] in ('IMAGE', 'VARIANCE'):
            reader = afwImage.ImageFitsReader(filename, hdu=idx+1)
        elif md['EXTNAME'] == 'MASK':
            reader = afwImage.MaskFitsReader(filename, hdu=idx+1)
        else:
            raise ValueError(f"Unknown extension type: {md['EXTNAME']}")
        stamp_parts.setdefault(md['EXTVER'], {})[md['EXTNAME'].lower()] = reader.read(**kwargs)
    if len(stamp_parts) != nStamps:
        raise ValueError(f'Number of stamps read ({len(stamp_parts)}) does not agree with the '
                         f'number of stamps recorded in the metadata ({nStamps}).')
    # construct stamps themselves
    stamps = []
    for k in range(nStamps):
        # Need to increment by one since EXTVER starts at 1
        maskedImage = afwImage.MaskedImageF(**stamp_parts[k+1])
        stamps.append(stamp_factory(maskedImage, metadata, k))

    return stamps, metadata


@dataclass
class AbstractStamp(abc.ABC):
    """Single abstract stamp

    Parameters
    ----------
    Inherit from this class to add metadata to the stamp
    """

    @classmethod
    @abc.abstractmethod
    def factory(cls, stamp_im, metadata, index):
        """This method is needed to service the FITS reader.
        We need a standard interface to construct objects like this.
        Parameters needed to construct this object are passed in via
        a metadata dictionary and then passed to the constructor of
        this class.

        Parameters
        ----------
        stamp : `lsst.afw.image.MaskedImage`
            Pixel data to pass to the constructor
        metadata : `dict`
            Dictionary containing the information
            needed by the constructor.
        idx : `int`
            Index into the lists in ``metadata``

        Returns
        -------
        stamp : `AbstractStamp`
            An instance of this class
        """
        raise NotImplementedError


@dataclass
class Stamp(AbstractStamp):
    """Single stamp

    Parameters
    ----------
    stamp_im : `lsst.afw.image.MaskedImageF`
        The actual pixel values for the postage stamp
    position : `lsst.geom.SpherePoint`
        Position of the center of the stamp.  Note the user
        must keep track of the coordinate system
    """
    stamp_im: afwImage.maskedImage.MaskedImageF
    position: SpherePoint

    @classmethod
    def factory(cls, stamp_im, metadata, index):
        """This method is needed to service the FITS reader.
        We need a standard interface to construct objects like this.
        Parameters needed to construct this object are passed in via
        a metadata dictionary and then passed to the constructor of
        this class.  If lists of values are passed with the following
        keys, they will be passed to the constructor, otherwise dummy
        values will be passed: RA_DEG, DEC_DEG.  They should
        each point to lists of values.

        Parameters
        ----------
        stamp : `lsst.afw.image.MaskedImage`
            Pixel data to pass to the constructor
        metadata : `dict`
            Dictionary containing the information
            needed by the constructor.
        idx : `int`
            Index into the lists in ``metadata``

        Returns
        -------
        stamp : `Stamp`
            An instance of this class
        """
        if 'RA_DEG' in metadata and 'DEC_DEG' in metadata:
            return cls(stamp_im=stamp_im,
                       position=SpherePoint(Angle(metadata.getArray('RA_DEG')[index], degrees),
                                            Angle(metadata.getArray('DEC_DEG')[index], degrees)))
        else:
            return cls(stamp_im=stamp_im, position=SpherePoint(Angle(numpy.nan), Angle(numpy.nan)))


class StampsBase(abc.ABC, Sequence):
    """Collection of  stamps and associated metadata.

    Parameters
    ----------
    stamps : iterable
        This should be an iterable of dataclass objects
        a la ``lsst.meas.algorithms.Stamp``.
    metadata : `lsst.daf.base.PropertyList`, optional
        Metadata associated with the bright stars.
    use_mask : `bool`, optional
        If ``True`` read and write the mask data.  Default ``True``.
    use_variance : `bool`, optional
        If ``True`` read and write the variance data. Default ``True``.

    Notes
    -----
    A butler can be used to read only a part of the stamps,
    specified by a bbox:

    >>> starSubregions = butler.get("brightStarStamps_sub", dataId, bbox=bbox)
    """

    def __init__(self, stamps, metadata=None, use_mask=True, use_variance=True):
        if not hasattr(stamps, '__iter__'):
            raise ValueError('The stamps parameter must be iterable.')
        for stamp in stamps:
            if not isinstance(stamp, AbstractStamp):
                raise ValueError('The entries in stamps must inherit from AbstractStamp.  '
                                 f'Got {type(stamp)}.')
        self._stamps = stamps
        self._metadata = PropertyList() if metadata is None else metadata.deepCopy()
        self.use_mask = use_mask
        self.use_variance = use_variance

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
        options : `PropertyList`
            Collection of metadata parameters
        """
        raise NotImplementedError

    @abc.abstractmethod
    def _refresh_metadata(self):
        """Make sure metadata is up to date since this object
        can be extende
        """
        raise NotImplementedError

    def writeFits(self, filename):
        """Write this object to a file.

        Parameters
        ----------
        filename : `str`
            Name of file to write
        """
        self._refresh_metadata()
        stamps_ims = self.getMaskedImages()
        writeFits(filename, stamps_ims, self._metadata, self.use_mask, self.use_variance)

    def __len__(self):
        return len(self._stamps)

    def __getitem__(self, index):
        return self._stamps[index]

    def __iter__(self):
        return iter(self._stamps)

    def getMaskedImages(self):
        """Retrieve star images.

        Returns
        -------
        maskedImages :
            `list` [`lsst.afw.image.maskedImage.maskedImage.MaskedImageF`]
        """
        return [stamp.stamp_im for stamp in self._stamps]

    @property
    def metadata(self):
        return self._metadata


class Stamps(StampsBase):
    def _refresh_metadata(self):
        positions = self.getPositions()
        self._metadata['RA_DEG'] = [p.getRa().asDegrees() for p in positions]
        self._metadata['DEC_DEG'] = [p.getDec().asDegrees() for p in positions]

    def getPositions(self):
        return [s.position for s in self._stamps]

    def append(self, item):
        """Add an additional stamp.

        Parameters
        ----------
        item : `Stamp`
            Stamp object to append.
        """
        if not isinstance(item, Stamp):
            raise ValueError("Ojbects added must be a Stamp object.")
        self._stamps.append(item)
        return None

    def extend(self, stamp_list):
        """Extend Stamps instance by appending elements from another instance.

        Parameters
        ----------
        stamps_list : `list` [`Stamp`]
            List of Stamp object to append.
        """
        for s in stamp_list:
            if not isinstance(s, Stamp):
                raise ValueError('Can only extend with Stamp objects')
        self._stamps += stamp_list

    @classmethod
    def readFits(cls, filename):
        """Build an instance of this class from a file.

        Parameters
        ----------
        filename : `str`
            Name of the file to read

        Returns
        -------
        object : `Stamps`
            An instance of this class
        """
        return cls.readFitsWithOptions(filename, None)

    @classmethod
    def readFitsWithOptions(cls, filename, options):
        """Build an instance of this class with options.

        Parameters
        ----------
        filename : `str`
            Name of the file to read
        options : `PropertyList`
            Collection of metadata parameters

        Returns
        -------
        object : `Stamps`
            An instance of this class
        """
        stamps, metadata = readFitsWithOptions(filename, Stamp.factory, options)
        return cls(stamps, metadata=metadata, use_mask=metadata['HAS_MASK'],
                   use_variance=metadata['HAS_VARIANCE'])
