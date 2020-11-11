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

__all__ = ["BrightStarStamp", "BrightStarStamps"]

import collections.abc
from dataclasses import dataclass

import lsst.afw.image as afwImage
import lsst.afw.fits as afwFits
from lsst.geom import Box2I, Point2I, Extent2I, Angle
from lsst.daf.base import PropertySet


class Stamp(dataclass):
    """Single stamp 
    """
    stamp: afwImage.maskedImage.MaskedImageF
    ra: Angle
    dec: Angle
    size: int

    @classmethod
    def factory(cls, stamp, metadata, index):
        if 'RA_DEG' in metadata.keys() and 'DEC_DEG' in metadata.keys() and 'SIZE' in metadata.keys():
            return cls(stamp, Angle(metadata['RA_DEG'][index]), Angle(metadata['DEC_DEG'][index]),
                                    metadata['SIZE'][index])
        else:
            return cls(stamp=stamp, ra=Angle(numpy.nan), dec=Angle(numpy.nan), size=-1)

class Stamps(collections.abc.Sequence):
    """Collection of  stamps and associated metadata.

    Parameters
    ----------
    stamps : `collections.abc.Sequence` [`Stamp`]
        Sequence of stamps.
    metadata : `lsst.daf.base.PropertyList`, optional
        Metadata associated with the bright stars.

    Raises
    ------
    ValueError
        Raises if the stamps are not the correct type

    Notes
    -----
    A (gen2) butler can be used to read only a part of the stamps,
    specified by a bbox:

    >>> starSubregions = butler.get("brightStarStamps_sub", dataId, bbox=bbox)
    """

    def __init__(self, stamps, metadata=None, has_mask=True, has_variance=True):
        for item in stamps:
            if not isinstance(item, Stamp):
                raise ValueError(f"Can only add instances of Stamp, got {type(item)}")
        self._stamps = stamps
        self._metadata = PropertySet() if metadata is None else metadata.deepCopy()
        self._metadata['HAS_MASK'] = has_variance
        self._metadata['HAS_VARIANCE'] = has_variance

    def __len__(self):
        return len(self._stamps)

    def __getitem__(self, index):
        return self._stamps[index]

    def __iter__(self):
        return iter(self._stamps)

    def _refresh_metadata(self):
        # ensure metadata contains current number of objects
        self._metadata["N_STAMPS"] = len(self)

    def append(self, item):
        """Add an additional bright star stamp.

        Parameters
        ----------
        item : `Stamp`
            Stamp to append.
        """
        if not isinstance(item, Stamp):
            raise ValueError(f"Can only add instances of Stamp, got {type(item)}.")
        self._stamps.append(item)
        return None

    def extend(self, s):
        """Extend Stamps instance by appending elements from another instance.

        Parameters
        ----------
        s : `Stamps`
            Other instance to concatenate.
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

    def writeFits(self, filename):
        """Write a single FITS file containing all stamps.
        """
        # Make sure we have up to date metadata
        self.refresh_metadata()
        # create primary HDU with global metadata
        fitsPrimary = afwFits.Fits(filename, "w")
        fitsPrimary.createEmpty()
        fitsPrimary.writeMetadata(self._metadata)
        fitsPrimary.closeFile()

        # add all stamps and mask planes
        for stamp in self.getMaskedImages():
            stamp.getImage().writeFits(filename, mode='a')
            if self._metadata['HAS_MASK']:
                stamp.getMask().writeFits(filename, mode='a')
            if self._metadata['HAS_VARIANCE']:
                stamp.getVariance().writeFits(filename, mode='a')
        return None

    @classmethod
    def readFits(cls, filename):
        """Read bright star stamps from FITS file.

        Returns
        -------
        stamps : `Stamps`
            Collection of stamps.
        """
        return cls.readFitsWithOptions(filename, None)

    @classmethod
    def readFitsWithOptions(cls, filename, options):
        """Read stamps from FITS file, allowing for only a
        subregion of the stamps to be read.

        Returns
        -------
        stamps : `Stamps`
            Collection of stamps.
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
        starStamps = []
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
                varReader = afwImage.VarianceFitsReader(filename, hdu=n_parts*idx + n_parts)
                parts['variance'] = varReader.read(**kwargs)
            maskedImage = afwImage.MaskedImageF(**parts)
            starStamps.append(Stamp.factory(maskedImage, metadata, idx))
        return cls(starStamps, metadata=visitMetadata)
