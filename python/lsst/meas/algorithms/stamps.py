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

"""Collection of small images (stamps)."""

__all__ = ["Stamp", "Stamps", "StampsBase", "writeFits", "readFitsWithOptions"]

import abc
from collections.abc import Sequence
import copy
from dataclasses import dataclass, field, fields

import numpy as np
from lsst.afw.fits import Fits, readMetadata
from lsst.afw.image import ImageFitsReader, MaskedImage, MaskedImageF, MaskFitsReader
from lsst.afw.table.io import InputArchive, OutputArchive, Persistable
from lsst.daf.base import PropertyList
from lsst.geom import Angle, Box2I, Extent2I, Point2I, SpherePoint, degrees
from lsst.utils import doImport
from lsst.utils.introspection import get_full_type_name


def writeFits(filename, stamps, metadata, type_name, write_mask, write_variance, write_archive=False):
    """Write a single FITS file containing all stamps.

    Parameters
    ----------
    filename : `str`
        A string indicating the output filename
    stamps : iterable of `BaseStamp`
        An iterable of Stamp objects
    metadata : `PropertyList`
        A collection of key, value metadata pairs to be
        written to the primary header
    type_name : `str`
        Python type name of the StampsBase subclass to use
    write_mask : `bool`
        Write the mask data to the output file?
    write_variance : `bool`
        Write the variance data to the output file?
    write_archive : `bool`, optional
        Write an archive to store Persistables along with each stamp?
        Default: ``False``.
    """
    metadata["HAS_MASK"] = write_mask
    metadata["HAS_VARIANCE"] = write_variance
    metadata["HAS_ARCHIVE"] = write_archive
    metadata["N_STAMPS"] = len(stamps)
    metadata["STAMPCLS"] = type_name
    # Record version number in case of future code changes
    metadata["VERSION"] = 1
    # create primary HDU with global metadata
    fitsFile = Fits(filename, "w")
    fitsFile.createEmpty()
    # Store Persistables in an OutputArchive and write it
    if write_archive:
        oa = OutputArchive()
        archive_ids = [oa.put(stamp.archive_element) for stamp in stamps]
        metadata["ARCHIVE_IDS"] = archive_ids
        fitsFile.writeMetadata(metadata)
        oa.writeFits(fitsFile)
    else:
        fitsFile.writeMetadata(metadata)
    fitsFile.closeFile()
    # add all pixel data optionally writing mask and variance information
    for i, stamp in enumerate(stamps):
        metadata = PropertyList()
        # EXTVER should be 1-based, the index from enumerate is 0-based
        metadata.update({"EXTVER": i + 1, "EXTNAME": "IMAGE"})
        stamp.stamp_im.getImage().writeFits(filename, metadata=metadata, mode="a")
        if write_mask:
            metadata = PropertyList()
            metadata.update({"EXTVER": i + 1, "EXTNAME": "MASK"})
            stamp.stamp_im.getMask().writeFits(filename, metadata=metadata, mode="a")
        if write_variance:
            metadata = PropertyList()
            metadata.update({"EXTVER": i + 1, "EXTNAME": "VARIANCE"})
            stamp.stamp_im.getVariance().writeFits(filename, metadata=metadata, mode="a")
    return None


def readFitsWithOptions(filename, stamp_factory, options):
    """Read stamps from FITS file, allowing for only a subregion of the stamps
    to be read.

    Parameters
    ----------
    filename : `str`
        A string indicating the file to read
    stamp_factory : classmethod
        A factory function defined on a dataclass for constructing
        stamp objects a la `~lsst.meas.algorithm.Stamp`
    options : `PropertyList` or `dict`
        A collection of parameters. If it contains a bounding box
        (``bbox`` key), or if certain other keys (``llcX``, ``llcY``,
        ``width``, ``height``) are available for one to be constructed,
        the bounding box is passed to the ``FitsReader`` in order to
        return a sub-image.

    Returns
    -------
    stamps : `list` of dataclass objects like `Stamp`, PropertyList
        A tuple of a list of `Stamp`-like objects
    metadata : `PropertyList`
        The metadata

    Notes
    -----
    The data are read using the data type expected by the
    `~lsst.afw.image.MaskedImage` class attached to the `AbstractStamp`
    dataclass associated with the factory method.
    """
    # extract necessary info from metadata
    metadata = readMetadata(filename, hdu=0)
    nStamps = metadata["N_STAMPS"]
    has_archive = metadata["HAS_ARCHIVE"]
    if has_archive:
        archive_ids = metadata.getArray("ARCHIVE_IDS")
    with Fits(filename, "r") as f:
        nExtensions = f.countHdus()
        # check if a bbox was provided
        kwargs = {}
        if options:
            # gen3 API
            if "bbox" in options.keys():
                kwargs["bbox"] = options["bbox"]
            # gen2 API
            elif "llcX" in options.keys():
                llcX = options["llcX"]
                llcY = options["llcY"]
                width = options["width"]
                height = options["height"]
                bbox = Box2I(Point2I(llcX, llcY), Extent2I(width, height))
                kwargs["bbox"] = bbox
        stamp_parts = {}

        # Determine the dtype from the factory. This allows a Stamp class
        # to be defined in terms of MaskedImageD or MaskedImageI without
        # forcing everything to floats.
        masked_image_cls = None
        for stamp_field in fields(stamp_factory.__self__):
            if issubclass(stamp_field.type, MaskedImage):
                masked_image_cls = stamp_field.type
                break
        else:
            raise RuntimeError("Stamp factory does not use MaskedImage.")
        default_dtype = np.dtype(masked_image_cls.dtype)
        variance_dtype = np.dtype(np.float32)  # Variance is always the same type.

        # We need to be careful because nExtensions includes the primary HDU.
        stamp_metadata = []
        for idx in range(nExtensions - 1):
            dtype = None
            _metadata = copy.copy(metadata)
            md = readMetadata(filename, hdu=idx + 1)
            _metadata.update(md)
            stamp_metadata.append(_metadata)
            # Skip binary tables that aren't images or archives.
            if md["XTENSION"] == "BINTABLE" and not ("ZIMAGE" in md and md["ZIMAGE"]):
                if md["EXTNAME"] != "ARCHIVE_INDEX":
                    continue
            if md["EXTNAME"] in ("IMAGE", "VARIANCE"):
                reader = ImageFitsReader(filename, hdu=idx + 1)
                if md["EXTNAME"] == "VARIANCE":
                    dtype = variance_dtype
                else:
                    dtype = default_dtype
            elif md["EXTNAME"] == "MASK":
                reader = MaskFitsReader(filename, hdu=idx + 1)
            elif md["EXTNAME"] == "ARCHIVE_INDEX":
                f.setHdu(idx + 1)
                archive = InputArchive.readFits(f)
                continue
            elif md["EXTTYPE"] == "ARCHIVE_DATA":
                continue
            else:
                raise ValueError(f"Unknown extension type: {md['EXTNAME']}")
            stamp_parts.setdefault(md["EXTVER"], {})[md["EXTNAME"].lower()] = reader.read(dtype=dtype,
                                                                                          **kwargs)
    if len(stamp_parts) != nStamps:
        raise ValueError(
            f"Number of stamps read ({len(stamp_parts)}) does not agree with the "
            f"number of stamps recorded in the metadata ({nStamps})."
        )
    # construct stamps themselves
    stamps = []
    for k in range(nStamps):
        # Need to increment by one since EXTVER starts at 1
        maskedImage = masked_image_cls(**stamp_parts[k + 1])
        archive_element = archive.get(archive_ids[k]) if has_archive else None
        stamps.append(stamp_factory(maskedImage, stamp_metadata[k], k, archive_element))

    return stamps, metadata


@dataclass
class AbstractStamp(abc.ABC):
    """Single abstract stamp.

    Parameters
    ----------
    Inherit from this class to add metadata to the stamp.
    """

    @classmethod
    @abc.abstractmethod
    def factory(cls, stamp_im, metadata, index, archive_element=None):
        """This method is needed to service the FITS reader. We need a standard
        interface to construct objects like this. Parameters needed to
        construct this object are passed in via a metadata dictionary and then
        passed to the constructor of this class.

        Parameters
        ----------
        stamp : `~lsst.afw.image.MaskedImage`
            Pixel data to pass to the constructor
        metadata : `dict`
            Dictionary containing the information
            needed by the constructor.
        idx : `int`
            Index into the lists in ``metadata``
        archive_element : `~lsst.afw.table.io.Persistable`, optional
            Archive element (e.g. Transform or WCS) associated with this stamp.

        Returns
        -------
        stamp : `AbstractStamp`
            An instance of this class
        """
        raise NotImplementedError


def _default_position():
    # SpherePoint is nominally mutable in C++ so we must use a factory
    # and return an entirely new SpherePoint each time a Stamps is created.
    return SpherePoint(Angle(np.nan), Angle(np.nan))


@dataclass
class Stamp(AbstractStamp):
    """Single stamp.

    Parameters
    ----------
    stamp_im : `~lsst.afw.image.MaskedImageF`
        The actual pixel values for the postage stamp.
    archive_element : `~lsst.afw.table.io.Persistable` or `None`, optional
        Archive element (e.g. Transform or WCS) associated with this stamp.
    position : `~lsst.geom.SpherePoint` or `None`, optional
        Position of the center of the stamp. Note the user must keep track of
        the coordinate system.
    """

    stamp_im: MaskedImageF
    archive_element: Persistable | None = None
    position: SpherePoint | None = field(default_factory=_default_position)
    metadata: PropertyList | None = None

    @classmethod
    def factory(cls, stamp_im, metadata, index, archive_element=None):
        """This method is needed to service the FITS reader. We need a standard
        interface to construct objects like this. Parameters needed to
        construct this object are passed in via a metadata dictionary and then
        passed to the constructor of this class. If lists of values are passed
        with the following keys, they will be passed to the constructor,
        otherwise dummy values will be passed: RA_DEG, DEC_DEG. They should
        each point to lists of values.

        Parameters
        ----------
        stamp : `~lsst.afw.image.MaskedImage`
            Pixel data to pass to the constructor
        metadata : `PropertyList`
            A collection of key, value metadata pairs containing the
            information needed by the constructor.
        idx : `int`
            Index into the lists in ``metadata``
        archive_element : `~lsst.afw.table.io.Persistable`, optional
            Archive element (e.g. Transform or WCS) associated with this stamp.

        Returns
        -------
        stamp : `Stamp`
            An instance of this class
        """
        if "RA_DEG" in metadata and "DEC_DEG" in metadata:
            return cls(
                stamp_im=stamp_im,
                archive_element=archive_element,
                position=SpherePoint(
                    Angle(metadata.getArray("RA_DEG")[index], degrees),
                    Angle(metadata.getArray("DEC_DEG")[index], degrees),
                ),
                metadata=metadata,
            )
        else:
            return cls(
                stamp_im=stamp_im,
                archive_element=archive_element,
                position=SpherePoint(Angle(np.nan), Angle(np.nan)),
                metadata=metadata,
            )


class StampsBase(abc.ABC, Sequence):
    """Collection of stamps and associated metadata.

    Parameters
    ----------
    stamps : iterable
        This should be an iterable of dataclass objects
        a la ``~lsst.meas.algorithms.Stamp``.
    metadata : `~lsst.daf.base.PropertyList`, optional
        Metadata associated with the objects within the stamps.
    use_mask : `bool`, optional
        If ``True`` read and write the mask data. Default ``True``.
    use_variance : `bool`, optional
        If ``True`` read and write the variance data. Default ``True``.
    use_archive : `bool`, optional
        If ``True``, read and write an Archive that contains a Persistable
        associated with each stamp, for example a Transform or a WCS.
        Default ``False``.

    Notes
    -----
    A butler can be used to read only a part of the stamps,
    specified by a bbox:

    >>> starSubregions = butler.get(
            "brightStarStamps",
            dataId,
            parameters={"bbox": bbox}
        )
    """

    def __init__(self, stamps, metadata=None, use_mask=True, use_variance=True, use_archive=False):
        for stamp in stamps:
            if not isinstance(stamp, AbstractStamp):
                raise ValueError(f"The entries in stamps must inherit from AbstractStamp. Got {type(stamp)}.")
        self._stamps = stamps
        self._metadata = PropertyList() if metadata is None else metadata.deepCopy()
        self.use_mask = use_mask
        self.use_variance = use_variance
        self.use_archive = use_archive

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
        options : `PropertyList`
            Collection of metadata parameters
        """
        # To avoid problems since this is no longer an abstract method.
        # TO-DO: Consider refactoring this method. This class check was added
        # to allow the butler formatter to use a generic type but still end up
        # giving the correct type back, ensuring that the abstract base class
        # is not used by mistake. Perhaps this logic can be optimised.
        if cls is not StampsBase:
            raise NotImplementedError(f"Please implement specific FITS reader for class {cls}")

        # Load metadata to get class
        metadata = readMetadata(filename, hdu=0)
        type_name = metadata.get("STAMPCLS")
        if type_name is None:
            raise RuntimeError(
                f"No class name in file {filename}. Unable to instantiate correct stamps subclass. "
                "Is this an old version format Stamps file?"
            )

        # Import class and override `cls`
        stamp_type = doImport(type_name)
        cls = stamp_type

        return cls.readFitsWithOptions(filename, options)

    @abc.abstractmethod
    def _refresh_metadata(self):
        """Make sure metadata is up to date, as this object can be extended."""
        raise NotImplementedError

    def writeFits(self, filename):
        """Write this object to a file.

        Parameters
        ----------
        filename : `str`
            Name of file to write.
        """
        self._refresh_metadata()
        type_name = get_full_type_name(self)
        writeFits(
            filename,
            self._stamps,
            self._metadata,
            type_name,
            self.use_mask,
            self.use_variance,
            self.use_archive,
        )

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
            `list` [`~lsst.afw.image.MaskedImageF`]
        """
        return [stamp.stamp_im for stamp in self._stamps]

    def getArchiveElements(self):
        """Retrieve archive elements associated with each stamp.

        Returns
        -------
        archiveElements :
            `list` [`~lsst.afw.table.io.Persistable`]
        """
        return [stamp.archive_element for stamp in self._stamps]

    @property
    def metadata(self):
        return self._metadata


class Stamps(StampsBase):
    def _refresh_metadata(self):
        positions = self.getPositions()
        self._metadata["RA_DEG"] = [p.getRa().asDegrees() for p in positions]
        self._metadata["DEC_DEG"] = [p.getDec().asDegrees() for p in positions]

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
            raise ValueError("Objects added must be a Stamp object.")
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
                raise ValueError("Can only extend with Stamp objects")
        self._stamps += stamp_list

    @classmethod
    def readFits(cls, filename):
        """Build an instance of this class from a file.

        Parameters
        ----------
        filename : `str`
            Name of the file to read.

        Returns
        -------
        object : `Stamps`
            An instance of this class.
        """
        return cls.readFitsWithOptions(filename, None)

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
        stamps, metadata = readFitsWithOptions(filename, Stamp.factory, options)
        return cls(
            stamps,
            metadata=metadata,
            use_mask=metadata["HAS_MASK"],
            use_variance=metadata["HAS_VARIANCE"],
            use_archive=metadata["HAS_ARCHIVE"],
        )
