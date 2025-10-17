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

"""Collection of small images (postage stamps) centered on bright stars."""

from __future__ import annotations

__all__ = ["BrightStarStamp", "BrightStarStamps"]

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from lsst.afw.detection import Psf
from lsst.afw.fits import Fits, readMetadata
from lsst.afw.geom import SkyWcs
from lsst.afw.image import ImageFitsReader, MaskedImageF, MaskFitsReader
from lsst.afw.table.io import InputArchive, OutputArchive
from lsst.daf.base import PropertyList
from lsst.geom import Angle, Point2D, degrees
from lsst.meas.algorithms.stamps import AbstractStamp
from lsst.utils.introspection import get_full_type_name


@dataclass
class BrightStarStamp(AbstractStamp):
    """A single postage stamp centered on a bright star.

    Attributes
    ----------
    stamp_im : `~lsst.afw.image.MaskedImageF`
        The pixel data for this stamp.
    psf : `~lsst.afw.detection.Psf`, optional
        The point-spread function for this star.
    wcs : `~lsst.afw.geom.SkyWcs`, optional
        World coordinate system associated with the stamp.
    visit : `int`
        Visit number of the observation.
    detector : `int`
        Detector ID within the visit.
    ref_id : `int`
        Reference catalog ID of the star.
    ref_mag : `float`
        Reference catalog magnitude of the star.
    position : `~lsst.geom.Point2D`
        Center position of the star on the detector in pixel coordinates.
    focal_plane_radius : `float`
        Radial distance from the focal plane center in tangent-plane pixels.
    focal_plane_angle : `~lsst.geom.Angle`
        Azimuthal angle on the focal plane (counterclockwise from +X).
    scale : `float`, optional
        Flux scaling factor applied to the PSF model.
    scale_err : `float`, optional
        Error in the flux scale.
    pedestal : `float`, optional
        Background pedestal level.
    pedestal_err : `float`, optional
        Error on the pedestal.
    pedestal_scale_cov : `float`, optional
        Covariance between pedestal and scale.
    gradient_x : `float`, optional
        Background gradient in the X direction.
    gradient_y : `float`, optional
        Background gradient in the Y direction.
    global_reduced_chi_squared : `float`, optional
        Reduced chi-squared for the global model fit.
    global_degrees_of_freedom : `int`, optional
        Degrees of freedom for the global model fit.
    psf_reduced_chi_squared : `float`, optional
        Reduced chi-squared for the PSF fit.
    psf_degrees_of_freedom : `int`, optional
        Degrees of freedom for the PSF fit.
    psf_masked_flux_fraction : `float`, optional
        Fraction of flux masked in the PSF.

    Notes
    -----
    This class is designed to be used with `BrightStarStamps`, which manages
    collections of these stamps and handles reading/writing them to FITS files.
    The `factory` class method provides a standard interface to construct
    instances from image data and metadata, while the `_getMetadata` method
    extracts metadata for storage in FITS headers.
    """

    stamp_im: MaskedImageF
    psf: Psf | None
    wcs: SkyWcs | None
    visit: int
    detector: int
    ref_id: int
    ref_mag: float
    position: Point2D
    focal_plane_radius: float | None
    focal_plane_angle: Angle | None
    scale: float | None
    scale_err: float | None
    pedestal: float | None
    pedestal_err: float | None
    pedestal_scale_cov: float | None
    gradient_x: float | None
    gradient_y: float | None
    global_reduced_chi_squared: float | None
    global_degrees_of_freedom: int | None
    psf_reduced_chi_squared: float | None
    psf_degrees_of_freedom: int | None
    psf_masked_flux_fraction: float | None

    # Mapping of metadata keys to attribute names
    _metadata_attribute_map = {
        "VISIT": "visit",
        "DETECTOR": "detector",
        "REF_ID": "ref_id",
        "REF_MAG": "ref_mag",
        "POSITION_X": "position.x",
        "POSITION_Y": "position.y",
        "FOCAL_PLANE_RADIUS": "focal_plane_radius",
        "FOCAL_PLANE_ANGLE_DEGREES": "focal_plane_angle",
        "SCALE": "scale",
        "SCALE_ERR": "scale_err",
        "PEDESTAL": "pedestal",
        "PEDESTAL_ERR": "pedestal_err",
        "PEDESTAL_SCALE_COV": "pedestal_scale_cov",
        "GRADIENT_X": "gradient_x",
        "GRADIENT_Y": "gradient_y",
        "GLOBAL_REDUCED_CHI_SQUARED": "global_reduced_chi_squared",
        "GLOBAL_DEGREES_OF_FREEDOM": "global_degrees_of_freedom",
        "PSF_REDUCED_CHI_SQUARED": "psf_reduced_chi_squared",
        "PSF_DEGREES_OF_FREEDOM": "psf_degrees_of_freedom",
        "PSF_MASKED_FLUX_FRACTION": "psf_masked_flux_fraction",
    }

    def _getMetadata(self) -> PropertyList:
        """Extract metadata from the stamp's attributes.

        This method constructs a `PropertyList` containing metadata
        extracted from the stamp's attributes. It is used when writing the
        stamp to a FITS file to store relevant metadata in the FITS headers.

        Returns
        -------
        metadata : `PropertyList`
            A `PropertyList` containing the metadata, or `None` if no
            metadata attributes are defined.
        """
        metadata = PropertyList()
        for metadata_key, attribute_name in self._metadata_attribute_map.items():
            if "." in attribute_name:
                top_attr, sub_attr = attribute_name.split(".")
                value = getattr(getattr(self, top_attr), sub_attr)
            elif metadata_key == "FOCAL_PLANE_ANGLE_DEGREES":
                value = getattr(self, attribute_name).asDegrees()
            else:
                value = getattr(self, attribute_name)
            metadata[metadata_key] = value
        return metadata

    @property
    def metadata(self) -> PropertyList:
        """Return the stamp's metadata as a PropertyList."""
        return self._getMetadata()

    @classmethod
    def factory(
        cls,
        stamp_im: MaskedImageF,
        psf: Psf | None,
        wcs: SkyWcs | None,
        metadata: PropertyList,
    ) -> BrightStarStamp:
        """Construct a `BrightStarStamp` from image data and metadata.

        This method provides a standard interface to create a `BrightStarStamp`
        from its image data, PSF, WCS, and associated metadata.
        It is used by the `BrightStarStamps.readFits` method to construct
        individual bright star stamps from FITS files.

        Parameters
        ----------
        stamp_im : `~lsst.afw.image.MaskedImageF`
            Masked image for the stamp.
        psf : `~lsst.afw.detection.Psf`, optional
            Point-spread function for the stamp.
        wcs : `~lsst.afw.geom.SkyWcs`, optional
            World coordinate system for the stamp.
        metadata : `PropertyList`
            Metadata associated with the stamp, containing keys for all
            required attributes.

        Returns
        -------
        brightStarStamp : `BrightStarStamp`
            The constructed `BrightStarStamp` instance.
        """
        kwargs = {}

        for metadata_key, attribute_name in cls._metadata_attribute_map.items():
            if "." in attribute_name:  # for nested attributes like position.x
                top_attr, sub_attr = attribute_name.split(".")
                if top_attr not in kwargs:  # avoid overwriting position
                    if top_attr == "position":  # make an initial Point2D
                        kwargs[top_attr] = Point2D(0, 0)
                setattr(kwargs[top_attr], sub_attr, metadata[metadata_key])
            elif attribute_name == "focal_plane_angle":
                kwargs[attribute_name] = Angle(metadata[metadata_key], degrees)
            else:
                kwargs[attribute_name] = metadata[metadata_key]

        return cls(stamp_im=stamp_im, psf=psf, wcs=wcs, **kwargs)


class BrightStarStamps(Sequence[BrightStarStamp]):
    """A collection of bright star stamps.

    Parameters
    ----------
    brightStarStamps : `Iterable` [`BrightStarStamp`]
        Collection of `BrightStarStamp` instances.
    metadata : `~lsst.daf.base.PropertyList`, optional
        Global metadata associated with the collection.
    """

    def __init__(
        self,
        brightStarStamps: Sequence[BrightStarStamp],
        metadata: PropertyList | None = None,
    ):
        self._stamps = list(brightStarStamps)
        self._metadata = PropertyList() if metadata is None else metadata.deepCopy()
        self.by_ref_id = {stamp.ref_id: stamp for stamp in self}

    def __len__(self):
        return len(self._stamps)

    def __getitem__(self, index):
        if isinstance(index, slice):
            return BrightStarStamps(self._stamps[index], metadata=self._metadata)
        return self._stamps[index]

    def __iter__(self):
        return iter(self._stamps)

    @property
    def metadata(self):
        """Return the collection's global metadata as a PropertyList."""
        return self._metadata

    @classmethod
    def readFits(cls, filename: str) -> BrightStarStamps:
        """Make a `BrightStarStamps` object from a FITS file.

        Parameters
        ----------
        filename : `str`
            Name of the FITS file to read.

        Returns
        -------
        brightStarStamps : `BrightStarStamps`
            The constructed `BrightStarStamps` instance.
        """
        return cls.readFitsWithOptions(filename, None)

    @classmethod
    def readFitsWithOptions(cls, filename: str, options: PropertyList | None) -> BrightStarStamps:
        """Make a `BrightStarStamps` object from a FITS file, with options.

        Parameters
        ----------
        filename : `str`
            Name of the FITS file to read.
        options : `~lsst.daf.base.PropertyList`, optional
            Options for reading the FITS file. Not currently used.

        Returns
        -------
        brightStarStamps : `BrightStarStamps`
            The constructed `BrightStarStamps` instance.
        """
        with Fits(filename, "r") as fits_file:
            stamp_planes = {}
            stamp_psf_ids = {}
            stamp_wcs_ids = {}
            stamp_metadata = {}
            archive = None

            for hdu_num in range(1, fits_file.countHdus()):  # Skip primary HDU
                metadata = readMetadata(filename, hdu=hdu_num)
                extname = metadata["EXTNAME"]
                stamp_id: int | None = metadata.get("EXTVER", None)

                # Skip non-image BINTABLEs (except ARCHIVE_INDEX)
                if metadata["XTENSION"] == "BINTABLE" and not metadata.get("ZIMAGE", False):
                    if extname != "ARCHIVE_INDEX":
                        continue

                # Handle the archive index separately
                if extname == "ARCHIVE_INDEX":
                    fits_file.setHdu(hdu_num)
                    archive = InputArchive.readFits(fits_file)
                    continue
                elif metadata.get("EXTTYPE") == "ARCHIVE_DATA":
                    continue

                # Select reader and dtype
                if extname == "IMAGE":
                    reader = ImageFitsReader(filename, hdu=hdu_num)
                    dtype = np.dtype(MaskedImageF.dtype)
                    stamp_psf_ids[stamp_id] = metadata.pop("PSF", None)
                    stamp_wcs_ids[stamp_id] = metadata.pop("WCS", None)
                    stamp_metadata[stamp_id] = metadata
                elif extname == "MASK":
                    reader = MaskFitsReader(filename, hdu=hdu_num)
                    dtype = None
                elif extname == "VARIANCE":
                    reader = ImageFitsReader(filename, hdu=hdu_num)
                    dtype = np.dtype("float32")
                else:
                    raise ValueError(f"Unknown extension type: {extname}")

                if stamp_id is not None:
                    stamp_planes.setdefault(stamp_id, {})[extname.lower()] = reader.read(dtype=dtype)

        primary_metadata = readMetadata(filename, hdu=0)
        num_stamps = primary_metadata["N_STAMPS"]

        if len(stamp_planes) != num_stamps:
            raise ValueError(
                f"Number of stamps read ({len(stamp_planes)}) does not agree with the "
                f"number of stamps recorded in the primary HDU metadata ({num_stamps})."
            )
        if archive is None:
            raise ValueError("No archive index was found in the FITS file; cannot read PSF or WCS.")

        brightStarStamps = []
        for stamp_id in range(1, num_stamps + 1):  # Need to increment by one as EXTVER starts at 1
            stamp = MaskedImageF(**stamp_planes[stamp_id])
            psf = archive.get(stamp_psf_ids[stamp_id])
            wcs = archive.get(stamp_wcs_ids[stamp_id])
            brightStarStamps.append(BrightStarStamp.factory(stamp, psf, wcs, stamp_metadata[stamp_id]))

        return cls(brightStarStamps, primary_metadata)

    def writeFits(self, filename: str):
        """Write this `BrightStarStamps` object to a FITS file.

        Parameters
        ----------
        filename : `str`
            Name of the FITS file to write.
        """
        metadata = self._metadata.deepCopy()

        # Store metadata in the primary HDU
        metadata["N_STAMPS"] = len(self._stamps)
        metadata["VERSION"] = 2  # Record version number in case of future code changes
        metadata["STAMPCLS"] = get_full_type_name(self)

        # Create and write to the FITS file within a context manager
        with Fits(filename, "w") as fits_file:
            fits_file.createEmpty()

            # Store Persistables in an OutputArchive
            output_archive = OutputArchive()
            stamp_psf_ids = []
            stamp_wcs_ids = []
            for stamp in self._stamps:
                stamp_psf_ids.append(output_archive.put(stamp.psf))
                stamp_wcs_ids.append(output_archive.put(stamp.wcs))

            # Write to the FITS file
            fits_file.writeMetadata(metadata)
            del metadata
            output_archive.writeFits(fits_file)

        # Add all pixel data to extension HDUs; note: EXTVER should be 1-based
        for stamp_id, (stamp, stamp_psf_id, stamp_wcs_id) in enumerate(
            zip(self._stamps, stamp_psf_ids, stamp_wcs_ids),
            start=1,
        ):
            metadata = PropertyList()
            metadata.update({"EXTVER": stamp_id, "EXTNAME": "IMAGE"})
            if stamp_metadata := stamp._getMetadata():
                metadata.update(stamp_metadata)
            metadata["PSF"] = stamp_psf_id
            metadata["WCS"] = stamp_wcs_id
            stamp.stamp_im.getImage().writeFits(filename, metadata=metadata, mode="a")

            metadata = PropertyList()
            metadata.update({"EXTVER": stamp_id, "EXTNAME": "MASK"})
            stamp.stamp_im.getMask().writeFits(filename, metadata=metadata, mode="a")

            metadata = PropertyList()
            metadata.update({"EXTVER": stamp_id, "EXTNAME": "VARIANCE"})
            stamp.stamp_im.getVariance().writeFits(filename, metadata=metadata, mode="a")
