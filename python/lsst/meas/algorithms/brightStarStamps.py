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
from typing import NamedTuple
from enum import Enum, auto

import lsst.afw.image as afwImage
import lsst.afw.fits as afwFits
from lsst.afw.geom import TransformPoint2ToPoint2
from lsst.geom import Box2I, Point2I, Extent2I
from lsst.daf.base import PropertySet


class RadiiEnum(Enum):
    INNER_RADIUS = auto()
    OUTER_RADIUS = auto()

    def __str__(self):
        return self.name


class BrightStarStamp(NamedTuple):
    """Single stamp centered on a bright star, normalized by its
    annularFlux.
    """
    starStamp: afwImage.maskedImage.MaskedImageF
    gaiaGMag: float
    gaiaId: int
    annularFlux: float
    XY0: Point2I
    transform: TransformPoint2ToPoint2


class BrightStarStamps(collections.abc.Sequence):
    """Collection of bright star stamps and associated metadata.

    Parameters
    ----------
    starStamps : `collections.abc.Sequence` [`BrightStarStamp`]
        Sequence of star stamps.
    innerRadius : `int`, optional
        Inner radius value, in pixels. This and ``outerRadius`` define the
        annulus used to compute the ``"annularFlux"`` values within each
        ``starStamp``. Must be provided if ``"INNER_RADIUS"`` and
        ``"OUTER_RADIUS"`` are not present in ``metadata``.
    outerRadius : `int`, optional
        Outer radius value, in pixels. This and ``innerRadius`` define the
        annulus used to compute the ``"annularFlux"`` values within each
        ``starStamp``. Must be provided if ``"INNER_RADIUS"`` and
        ``"OUTER_RADIUS"`` are not present in ``metadata``.
    nb90Rots : `int`, optional
        Number of 90 degree rotations required to compensate for detector
        orientation.
    metadata : `lsst.daf.base.PropertyList`, optional
        Metadata associated with the bright stars.

    Raises
    ------
    ValueError
        Raised if one of the star stamps provided does not contain the
        required keys.
    AttributeError
        Raised if the definition of the annulus used to compute each star's
        normalization factor are not provided, that is, if ``"INNER_RADIUS"``
        and ``"OUTER_RADIUS"`` are not present in ``metadata`` _and_
        ``innerRadius`` and ``outerRadius`` are not provided.

    Notes
    -----
    A (gen2) butler can be used to read only a part of the stamps,
    specified by a bbox:

    >>> starSubregions = butler.get("brightStarStamps_sub", dataId, bbox=bbox)
    """

    def __init__(self, starStamps, innerRadius=None, outerRadius=None,
                 nb90Rots=None, metadata=None,):
        for item in starStamps:
            if not isinstance(item, BrightStarStamp):
                raise ValueError(f"Can only add instances of BrightStarStamp, got {type(item)}")
        self._starStamps = starStamps
        self.nb90Rots = nb90Rots
        self._metadata = PropertySet() if metadata is None else metadata.deepCopy()
        # Add inner and outer radii to metadata
        self._checkRadius(innerRadius, RadiiEnum.INNER_RADIUS)
        self._innerRadius = innerRadius
        self._checkRadius(outerRadius, RadiiEnum.OUTER_RADIUS)
        self._outerRadius = outerRadius

    def __len__(self):
        return len(self._starStamps)

    def __getitem__(self, index):
        return self._starStamps[index]

    def __iter__(self):
        return iter(self._starStamps)

    def append(self, item, innerRadius, outerRadius):
        """Add an additional bright star stamp.

        Parameters
        ----------
        item : `BrightStarStamp`
            Bright star stamp to append.
        innerRadius : `int`
            Inner radius value, in4 pixels. This and ``outerRadius`` define the
            annulus used to compute the ``"annularFlux"`` values within each
            ``starStamp``.
        outerRadius : `int`, optional
            Outer radius value, in pixels. This and ``innerRadius`` define the
            annulus used to compute the ``"annularFlux"`` values within each
            ``starStamp``.
        """
        if not isinstance(item, BrightStarStamp):
            raise ValueError(f"Can only add instances of BrightStarStamp, got {type(item)}.")
        self._checkRadius(innerRadius, RadiiEnum.INNER_RADIUS)
        self._checkRadius(outerRadius, RadiiEnum.OUTER_RADIUS)
        self._starStamps.append(item)
        return None

    def extend(self, bss):
        """Extend BrightStarStamps instance by appending elements from another
        instance.

        Parameters
        ----------
        bss : `BrightStarStamps`
            Other instance to concatenate.
        """
        self._checkRadius(bss._innerRadius, RadiiEnum.INNER_RADIUS)
        self._checkRadius(bss._outerRadius, RadiiEnum.OUTER_RADIUS)
        self._starStamps += bss._starStamps

    def getMaskedImages(self):
        """Retrieve star images.

        Returns
        -------
        maskedImages :
            `list` [`lsst.afw.image.maskedImage.maskedImage.MaskedImageF`]
        """
        return [stamp.starStamp for stamp in self._starStamps]

    def getMagnitudes(self):
        """Retrieve Gaia G magnitudes for each star.

        Returns
        -------
        gaiaGMags : `list` [`float`]
        """
        return [stamp.gaiaGMag for stamp in self._starStamps]

    def getGaiaIds(self):
        """Retrieve Gaia IDs for each star.

        Returns
        -------
        gaiaIds : `list` [`int`]
        """
        return [stamp.gaiaId for stamp in self._starStamps]

    def getAnnularFluxes(self):
        """Retrieve normalization factors for each star.

        These are computed by integrating the flux in annulus centered on the
        bright star, far enough from center to be beyond most severe ghosts and
        saturation. The inner and outer radii that define the annulus can be
        recovered from the metadata.

        Returns
        -------
        annularFluxes : list[`float`]
        """
        return [stamp.annularFlux for stamp in self._starStamps]

    def getXY0s(self):
        """Retrieve the coordinates of the bottom-left pixel for each star.

        These correspond to that quantity before warping and rotations are
        applied.

        Returns
        -------
        XY0s : list[`tuple`]
        """
        return [stamp.XY0 for stamp in self._starStamps]

    def getTransforms(self):
        """Retrieve Transform from each star's initial stamp to the common
        model grid.

        Returns
        -------
        transforms : `list` [`TransformPoint2toPoint2`]
        """
        return [stamp.transform for stamp in self._starStamps]

    def selectByMag(self, magMin=None, magMax=None):
        """Return the subset of bright star stamps for objects with specified
        magnitude cuts (in Gaia G).

        Parameters
        ----------
        magMin : `float`, optional
            Keep only stars fainter than this value.
        magMax : `float`, optional
            Keep only stars brighter than this value.
        """
        subset = [stamp for stamp in self._starStamps
                  if (magMin is None or stamp.gaiaGMag > magMin)
                  and (magMax is None or stamp.gaiaGMag < magMax)]
        # This is an optimization to save looping over the init argument when
        # it is already guaranteed to be the correct type
        instance = BrightStarStamps((), metadata=self._metadata)
        instance._starStamps = subset
        return instance

    @property
    def metadata(self):
        return self._metadata.deepCopy()

    def _checkRadius(self, radiusValue, metadataEnum):
        """Ensure provided annulus radius is consistent with that present
        in metadata. If metadata does not contain annulus radius, add it.
        """
        # if a radius value is already present in metadata, ensure it matches
        # the one given
        metadataName = str(metadataEnum)
        if self._metadata.exists(metadataName):
            if radiusValue is not None:
                if self._metadata[metadataName] != radiusValue:
                    raise AttributeError("BrightStarStamps instance already contains different annulus radii "
                                         + f"values ({metadataName}).")
        # if not already in metadata, a value must be provided
        elif radiusValue is None:
            raise AttributeError("No radius value provided for the AnnularFlux measurement "
                                 + f"({metadataName}), and none present in metadata.")
        else:
            self._metadata[metadataName] = radiusValue
        return None

    def writeFits(self, filename):
        """Write a single FITS file containing all bright star stamps.
        """
        # ensure metadata contains current number of objects
        self._metadata["N_STARS"] = len(self)

        # if class instance contains number of rotations, save it to header
        if self.nb90Rots is not None:
            self._metadata["NB_90_ROTS"] = self.nb90Rots

        # add full list of Gaia magnitudes, IDs, annularFluxes and transforms
        # to shared metadata
        self._metadata["G_MAGS"] = self.getMagnitudes()
        self._metadata["GAIA_IDS"] = self.getGaiaIds()
        self._metadata["X0S"] = [XY0[0] for XY0 in self.getXY0s()]
        self._metadata["Y0S"] = [XY0[1] for XY0 in self.getXY0s()]
        self._metadata["ANNULAR_FLUXES"] = self.getAnnularFluxes()
        # create primary HDU with global metadata
        fitsPrimary = afwFits.Fits(filename, "w")
        fitsPrimary.createEmpty()
        fitsPrimary.writeMetadata(self._metadata)
        fitsPrimary.closeFile()

        # add all stamps and mask planes
        for stamp, transform in zip(self.getMaskedImages(), self.getTransforms()):
            stamp.getImage().writeFits(filename, mode='a')
            stamp.getMask().writeFits(filename, mode='a')
            transform.writeFits(filename, mode='a')
        return None

    @classmethod
    def readFits(cls, filename):
        """Read bright star stamps from FITS file.

        Returns
        -------
        bss : `BrightStarStamps`
            Collection of bright star stamps.
        """
        bss = cls.readFitsWithOptions(filename, None)
        return bss

    @classmethod
    def readFitsWithOptions(cls, filename, options):
        """Read bright star stamps from FITS file, allowing for only a
        subregion of the stamps to be read.

        Returns
        -------
        bss : `BrightStarStamps`
            Collection of bright star stamps.
        """
        # extract necessary info from metadata
        visitMetadata = afwFits.readMetadata(filename, hdu=0)
        nbStarStamps = visitMetadata["N_STARS"]
        gaiaGMags = visitMetadata.getArray("G_MAGS")
        gaiaIds = visitMetadata.getArray("GAIA_IDS")
        XY0s = [Point2I(x0, y0) for x0, y0 in zip(visitMetadata.getArray("X0S"),
                                                  visitMetadata.getArray("Y0S"))]
        annularFluxes = visitMetadata.getArray("ANNULAR_FLUXES")
        try:
            nb90Rots = visitMetadata["NB_90_ROTS"]
        except KeyError:
            nb90Rots = None
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
        for bStarIdx in range(nbStarStamps):
            # read and reassign HDUs. Note TransformPoint2ToPoint2 writes two
            # HDUs per call of the write method, hence the 4 here.
            imReader = afwImage.ImageFitsReader(filename, hdu=4*bStarIdx + 1)
            maskReader = afwImage.MaskFitsReader(filename, hdu=4*bStarIdx + 2)
            maskedImage = afwImage.MaskedImageF(image=imReader.read(**kwargs),
                                                mask=maskReader.read(**kwargs))
            transform = TransformPoint2ToPoint2.readFits(filename, hdu=4*bStarIdx + 3)
            starStamps.append(BrightStarStamp(starStamp=maskedImage,
                                              gaiaGMag=gaiaGMags[bStarIdx],
                                              gaiaId=gaiaIds[bStarIdx],
                                              XY0=XY0s[bStarIdx],
                                              annularFlux=annularFluxes[bStarIdx],
                                              transform=transform))
        bss = cls(starStamps, nb90Rots=nb90Rots, metadata=visitMetadata)
        return bss
