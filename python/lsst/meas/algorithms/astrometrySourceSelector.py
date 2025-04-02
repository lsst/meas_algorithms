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

"""Select sources that are useful for astrometry.

Such sources have good signal-to-noise, are well centroided, not blended,
and not flagged with a handful of "bad" flags.
"""

__all__ = ["AstrometrySourceSelectorConfig", "AstrometrySourceSelectorTask"]

from deprecated.sphinx import deprecated

import numpy as np

import lsst.pex.config as pexConfig
from lsst.pex.exceptions import NotFoundError
from .sourceSelector import BaseSourceSelectorConfig, BaseSourceSelectorTask, sourceSelectorRegistry
from lsst.pipe.base import Struct
from functools import reduce


class AstrometrySourceSelectorConfig(BaseSourceSelectorConfig):
    badFlags = pexConfig.ListField(
        doc="List of flags which cause a source to be rejected as bad",
        dtype=str,
        default=[
            "base_PixelFlags_flag_edge",
            "base_PixelFlags_flag_nodata",
            "base_PixelFlags_flag_interpolatedCenter",
            "base_PixelFlags_flag_saturatedCenter",
            "base_PixelFlags_flag_crCenter",
            "base_PixelFlags_flag_bad",
        ],
    )
    sourceFluxType = pexConfig.Field(
        doc="Type of source flux; typically one of Ap or Psf",
        dtype=str,
        default="Ap",
    )
    minSnr = pexConfig.Field(
        dtype=float,
        doc="Minimum allowed signal-to-noise ratio for sources used for matching "
        "(in the flux specified by sourceFluxType); <= 0 for no limit",
        default=10,
    )


# remove this file on DM-41146
@deprecated(reason=("This Task has been replaced by an appropriately configured ScienceSourceSelector."
                    " See `AstrometryConfig.setDefaults` in meas_astrom for an example config that was "
                    "made to closely match this Task. Will be removed after v27."),
            version="v27.0", category=FutureWarning)
@pexConfig.registerConfigurable("astrometry", sourceSelectorRegistry)
class AstrometrySourceSelectorTask(BaseSourceSelectorTask):
    """Select sources that are useful for astrometry.

    Good astrometry sources have high signal/noise, are non-blended, and
    did not have certain "bad" flags set during source extraction. They need not
    be PSF sources, just have reliable centroids.
    """
    ConfigClass = AstrometrySourceSelectorConfig

    def __init__(self, *args, **kwargs):
        BaseSourceSelectorTask.__init__(self, *args, **kwargs)

    def selectSources(self, sourceCat, matches=None, exposure=None):
        """Return a selection of sources that are useful for astrometry.

        Parameters
        ----------
        sourceCat : `lsst.afw.table.SourceCatalog`
            Catalog of sources to select from.
            This catalog must be contiguous in memory.
        matches : `list` of `lsst.afw.table.ReferenceMatch` or None
            Ignored in this SourceSelector.
        exposure : `lsst.afw.image.Exposure` or None
            The exposure the catalog was built from; used for debug display.

        Returns
        -------
        struct : `lsst.pipe.base.Struct`
            The struct contains the following data:

            ``selected``
                Boolean array of sources that were selected, same length as
                sourceCat. (`numpy.ndarray` of `bool`)
        """
        self._getSchemaKeys(sourceCat.schema)

        bad = reduce(lambda x, y: np.logical_or(x, sourceCat[y]), self.config.badFlags, False)
        good = self._isGood(sourceCat)
        return Struct(selected=good & ~bad)

    def _getSchemaKeys(self, schema):
        """Extract and save the necessary keys from schema with asKey.
        """
        self.parentKey = schema["parent"].asKey()
        self.nChildKey = schema["deblend_nChild"].asKey()
        self.centroidXKey = schema["slot_Centroid_x"].asKey()
        self.centroidYKey = schema["slot_Centroid_y"].asKey()
        self.centroidXErrKey = schema["slot_Centroid_xErr"].asKey()
        self.centroidYErrKey = schema["slot_Centroid_yErr"].asKey()
        self.centroidFlagKey = schema["slot_Centroid_flag"].asKey()
        try:
            self.primaryKey = schema["detect_isPrimary"].asKey()
        except NotFoundError:
            self.primaryKey = None

        self.edgeKey = schema["base_PixelFlags_flag_edge"].asKey()
        self.interpolatedCenterKey = schema["base_PixelFlags_flag_interpolatedCenter"].asKey()
        self.saturatedKey = schema["base_PixelFlags_flag_saturated"].asKey()

        fluxPrefix = "slot_%sFlux_" % (self.config.sourceFluxType,)
        self.instFluxKey = schema[fluxPrefix + "instFlux"].asKey()
        self.fluxFlagKey = schema[fluxPrefix + "flag"].asKey()
        self.instFluxErrKey = schema[fluxPrefix + "instFluxErr"].asKey()

    def _isMultiple(self, sourceCat):
        """Return True for each source that is likely multiple sources.
        """
        test = (sourceCat[self.parentKey] != 0) | (sourceCat[self.nChildKey] != 0)
        # have to currently manage footprints on a source-by-source basis.
        for i, cat in enumerate(sourceCat):
            footprint = cat.getFootprint()
            test[i] |= (footprint is not None) and (len(footprint.getPeaks()) > 1)
        return test

    def _hasCentroid(self, sourceCat):
        """Return True for each source that has a valid centroid
        """
        def checkNonfiniteCentroid():
            """Return True for sources with non-finite centroids.
            """
            return ~np.isfinite(sourceCat[self.centroidXKey]) | \
                ~np.isfinite(sourceCat[self.centroidYKey])
        assert ~checkNonfiniteCentroid().any(), \
            "Centroids not finite for %d unflagged sources." % (checkNonfiniteCentroid().sum())
        return np.isfinite(sourceCat[self.centroidXErrKey]) \
            & np.isfinite(sourceCat[self.centroidYErrKey]) \
            & ~sourceCat[self.centroidFlagKey]

    def _goodSN(self, sourceCat):
        """Return True for each source that has Signal/Noise > config.minSnr.
        """
        if self.config.minSnr <= 0:
            return True
        else:
            with np.errstate(invalid="ignore"):  # suppress NAN warnings
                return sourceCat[self.instFluxKey]/sourceCat[self.instFluxErrKey] > self.config.minSnr

    def _isUsable(self, sourceCat):
        """Return True for each source that is usable for matching, even if it may
        have a poor centroid.

        For a source to be usable it must:
        - have a valid centroid
        - not be deblended
        - have a valid flux (of the type specified in this object's constructor)
        - have adequate signal-to-noise
        """

        return self._hasCentroid(sourceCat) \
            & ~self._isMultiple(sourceCat) \
            & self._goodSN(sourceCat) \
            & ~sourceCat[self.fluxFlagKey]

    def _isPrimary(self, sourceCat):
        """Return True if this is a primary source.
        """
        if self.primaryKey:
            return sourceCat[self.primaryKey]
        else:
            return np.ones(len(sourceCat), dtype=bool)

    def _isGood(self, sourceCat):
        """Return True for each source that is usable for matching and likely has a
        good centroid.

        The additional tests for a good centroid, beyond isUsable, are:
        - not interpolated in the center
        - not saturated
        - not near the edge
        """

        return self._isUsable(sourceCat) \
            & self._isPrimary(sourceCat) \
            & ~sourceCat[self.saturatedKey] \
            & ~sourceCat[self.interpolatedCenterKey] \
            & ~sourceCat[self.edgeKey]

    def _isBadFlagged(self, source):
        """Return True if any of config.badFlags are set for this source.
        """
        return any(source[flag] for flag in self.config.badFlags)
