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

__all__ = ["MatcherSourceSelectorConfig", "MatcherSourceSelectorTask"]

from deprecated.sphinx import deprecated

import numpy as np

import lsst.pex.config as pexConfig
from .sourceSelector import BaseSourceSelectorConfig, BaseSourceSelectorTask, sourceSelectorRegistry
from lsst.pipe.base import Struct


class MatcherSourceSelectorConfig(BaseSourceSelectorConfig):
    sourceFluxType = pexConfig.Field(
        doc="Type of source flux; typically one of Ap or Psf",
        dtype=str,
        default="Ap",
    )
    minSnr = pexConfig.Field(
        dtype=float,
        doc="Minimum allowed signal-to-noise ratio for sources used for matching "
        "(in the flux specified by sourceFluxType); <= 0 for no limit",
        default=40,
    )
    excludePixelFlags = pexConfig.Field(
        dtype=bool,
        doc="Exclude objects that have saturated, interpolated, or edge "
            "pixels using PixelFlags. For matchOptimisticB set this to False "
            "to recover previous matcher selector behavior.",
        default=True,
    )


# remove this file on DM-41146
@deprecated(reason=("This Task has been replaced by an appropriately configured ScienceSourceSelector."
                    " Will be removed after v27."),
            version="v27.0", category=FutureWarning)
@pexConfig.registerConfigurable("matcher", sourceSelectorRegistry)
class MatcherSourceSelectorTask(BaseSourceSelectorTask):
    """Select sources that are useful for matching.

    Good matching sources have high signal/noise, are non-blended. They need not
    be PSF sources, just have reliable centroids.

    Distinguished from astrometrySourceSelector because it is more lenient
    (i.e. not checking footprints or bad flags).
    """
    ConfigClass = MatcherSourceSelectorConfig

    def __init__(self, *args, **kwargs):
        BaseSourceSelectorTask.__init__(self, *args, **kwargs)

    def selectSources(self, sourceCat, matches=None, exposure=None):
        """Return a selection of sources that are useful for matching.

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

        good = self._isUsable(sourceCat)
        if self.config.excludePixelFlags:
            good = good & self._isGood(sourceCat)
        return Struct(selected=good)

    def _getSchemaKeys(self, schema):
        """Extract and save the necessary keys from schema with asKey.
        """
        self.parentKey = schema["parent"].asKey()
        self.centroidXKey = schema["slot_Centroid_x"].asKey()
        self.centroidYKey = schema["slot_Centroid_y"].asKey()
        self.centroidFlagKey = schema["slot_Centroid_flag"].asKey()

        fluxPrefix = "slot_%sFlux_" % (self.config.sourceFluxType,)
        self.fluxField = fluxPrefix + "instFlux"
        self.fluxKey = schema[fluxPrefix + "instFlux"].asKey()
        self.fluxFlagKey = schema[fluxPrefix + "flag"].asKey()
        self.fluxErrKey = schema[fluxPrefix + "instFluxErr"].asKey()

        self.edgeKey = schema["base_PixelFlags_flag_edge"].asKey()
        self.interpolatedCenterKey = schema["base_PixelFlags_flag_interpolatedCenter"].asKey()
        self.saturatedKey = schema["base_PixelFlags_flag_saturated"].asKey()

    def _isParent(self, sourceCat):
        """Return True for each source that is the parent source.
        """
        test = (sourceCat[self.parentKey] == 0)
        return test

    def _hasCentroid(self, sourceCat):
        """Return True for each source that has a valid centroid
        """
        return np.isfinite(sourceCat[self.centroidXKey]) \
            & np.isfinite(sourceCat[self.centroidYKey]) \
            & ~sourceCat[self.centroidFlagKey]

    def _goodSN(self, sourceCat):
        """Return True for each source that has Signal/Noise > config.minSnr.
        """
        if self.config.minSnr <= 0:
            return True
        else:
            with np.errstate(invalid="ignore"):  # suppress NAN warnings
                return sourceCat[self.fluxKey]/sourceCat[self.fluxErrKey] > self.config.minSnr

    def _isUsable(self, sourceCat):
        """
        Return True for each source that is usable for matching, even if it may
        have a poor centroid.

        For a source to be usable it must:

        - have a valid centroid
        - not be deblended
        - have a valid instFlux (of the type specified in this object's constructor)
        - have adequate signal-to-noise
        """
        return self._hasCentroid(sourceCat) \
            & self._isParent(sourceCat) \
            & self._goodSN(sourceCat) \
            & ~sourceCat[self.fluxFlagKey]

    def _isGood(self, sourceCat):
        """
        Return True for each source that is usable for matching, even if it may
        have a poor centroid.

        For a source to be usable it must:

        - Not be on a CCD edge.
        - Not have an interpolated pixel within 3x3 around their centroid.
        - Not have a saturated pixel in their footprint.
        """
        return ~sourceCat[self.edgeKey] & \
               ~sourceCat[self.interpolatedCenterKey] & \
               ~sourceCat[self.saturatedKey]
