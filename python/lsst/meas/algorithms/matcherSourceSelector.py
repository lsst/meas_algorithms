#
# LSST Data Management System
#
# Copyright 2008-2017  AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#

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
        """Return a catalog of sources: a subset of sourceCat.

        The input catalog must be contiguous in memory.

        Parameters:
        -----------
        sourceCat : `lsst.afw.table.SourceCatalog`
            Catalog of sources to select from.
        matches : `list` of `lsst.afw.table.ReferenceMatch` or None
            Ignored in this SourceSelector.
        exposure : `lsst.afw.image.Exposure` or None
            The exposure the catalog was built from; used for debug display.

        Return
        ------
        struct : `lsst.pipe.base.Struct`
            The struct contains the following data:

            - selected : `array` of `bool``
                Boolean array of sources that were selected, same length as
                sourceCat.
        """
        self._getSchemaKeys(sourceCat.schema)

        good = self._isUsable(sourceCat)
        return Struct(selected=good)

    def _getSchemaKeys(self, schema):
        """Extract and save the necessary keys from schema with asKey."""
        self.parentKey = schema["parent"].asKey()
        self.centroidXKey = schema["slot_Centroid_x"].asKey()
        self.centroidYKey = schema["slot_Centroid_y"].asKey()
        self.centroidFlagKey = schema["slot_Centroid_flag"].asKey()

        fluxPrefix = "slot_%sFlux_" % (self.config.sourceFluxType,)
        self.fluxField = fluxPrefix + "flux"
        self.fluxKey = schema[fluxPrefix + "flux"].asKey()
        self.fluxFlagKey = schema[fluxPrefix + "flag"].asKey()
        self.fluxSigmaKey = schema[fluxPrefix + "fluxSigma"].asKey()

    def _isParent(self, sourceCat):
        """Return True for each source that is the parent source."""
        test = (sourceCat.get(self.parentKey) == 0)
        return test

    def _hasCentroid(self, sourceCat):
        """Return True for each source that has a valid centroid"""
        return np.isfinite(sourceCat.get(self.centroidXKey)) \
            & np.isfinite(sourceCat.get(self.centroidYKey)) \
            & ~sourceCat.get(self.centroidFlagKey)

    def _goodSN(self, sourceCat):
        """Return True for each source that has Signal/Noise > config.minSnr."""
        if self.config.minSnr <= 0:
            return True
        else:
            with np.errstate(invalid="ignore"):  # suppress NAN warnings
                return sourceCat.get(self.fluxKey)/sourceCat.get(self.fluxSigmaKey) > self.config.minSnr

    def _isUsable(self, sourceCat):
        """
        Return True for each source that is usable for matching, even if it may
        have a poor centroid.

        For a source to be usable it must:
        - have a valid centroid
        - not be deblended
        - have a valid flux (of the type specified in this object's constructor)
        - have adequate signal-to-noise
        """
        return self._hasCentroid(sourceCat) \
            & self._isParent(sourceCat) \
            & self._goodSN(sourceCat) \
            & ~sourceCat.get(self.fluxFlagKey)


@pexConfig.registerConfigurable("matcherPessimistic", sourceSelectorRegistry)
class MatcherPessimisticSourceSelectorTask(MatcherSourceSelectorTask):
    """Select sources that are useful for matching.

    Good matching sources have high signal/noise, are non-blended. They need not
    be PSF sources, just have reliable centroids. This inherited class adds
    the removal of saturated, interpolated, and edge_key objects to the set of
    bad flags. It is a temporary addition designed preserve the source selction
    used in matchOptimisticB. Once matchPessimisticB is adopted as the default
    source selector the class will be removed and the saturated, interpoalted, and
    edge_key flags will be added to the matcherSourceSelector class.

    TODO: Once DM-10399 is complete an RFC will be filed to make matchPessimisticB
    the default matcher this class will replace matcherSourceSelector with this source
    selector resulting in only one matcherSourceSeletor. The ticket describing
    this work is DM-10800.
    """
    def _getSchemaKeys(self, schema):
        """Extract and save the necessary keys from schema with asKey."""
        MatcherSourceSelectorTask._getSchemaKeys(self, schema)

        self.edgeKey = schema["base_PixelFlags_flag_edge"].asKey()
        self.interpolatedCenterKey = schema["base_PixelFlags_flag_interpolatedCenter"].asKey()
        self.saturatedKey = schema["base_PixelFlags_flag_saturated"].asKey()

    def _isUsable(self, sourceCat):
        """
        Return True for each source that is usable for matching, even if it may
        have a poor centroid.

        For a source to be usable it must:
        - have a valid centroid
        - not be deblended
        - have a valid flux (of the type specified in this object's constructor)
        - have adequate signal-to-noise
        """
        result = MatcherSourceSelectorTask._isUsable(self, sourceCat)

        return result \
            & ~sourceCat.get(self.edgeKey) \
            & ~sourceCat.get(self.interpolatedCenterKey) \
            & ~sourceCat.get(self.saturatedKey)
