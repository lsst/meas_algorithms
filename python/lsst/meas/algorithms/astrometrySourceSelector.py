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
from __future__ import absolute_import, division, print_function

import numpy as np

from lsst.afw import table
import lsst.pex.config as pexConfig
from .sourceSelector import BaseSourceSelectorConfig, BaseSourceSelectorTask, sourceSelectorRegistry
from lsst.pipe.base import Struct
from functools import reduce


class AstrometrySourceSelectorConfig(BaseSourceSelectorConfig):
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


class AstrometrySourceSelectorTask(BaseSourceSelectorTask):
    """
    !Select sources that are useful for astrometry.

    Good astrometry sources have high signal/noise, are non-blended, and
    did not have certain "bad" flags set during source extraction. They need not
    be PSF sources, just have reliable centroids.
    """
    ConfigClass = AstrometrySourceSelectorConfig

    def __init__(self, *args, **kwargs):
        BaseSourceSelectorTask.__init__(self, *args, **kwargs)

    def selectSources(self, sourceCat, matches=None):
        """
        !Return a catalog of sources: a subset of sourceCat.

        If sourceCat is cotiguous in memory, will use vectorized tests for ~100x
        execution speed advantage over non-contiguous catalogs. This would be
        even faster if we didn't have to check footprints for multiple peaks.

        @param[in] sourceCat  catalog of sources that may be sources
                                (an lsst.afw.table.SourceCatalog)

        @return a pipeBase.Struct containing:
        - sourceCat  a catalog of sources
        """
        self._getSchemaKeys(sourceCat.schema)

        if sourceCat.isContiguous():
            bad = reduce(lambda x, y: np.logical_or(x, sourceCat.get(y)), self.config.badFlags, False)
            good = self._isGood_vector(sourceCat)
            result = sourceCat[good & ~bad]
        else:
            result = table.SourceCatalog(sourceCat.table)
            for i, source in enumerate(sourceCat):
                if self._isGood(source) and not self._isBad(source):
                    result.append(source)
        return Struct(sourceCat=result)

    def _getSchemaKeys(self, schema):
        """Extract and save the necessary keys from schema with asKey."""
        self.parentKey = schema["parent"].asKey()
        self.nChildKey = schema["deblend_nChild"].asKey()
        self.centroidXKey = schema["slot_Centroid_x"].asKey()
        self.centroidYKey = schema["slot_Centroid_y"].asKey()
        self.centroidFlagKey = schema["slot_Centroid_flag"].asKey()

        self.edgeKey = schema["base_PixelFlags_flag_edge"].asKey()
        self.interpolatedCenterKey = schema["base_PixelFlags_flag_interpolatedCenter"].asKey()
        self.saturatedKey = schema["base_PixelFlags_flag_saturated"].asKey()

        fluxPrefix = "slot_%sFlux_" % (self.config.sourceFluxType,)
        self.fluxKey = schema[fluxPrefix + "flux"].asKey()
        self.fluxFlagKey = schema[fluxPrefix + "flag"].asKey()
        self.fluxSigmaKey = schema[fluxPrefix + "fluxSigma"].asKey()

    def _isMultiple_vector(self, sourceCat):
        """Return True for each source that is likely multiple sources."""
        test = (sourceCat.get(self.parentKey) != 0) | (sourceCat.get(self.nChildKey) != 0)
        # have to currently manage footprints on a source-by-source basis.
        for i, cat in enumerate(sourceCat):
            footprint = cat.getFootprint()
            test[i] |= (footprint is not None) and (len(footprint.getPeaks()) > 1)
        return test

    def _isMultiple(self, source):
        """Return True if source is likely multiple sources."""
        if (source.get(self.parentKey) != 0) or (source.get(self.nChildKey) != 0):
            return True
        footprint = source.getFootprint()
        return footprint is not None and len(footprint.getPeaks()) > 1

    def _hasCentroid_vector(self, sourceCat):
        """Return True for each source that has a valid centroid"""
        return np.isfinite(sourceCat.get(self.centroidXKey)) \
            & np.isfinite(sourceCat.get(self.centroidYKey)) \
            & ~sourceCat.get(self.centroidFlagKey)

    def _hasCentroid(self, source):
        """Return True if the source has a valid centroid"""
        centroid = source.getCentroid()
        return np.all(np.isfinite(centroid)) and not source.getCentroidFlag()

    def _goodSN_vector(self, sourceCat):
        """Return True for each source that has Signal/Noise > config.minSnr."""
        if self.config.minSnr <= 0:
            return True
        else:
            return sourceCat.get(self.fluxKey)/sourceCat.get(self.fluxSigmaKey) > self.config.minSnr

    def _goodSN(self, source):
        """Return True if source has Signal/Noise > config.minSnr."""
        return (self.config.minSnr <= 0 or
                (source.get(self.fluxKey)/source.get(self.fluxSigmaKey) > self.config.minSnr))

    def _isUsable_vector(self, sourceCat):
        """
        Return True for each source that is usable for matching, even if it may
        have a poor centroid.

        For a source to be usable it must:
        - have a valid centroid
        - not be deblended
        - have a valid flux (of the type specified in this object's constructor)
        - have adequate signal-to-noise
        """

        return self._hasCentroid_vector(sourceCat) \
            & ~self._isMultiple_vector(sourceCat) \
            & self._goodSN_vector(sourceCat) \
            & ~sourceCat.get(self.fluxFlagKey)

    def _isUsable(self, source):
        """
        Return True if the source is usable for matching, even if it may have a
        poor centroid.

        For a source to be usable it must:
        - have a valid centroid
        - not be deblended
        - have a valid flux (of the type specified in this object's constructor)
        - have adequate signal-to-noise
        """
        return self._hasCentroid(source) \
            and not self._isMultiple(source) \
            and not source.get(self.fluxFlagKey) \
            and self._goodSN(source)

    def _isGood_vector(self, sourceCat):
        """
        Return True for each source that is usable for matching and likely has a
        good centroid.

        The additional tests for a good centroid, beyond isUsable, are:
        - not interpolated in the center
        - not saturated
        - not near the edge
        """

        return self._isUsable_vector(sourceCat) \
            & ~sourceCat.get(self.saturatedKey) \
            & ~sourceCat.get(self.interpolatedCenterKey) \
            & ~sourceCat.get(self.edgeKey)

    def _isGood(self, source):
        """
        Return True if source is usable for matching and likely has a good centroid.

        The additional tests for a good centroid, beyond isUsable, are:
        - not interpolated in the center
        - not saturated
        - not near the edge
        """
        return self._isUsable(source) \
            and not source.get(self.saturatedKey) \
            and not source.get(self.interpolatedCenterKey) \
            and not source.get(self.edgeKey)


sourceSelectorRegistry.register("astrometry", AstrometrySourceSelectorTask)
