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

__all__ = ["BaseSourceSelectorConfig", "BaseSourceSelectorTask", "sourceSelectorRegistry"]

import abc

import lsst.afw.table as afwTable
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from future.utils import with_metaclass


class BaseSourceSelectorConfig(pexConfig.Config):
    badFlags = pexConfig.ListField(
        doc="List of flags which cause a source to be rejected as bad",
        dtype=str,
        default=[
            "base_PixelFlags_flag_edge",
            "base_PixelFlags_flag_interpolatedCenter",
            "base_PixelFlags_flag_saturatedCenter",
            "base_PixelFlags_flag_crCenter",
            "base_PixelFlags_flag_bad",
            "base_PixelFlags_flag_interpolated",
        ],
    )


class BaseSourceSelectorTask(with_metaclass(abc.ABCMeta, pipeBase.Task)):
    """!Base class for source selectors

    Register all source selectors with the sourceSelectorRegistry using:
        sourceSelectorRegistry.register(name, class)
    """

    ConfigClass = BaseSourceSelectorConfig
    _DefaultName = "sourceSelector"

    def __init__(self, **kwargs):
        """!Initialize a source selector."""
        pipeBase.Task.__init__(self, **kwargs)

    def run(self, sourceCat, maskedImage=None, **kwargs):
        """!Select sources and return them.

        @param[in] sourceCat  catalog of sources that may be sources (an lsst.afw.table.SourceCatalog)
        @param[in] maskedImage  the maskedImage containing the sources, for plotting.

        @return an lsst.pipe.base.Struct containing:
        - sourceCat  catalog of sources that were selected
        """
        return self.selectSources(maskedImage=maskedImage, sourceCat=sourceCat, **kwargs)

    @abc.abstractmethod
    def selectSources(self, sourceCat, matches=None):
        """!Return a catalog of sources: a subset of sourceCat.

        @param[in] sourceCat  catalog of sources that may be sources (an lsst.afw.table.SourceCatalog)

        @return a pipeBase.Struct containing:
        - sourceCat  a catalog of sources
        """

        # NOTE: example implementation, returning all sources that have no bad flags set.
        result = afwTable.SourceCatalog(sourceCat.table)
        for source in sourceCat:
            if not self._isBad(source):
                result.append(source)
        return pipeBase.Struct(sourceCat=result)

    def _isBad(self, source):
        """Return True if any of config.badFlags are set for this source."""
        return any(source.get(flag) for flag in self.config.badFlags)


sourceSelectorRegistry = pexConfig.makeRegistry(
    doc="A registry of source selectors (subclasses of BaseSourceSelectorTask)",
)
