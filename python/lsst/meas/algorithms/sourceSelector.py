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

__all__ = ["BaseSourceSelectorConfig", "BaseSourceSelectorTask",
           "sourceSelectorRegistry"]

import abc

import lsst.afw.table as afwTable
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from future.utils import with_metaclass


class BaseSourceSelectorConfig(pexConfig.Config):
    bad_flags = pexConfig.ListField(
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


@with_metaclass(abc.ABCMeta)
class BaseSourceSelectorTask(pipeBase.Task):
    """Base class for source selectors

    Write paragraph about what a source selector is/does.

    Register all source selectors with the sourceSelectorRegistry using:
        sourceSelectorRegistry.register(name, class)

    Attributes
    ----------
    usesMatches : bool
        A boolean variable specifiy if the inherited source selector uses
        matches.
    """

    ConfigClass = BaseSourceSelectorConfig
    _DefaultName = "source_selector"
    uses_matches = False

    def __init__(self, **kwargs):
        """Create a source selector.

        Source selectors may have their own special arguments plus
        any kwargs for lsst.pipe.base.Task
        """
        pipeBase.Task.__init__(self, **kwargs)

    def run(self, source_cat, source_selected_field=None, masked_image=None,
            matches=None):
        """Select sources and return them.

        Parameters:
        -----------
        source_cat : lsst.afw.table.SourceCatalog
            Catalog of sources that may be sources
        source_selected_field : {None} lsst.afw.table.Key
            Key specifying a location field to write to. If None then no data
            is written to the output source cat.
        masked_image : {None} lsst.afw.image.MaskedImage
            Masked image containing the sources. A few source selectors may use
            this in their selection but most will just use it for plotting.
        matches : {None} list of lsst.afw.table.ReferenceMatch
            A list of lsst.afw.table.ReferenceMatch objects. If use_matches set
            in source selector, this field is required otherwise ignored.

        Return
        ------
        lsst.pipe.base.Struct
            The struct contains the following data:

            source_cat : lsst.afw.table.SourceCatalog
                The catalog of sources that were selected.
        """
        return self.select_sources(source_cat=source_cat,
                                   source_selected_field=None,
                                   masked_image=masked_image,
                                   matches=matches)

    @abc.abstractmethod
    def select_sources(self, source_cat, source_selected_field=None,
                       masked_image=None, matches=None):
        """Return a catalog of sources selected by specified criteria.

        We would prefer the input catalog to be contiguous in memory
        for speed and simplicity of the code. However, the code does
        allows for the input of non contiguous catalogs. Through
        a different code path.

        Parameters
        ----------
        source_cat : lsst.afw.table.SourceCatalog
            catalog of sources that may be sources
        source_selected_field : {None} lsst.afw.table.Key
            Key specifying a location field to write to
            if the key was set.
        masked_image : {None} lsst.afw.image
            An image containing the sources tests or for plotting.
        matches : {None} list of lsst.afw.table.ReferenceMatch
            A list of lsst.afw.table.ReferenceMatch objects

        Return
        ------
        lsst.pipe.base.Struct
            The struct contains the following data:

            result : lsst.afw.table.SourceCatalog
                The catalog of sources that were selected.
        """

        # NOTE: example implementation, returning all sources that have no bad
        # flags set.
        result = afwTable.SourceCatalog(source_cat.table)

        is_bad_array = self._isBad(source_cat)
        if source_selected_field is not None:
            result.get(source_selected_field) = np.logical_not(
                is_bad_array)
            return pipeBase.Struct(source_cat=result)
        return pipeBase.Struct(source_cat=result[np.logical_not(is_bad_array)])

    def _is_bad(self, source):
        """Return True if any of config.badFlags are set for this source.

        Parameters
        ----------
        source : lsst.afw.tabel.SourceRecord
            Source record object to test again the selection.

        Return
        ------
        bool array
            True if source fails any of the selections.
        """
        return any(source.get(flag) for flag in self.config.bad_flags)



sourceSelectorRegistry = pexConfig.makeRegistry(
    doc="A registry of source selectors (subclasses of BaseSourceSelectorTask)",
)
