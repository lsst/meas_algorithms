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

import lsst.pex.config
import lsst.afw.table
import lsst.pipe.base as pipeBase

from .sourceSelector import BaseSourceSelectorTask, sourceSelectorRegistry

__all__ = ["FlaggedSourceSelectorConfig", "FlaggedSourceSelectorTask"]


class FlaggedSourceSelectorConfig(BaseSourceSelectorTask.ConfigClass):
    field = lsst.pex.config.Field(
        dtype=str, default="calib_psfUsed",
        doc="Name of a flag field that is True for Sources that should be "
            "used.",
    )


class FlaggedSourceSelectorTask(BaseSourceSelectorTask):
    """!
    A trivial SourceSelector that simply uses an existing flag field to filter
    a SourceCatalog.

    This is most frequently used in steps that occur after the a PSF model has
    been built, to allow other procedures that need Sources to use the set of
    Sources used to determine the PSF.

    Attributes
    ----------
    usesMatches : bool
        A boolean variable specify if the inherited source selector uses
        matches.
    key : lsst.afw.table.Key
        Schema key specifying which catalog column flag to select on.
    """

    ConfigClass = FlaggedSourceSelectorConfig
    _DefaultName = "flagged"

    def __init__(self, schema, **kwds):
        BaseSourceSelectorTask.__init__(self, **kwds)
        self.key = schema.find(self.config.field).key

    def selectSources(self, sourceCat, maskedImage=None, matches=None):
        """Return a bool array representing which sources to select from
        sourceCat.

        The input catalog must be contiguous in memory.

        Parameters
        ----------
        sourceCat : lsst.afw.table.SourceCatalog
            Catalog of sources to select from.
        maskedImage : {None} lsst.afw.image
            An image containing the sources for use in selection tests or for
            plotting.
        matches : {None} list of lsst.afw.table.ReferenceMatch
            A list of lsst.afw.table.ReferenceMatch objects

        Return
        ------
        lsst.pipe.base.Struct
            The struct contains the following data:

            selected : bool array
                Boolean array of sources that were selected, same length as
                sourceCat.
        """
        return pipeBase.Struct(
            selected=sourceCat.get(self.key),)

sourceSelectorRegistry.register("flagged", FlaggedSourceSelectorTask)
