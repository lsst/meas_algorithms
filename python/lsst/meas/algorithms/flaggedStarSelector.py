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

from .starSelector import BaseStarSelectorTask, starSelectorRegistry

__all__ = ["FlaggedStarSelectorConfig", "FlaggedStarSelectorTask"]


class FlaggedStarSelectorConfig(BaseStarSelectorTask.ConfigClass):
    field = lsst.pex.config.Field(
        dtype=str, default="calib_psfUsed",
        doc="Name of a flag field that is True for stars that should be "
            "used.",
    )


class FlaggedStarSelectorTask(BaseStarSelectorTask):
    """A trivial StarSelector that simply uses an existing flag field to filter
    a SourceCatalog.

    This is most frequently used in steps that occur after the a PSF model has
    been built, to allow other procedures that need Stars to use the set of
    Stars used to determine the PSF.

    Attributes
    ----------
    usesMatches : bool
        A boolean variable specify if the inherited star selector uses matches.
    key : lsst.afw.table.Key
        Schema key specifying which catalog column flag to select on.
    """

    ConfigClass = FlaggedStarSelectorConfig
    _DefaultName = "flagged"

    def __init__(self, schema, **kwds):
        BaseStarSelectorTask.__init__(self, schema, **kwds)
        self.key = schema.find(self.config.field).key

    def selectStars(self, exposure, sourceCat, maskedImage=None, matches=None):
        """Return a bool array representing which stars to select from
        sourceCat.

        The input catalog must be contiguous in memory.

        Parameters
        ----------
        sourceCat : lsst.afw.table.SourceCatalog
            Catalog of stars to select from.
        maskedImage : {None} lsst.afw.image
            An image containing the stars for use in selection tests or for
            plotting.
        matches : {None} list of lsst.afw.table.ReferenceMatch
            A list of lsst.afw.table.ReferenceMatch objects

        Return
        ------
        lsst.pipe.base.Struct
            The struct contains the following data:

            starCat : bool array
                Boolean array of stars that were selected, same length as
                sourceCat.
        """
        return pipeBase.Struct(
            starCat=sourceCat[sourceCat.get(self.key)])

starSelectorRegistry.register("flagged", FlaggedStarSelectorTask)
