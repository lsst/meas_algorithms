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
__all__ = ("FlaggedStarSelectorConfig", "FlaggedStarSelectorTask")

import lsst.pex.config
import lsst.afw.table
import lsst.pipe.base

from .starSelector import BaseStarSelectorTask, starSelectorRegistry


class FlaggedStarSelectorConfig(BaseStarSelectorTask.ConfigClass):
    field = lsst.pex.config.Field(
        dtype=str, default="calib_psfUsed",
        doc="Name of a flag field that is True for stars that should be used."
    )


class FlaggedStarSelectorTask(BaseStarSelectorTask):
    """!
    A trivial StarSelector that simply uses an existing flag field to filter a SourceCatalog.

    This is most frequently used in steps that occur after the a PSF model has
    been built, to allow other procedures that need stars to use the set of
    stars used to determine the PSF.
    """

    usesMatches = False  # This selector does not require a match to an external catalog
    ConfigClass = FlaggedStarSelectorConfig

    def __init__(self, schema, **kwds):
        BaseStarSelectorTask.__init__(self, schema=schema, **kwds)
        self.key = schema.find(self.config.field).key

    def selectStars(self, exposure, sourceCat, matches=None):
        starCat = lsst.afw.table.SourceCatalog(sourceCat.table)
        for record in sourceCat:
            if record.get(self.key):
                starCat.append(record)
        return lsst.pipe.base.Struct(starCat=starCat)


starSelectorRegistry.register("flagged", FlaggedStarSelectorTask)
