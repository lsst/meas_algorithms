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
from __future__ import absolute_import, division

__all__ = ["IndexerRegistry"]

from lsst.pex.config import Config, makeRegistry, Field
from .htmIndexer import HtmIndexer

IndexerRegistry = makeRegistry(
    """Registry of indexing algorithms
    """
)


class HtmIndexerConfig(Config):
    depth = Field(
        doc = """Depth of the HTM tree to make.  Default is depth=7 which gives
              ~ 0.3 sq. deg. per trixel.""",
        dtype = int,
        default = 7,
    )


def makeHtmIndexer(config):
    """Make an HtmIndexer
    """
    return HtmIndexer(depth=config.depth)


makeHtmIndexer.ConfigClass = HtmIndexerConfig
IndexerRegistry.register("HTM", makeHtmIndexer)
