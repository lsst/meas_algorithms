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
import esutil


class HtmIndexer:
    """Manage a spatial index of hierarchical triangular mesh (HTM)
    shards.

    Parameters
    ----------
    depth : `int`
        Depth of the HTM hierarchy to construct.
    """
    def __init__(self, depth=8):
        self.htm = esutil.htm.HTM(depth)

    def getShardIds(self, ctrCoord, radius):
        """Get the IDs of all shards that touch a circular aperture.

        Parameters
        ----------
        ctrCoord : `lsst.geom.SpherePoint`
            ICRS center of the aperture

        radius : `lsst.geom.Angle`
            object of the aperture radius

        Returns
        -------
        shardIdList : `pipeBase.Struct`
            Struct with the list of shards, shards, and a boolean arry, boundary_mask,
            indicating whether the shard touches the boundary (True) or is fully contained (False).

        isOnBoundary : `pipeBase.Struct`
            Struct with the list of shards, shards, and a boolean arry, boundary_mask,
            indicating whether the shard touches the boundary (True) or is fully contained (False).
        """
        shardIdList = self.htm.intersect(ctrCoord.getLongitude().asDegrees(),
                                         ctrCoord.getLatitude().asDegrees(),
                                         radius.asDegrees(), inclusive=True)
        coveredShardIdList = self.htm.intersect(ctrCoord.getLongitude().asDegrees(),
                                                ctrCoord.getLatitude().asDegrees(),
                                                radius.asDegrees(), inclusive=False)
        isOnBoundary = (shardId not in coveredShardIdList for shardId in shardIdList)
        return shardIdList, isOnBoundary

    def indexPoints(self, raList, decList):
        """Generate trixel ids for each row in an input file

        Parameters
        ----------
        ra_list: `float`
        List of RA coordinate in degrees
        dec_list: `float`
        List of Dec coordinate in degrees

        Returns
        -------
        htm.lookup_id : `int`
        A list of pixel ids using `ra_list` and `dec_list` as parameters.
        """
        return self.htm.lookup_id(raList, decList)

    @staticmethod
    def makeDataId(shardId, datasetName):
        """Make a data id from a shard ID.

        Parameters
        ----------
        shardId : `int`
            ID of shard in question.
        datasetName : `str`
            Name of dataset to use.

        Returns
        -------
        dataId : `dict`
            Data ID for shard.
        """
        if shardId is None:
            # NoneType doesn't format, so make dummy pixel
            shardId = 0
        return {'pixel_id': shardId, 'name': datasetName}
