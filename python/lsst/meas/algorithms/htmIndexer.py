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
    """Manage a spatial index using a hierarchical triangular mesh (HTM).

    Parameters
    ----------
    depth : `int`
        Depth of the HTM hierarchy to construct.
    """
    def __init__(self, depth=8):
        self.htm = esutil.htm.HTM(depth)

    def getShardIds(self, ctrCoord, radius):
        """Get all shards that touch a circular aperture

        Parameters
        ----------
        ctrCoord : `lsst.geom.SpherePoint`
            ICRS center of search region.
        radius : `lsst.geom.Angle`
            Radius of search region.

        Returns
        -------
        results : `tuple`
            A tuple containing:

            - shardIdList : `list` of `int`
                List of shard IDs
            - isOnBoundary : `list` of `bool`
                For each shard in ``shardIdList`` is the shard on the
                boundary?
        """
        pixel_id_list = self.htm.intersect(ctrCoord.getLongitude().asDegrees(),
                                           ctrCoord.getLatitude().asDegrees(),
                                           radius.asDegrees(), inclusive=True)
        covered_pixel_id_list = self.htm.intersect(ctrCoord.getLongitude().asDegrees(),
                                                   ctrCoord.getLatitude().asDegrees(),
                                                   radius.asDegrees(), inclusive=False)
        is_on_boundary = (pixel_id not in covered_pixel_id_list for pixel_id in pixel_id_list)
        return pixel_id_list, is_on_boundary

    def indexPoints(self, ra_list, dec_list):
        """Generate shard IDs for sky positions.

        Parameters
        ----------
        ra_list : `list` of `float`
            List of right ascensions, in degrees.
        dec_list : `list` of `float`
            List of declinations, in degrees.

        Returns
        -------
        shardIds : `list` of `int`
            List of shard IDs
        """
        return self.htm.lookup_id(ra_list, dec_list)

    @staticmethod
    def makeDataId(pixel_id, dataset_name):
        """Make a data id from a shard ID.

        Meant to be overridden.

        Parameters
        ----------
        pixel_id : `int`
            ID of shard in question.
        dataset_name : `str`
            Name of dataset to use.

        Returns
        -------
        dataId : `dict`
            Data ID for shard.
        """
        if pixel_id is None:
            # NoneType doesn't format, so make dummy pixel
            pixel_id = 0
        return {'pixel_id': pixel_id, 'name': dataset_name}
