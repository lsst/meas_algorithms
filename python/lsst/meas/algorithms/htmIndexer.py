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
from builtins import object
import esutil


class HtmIndexer(object):

    def __init__(self, depth=8):
        """!Construct the indexer object

        @param[in] depth  depth of the hierarchy to construct
        """
        self.htm = esutil.htm.HTM(depth)

    def get_pixel_ids(self, ctrCoord, radius):
        """!Get all shards that touch a circular aperture

        @param[in] ctrCoord  afwCoord.Coord object of the center of the aperture
        @param[in] radius  afwGeom.Angle object of the aperture radius
        @param[out] A pipeBase.Struct with the list of shards, shards, and a boolean arry, boundary_mask,
                    indicating whether the shard touches the boundary (True) or is fully contained (False).
        """
        pixel_id_list = self.htm.intersect(ctrCoord.getLongitude().asDegrees(),
                                           ctrCoord.getLatitude().asDegrees(),
                                           radius.asDegrees(), inclusive=True)
        covered_pixel_id_list = self.htm.intersect(ctrCoord.getLongitude().asDegrees(),
                                                   ctrCoord.getLatitude().asDegrees(),
                                                   radius.asDegrees(), inclusive=False)
        is_on_boundary = (pixel_id not in covered_pixel_id_list for pixel_id in pixel_id_list)
        return pixel_id_list, is_on_boundary

    def index_points(self, ra_list, dec_list):
        """!Generate trixel ids for each row in an input file

        @param[in] ra_list  List of RA coordinate in degrees
        @param[in] dec_list  List of Dec coordinate in degrees
        @param[out] A list of pixel ids
        """
        return self.htm.lookup_id(ra_list, dec_list)

    @staticmethod
    def make_data_id(pixel_id, dataset_name):
        """!Make a data id.  Meant to be overridden.
        @param[in] pixel_id  An identifier for the pixel in question.
        @param[in] dataset_name  Name of the dataset to use.
        @param[out] dataId (dictionary)
        """
        if pixel_id is None:
            # NoneType doesn't format, so make dummy pixel
            pixel_id = 0
        return {'pixel_id': pixel_id, 'name': dataset_name}
