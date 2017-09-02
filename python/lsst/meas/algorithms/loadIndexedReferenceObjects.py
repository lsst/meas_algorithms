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

__all__ = ["LoadIndexedReferenceObjectsConfig", "LoadIndexedReferenceObjectsTask"]

from builtins import zip
from lsst.meas.algorithms import getRefFluxField, LoadReferenceObjectsTask, LoadReferenceObjectsConfig
import lsst.afw.table as afwTable
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from .indexerRegistry import IndexerRegistry


class LoadIndexedReferenceObjectsConfig(LoadReferenceObjectsConfig):
    ref_dataset_name = pexConfig.Field(
        dtype=str,
        default='cal_ref_cat',
        doc='Name of the ingested reference dataset'
    )


class LoadIndexedReferenceObjectsTask(LoadReferenceObjectsTask):
    ConfigClass = LoadIndexedReferenceObjectsConfig
    _DefaultName = 'LoadIndexedReferenceObjectsTask'

    def __init__(self, butler, *args, **kwargs):
        LoadReferenceObjectsTask.__init__(self, *args, **kwargs)
        dataset_config = butler.get("ref_cat_config", name=self.config.ref_dataset_name, immediate=True)
        self.indexer = IndexerRegistry[dataset_config.indexer.name](dataset_config.indexer.active)
        # This needs to come from the loader config, not the dataset_config since directory aliases can
        # change the path where the shards are found.
        self.ref_dataset_name = self.config.ref_dataset_name
        self.butler = butler

    @pipeBase.timeMethod
    def loadSkyCircle(self, ctrCoord, radius, filterName=None):
        """!Load reference objects that overlap a circular sky region

        @param[in] ctrCoord  center of search region (an lsst.afw.geom.Coord)
        @param[in] radius  radius of search region (an lsst.afw.geom.Angle)
        @param[in] filterName  name of filter, or None for the default filter;
            used for flux values in case we have flux limits (which are not yet implemented)

        @return an lsst.pipe.base.Struct containing:
        - refCat a catalog of reference objects with the
            \link meas_algorithms_loadReferenceObjects_Schema standard schema \endlink
            as documented in LoadReferenceObjects, including photometric, resolved and variable;
            hasCentroid is False for all objects.
        - fluxField = name of flux field for specified filterName.  None if refCat is None.
        """
        id_list, boundary_mask = self.indexer.get_pixel_ids(ctrCoord, radius)
        shards = self.get_shards(id_list)
        refCat = self.butler.get('ref_cat',
                                 dataId=self.indexer.make_data_id('master_schema', self.ref_dataset_name),
                                 immediate=True)
        self._addFluxAliases(refCat.schema)
        fluxField = getRefFluxField(schema=refCat.schema, filterName=filterName)
        for shard, is_on_boundary in zip(shards, boundary_mask):
            if shard is None:
                continue
            if is_on_boundary:
                refCat.extend(self._trim_to_circle(shard, ctrCoord, radius))
            else:
                refCat.extend(shard)

        # make sure catalog is contiguous
        if not refCat.isContiguous():
            refCat = refCat.copy()

        # add and initialize centroid and hasCentroid fields (these are added
        # after loading to avoid wasting space in the saved catalogs)
        # the new fields are automatically initialized to (nan, nan) and False
        # so no need to set them explicitly
        mapper = afwTable.SchemaMapper(refCat.schema, True)
        mapper.addMinimalSchema(refCat.schema, True)
        mapper.editOutputSchema().addField("centroid_x", type=float)
        mapper.editOutputSchema().addField("centroid_y", type=float)
        mapper.editOutputSchema().addField("hasCentroid", type="Flag")
        expandedCat = afwTable.SimpleCatalog(mapper.getOutputSchema())
        expandedCat.extend(refCat, mapper=mapper)
        del refCat  # avoid accidentally returning the unexpanded reference catalog

        # return reference catalog
        return pipeBase.Struct(
            refCat=expandedCat,
            fluxField=fluxField,
        )

    def get_shards(self, id_list):
        """!Get all shards that touch a circular aperture

        @param[in] id_list  A list of integer pixel ids
        @param[out] a list of SourceCatalogs for each pixel, None if not data exists
        """
        shards = []
        for pixel_id in id_list:
            if self.butler.datasetExists('ref_cat',
                                         dataId=self.indexer.make_data_id(pixel_id, self.ref_dataset_name)):
                shards.append(self.butler.get('ref_cat',
                                              dataId=self.indexer.make_data_id(pixel_id,
                                                                               self.ref_dataset_name),
                                              immediate=True))
        return shards

    def _trim_to_circle(self, catalog_shard, ctrCoord, radius):
        """!Trim a catalog to a circular aperture.

        @param[in] catalog_shard  SourceCatalog to be trimmed
        @param[in] ctrCoord  afw.Coord to compare each record to
        @param[in] radius  afwGeom.Angle indicating maximume separation
        @param[out] a SourceCatalog constructed from records that fall in the circular aperture
        """
        temp_cat = type(catalog_shard)(catalog_shard.schema)
        for record in catalog_shard:
            if record.getCoord().angularSeparation(ctrCoord) < radius:
                temp_cat.append(record)
        return temp_cat
