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

__all__ = ["LoadIndexedReferenceObjectsConfig", "LoadIndexedReferenceObjectsTask"]

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
    """Load reference objects from an indexed catalog ingested by
    IngestIndexReferenceTask.

    Parameters
    ----------
    butler : `lsst.daf.persistence.Butler`
        Data butler for reading catalogs
    """
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
    def loadSkyCircle(self, ctrCoord, radius, filterName=None, epoch=None):
        idList, boundary_mask = self.indexer.get_pixel_ids(ctrCoord, radius)
        shards = self.getShards(idList)
        refCat = self.butler.get('ref_cat',
                                 dataId=self.indexer.make_data_id('master_schema', self.ref_dataset_name),
                                 immediate=True)
        self._addFluxAliases(refCat.schema)
        fluxField = getRefFluxField(schema=refCat.schema, filterName=filterName)
        for shard, is_on_boundary in zip(shards, boundary_mask):
            if shard is None:
                continue
            if is_on_boundary:
                refCat.extend(self._trimToCircle(shard, ctrCoord, radius))
            else:
                refCat.extend(shard)

        if epoch is not None:
            self.applyProperMotions(refCat, epoch)

        # add and initialize centroid and hasCentroid fields (these are
        # added after loading to avoid wasting space in the saved catalogs)
        # the new fields are automatically initialized to (nan, nan) and
        # False so no need to set them explicitly
        mapper = afwTable.SchemaMapper(refCat.schema, True)
        mapper.addMinimalSchema(refCat.schema, True)
        mapper.editOutputSchema().addField("centroid_x", type=float)
        mapper.editOutputSchema().addField("centroid_y", type=float)
        mapper.editOutputSchema().addField("hasCentroid", type="Flag")
        expandedCat = afwTable.SimpleCatalog(mapper.getOutputSchema())
        expandedCat.extend(refCat, mapper=mapper)
        del refCat  # avoid accidentally returning the unexpanded ref cat

        # make sure catalog is contiguous
        if not expandedCat.isContiguous():
            expandedCat = expandedCat.copy(True)

        # return reference catalog
        return pipeBase.Struct(
            refCat=expandedCat,
            fluxField=fluxField,
        )

    def getShards(self, idList):
        """Get shards by ID.

        Parameters
        ----------
        idList : `list` of `int`
            A list of integer shard ids.

        Returns
        -------
        catalogs : `list` of `lsst.afw.table.SimpleCatalog`
            A list of reference catalogs, one for each entry in idList.
        """
        shards = []
        for shardId in idList:
            if self.butler.datasetExists('ref_cat',
                                         dataId=self.indexer.make_data_id(shardId, self.ref_dataset_name)):
                shards.append(self.butler.get('ref_cat',
                                              dataId=self.indexer.make_data_id(shardId,
                                                                               self.ref_dataset_name),
                                              immediate=True))
        return shards

    def _trimToCircle(self, refCat, ctrCoord, radius):
        """Trim a reference catalog to a circular aperture.

        Parameters
        ----------
        refCat : `lsst.afw.table.SimpleCatalog`
            Reference catalog to be trimmed.
        ctrCoord : `lsst.geom.SpherePoint`
            ICRS center of search region.
        radius : `lsst.geom.Angle`
            Radius of search region.

        Returns
        -------
        catalog : `lsst.afw.table.SimpleCatalog`
            Catalog containing objects that fall in the circular aperture.
        """
        temp_cat = type(refCat)(refCat.schema)
        for record in refCat:
            if record.getCoord().separation(ctrCoord) < radius:
                temp_cat.append(record)
        return temp_cat
