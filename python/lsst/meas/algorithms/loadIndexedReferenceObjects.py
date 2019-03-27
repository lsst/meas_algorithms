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

from .loadReferenceObjects import hasNanojanskyFluxUnits, convertToNanojansky, getFormatVersionFromRefCat
from lsst.meas.algorithms import getRefFluxField, LoadReferenceObjectsTask, LoadReferenceObjectsConfig
import lsst.afw.table as afwTable
import lsst.geom
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
        self.dataset_config = butler.get("ref_cat_config", name=self.config.ref_dataset_name, immediate=True)
        self.indexer = IndexerRegistry[self.dataset_config.indexer.name](self.dataset_config.indexer.active)
        # This needs to come from the loader config, not the dataset_config since directory aliases can
        # change the path where the shards are found.
        self.ref_dataset_name = self.config.ref_dataset_name
        self.butler = butler

    @pipeBase.timeMethod
    def loadSkyCircle(self, ctrCoord, radius, filterName=None, epoch=None):
        shardIdList, isOnBoundaryList = self.indexer.getShardIds(ctrCoord, radius)
        shards = self.getShards(shardIdList)
        refCat = self.butler.get('ref_cat',
                                 dataId=self.indexer.makeDataId('master_schema', self.ref_dataset_name),
                                 immediate=True)

        # load the catalog, one shard at a time
        for shard, isOnBoundary in zip(shards, isOnBoundaryList):
            if shard is None:
                continue
            if isOnBoundary:
                refCat.extend(self._trimToCircle(shard, ctrCoord, radius))
            else:
                refCat.extend(shard)

        # apply proper motion corrections
        if epoch is not None and "pm_ra" in refCat.schema:
            # check for a catalog in a non-standard format
            if isinstance(refCat.schema["pm_ra"].asKey(), lsst.afw.table.KeyAngle):
                self.applyProperMotions(refCat, epoch)
            else:
                self.log.warn("Catalog pm_ra field is not an Angle; not applying proper motion")

        # update version=0 style refcats to have nJy fluxes
        if self.dataset_config.format_version == 0 or not hasNanojanskyFluxUnits(refCat.schema):
            self.log.warn("Found version 0 reference catalog with old style units in schema.")
            self.log.warn("run `meas_algorithms/bin/convert_refcat_to_nJy.py` to convert fluxes to nJy.")
            self.log.warn("See RFC-575 for more details.")
            refCat = convertToNanojansky(refCat, self.log)
        else:
            # For version >= 1, the version should be in the catalog header,
            # too, and should be consistent with the version in the config.
            catVersion = getFormatVersionFromRefCat(refCat)
            if catVersion != self.dataset_config.format_version:
                raise RuntimeError(f"Format version in reference catalog ({catVersion}) does "
                                   f"not match configuration ({self.dataset_config.format_version})")

        self._addFluxAliases(refCat.schema)
        fluxField = getRefFluxField(schema=refCat.schema, filterName=filterName)

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

    def getShards(self, shardIdList):
        """Get shards by ID.

        Parameters
        ----------
        shardIdList : `list` of `int`
            A list of integer shard ids.

        Returns
        -------
        catalogs : `list` of `lsst.afw.table.SimpleCatalog`
            A list of reference catalogs, one for each entry in shardIdList.
        """
        shards = []
        for shardId in shardIdList:
            if self.butler.datasetExists('ref_cat',
                                         dataId=self.indexer.makeDataId(shardId, self.ref_dataset_name)):
                shards.append(self.butler.get('ref_cat',
                                              dataId=self.indexer.makeDataId(shardId, self.ref_dataset_name),
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
        tempCat = type(refCat)(refCat.schema)
        for record in refCat:
            if record.getCoord().separation(ctrCoord) < radius:
                tempCat.append(record)
        return tempCat
