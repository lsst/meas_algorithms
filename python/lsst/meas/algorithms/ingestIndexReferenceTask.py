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

__all__ = ["IngestIndexedReferenceConfig", "IngestIndexedReferenceTask", "DatasetConfig"]

import os.path

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.geom
import lsst.sphgeom
import lsst.afw.table as afwTable
from lsst.daf.base import PropertyList
from .indexerRegistry import IndexerRegistry
from .readTextCatalogTask import ReadTextCatalogTask
from .loadReferenceObjects import LoadReferenceObjectsTask
from .ingestIndexManager import IngestIndexManager

# The most recent Indexed Reference Catalog on-disk format version.
LATEST_FORMAT_VERSION = 1


def addRefCatMetadata(catalog):
    """Add metadata to a new (not yet populated) reference catalog.

    Parameters
    ----------
    catalog : `lsst.afw.table.SimpleCatalog`
        Catalog to which metadata should be attached.  Will be modified
        in-place.
    """
    md = catalog.getMetadata()
    if md is None:
        md = PropertyList()
    md.set("REFCAT_FORMAT_VERSION", LATEST_FORMAT_VERSION)
    catalog.setMetadata(md)


class IngestReferenceRunner(pipeBase.TaskRunner):
    """Task runner for the reference catalog ingester

    Data IDs are ignored so the runner should just run the task on the parsed command.
    """

    def run(self, parsedCmd):
        """Run the task.

        Several arguments need to be collected to send on to the task methods.

        Parameters
        ----------
        parsedCmd : `argparse.Namespace`
            Parsed command.

        Returns
        -------
        results : `lsst.pipe.base.Struct` or `None`
            A empty struct if self.doReturnResults, else None
        """
        files = parsedCmd.files
        butler = parsedCmd.butler
        task = self.TaskClass(config=self.config, log=self.log, butler=butler)
        task.writeConfig(parsedCmd.butler, clobber=self.clobberConfig, doBackup=self.doBackup)

        task.createIndexedCatalog(files)
        if self.doReturnResults:
            return pipeBase.Struct()


class DatasetConfig(pexConfig.Config):
    """The description of the on-disk storage format for the persisted
    reference catalog.
    """
    format_version = pexConfig.Field(
        dtype=int,
        doc="Version number of the persisted on-disk storage format."
        "\nVersion 0 had Jy as flux units (default 0 for unversioned catalogs)."
        "\nVersion 1 had nJy as flux units.",
        default=0  # This needs to always be 0, so that unversioned catalogs are interpreted as version 0.
    )
    ref_dataset_name = pexConfig.Field(
        dtype=str,
        default='cal_ref_cat',
        doc='String to pass to the butler to retrieve persisted files.',
    )
    indexer = IndexerRegistry.makeField(
        default='HTM',
        doc='Name of indexer algoritm to use.  Default is HTM',
    )


class IngestIndexedReferenceConfig(pexConfig.Config):
    dataset_config = pexConfig.ConfigField(
        dtype=DatasetConfig,
        doc="Configuration for reading the ingested data",
    )
    n_processes = pexConfig.Field(
        dtype=int,
        doc=("Number of python processes to use when ingesting."),
        default=1
    )
    file_reader = pexConfig.ConfigurableField(
        target=ReadTextCatalogTask,
        doc='Task to use to read the files.  Default is to expect text files.'
    )
    ra_name = pexConfig.Field(
        dtype=str,
        doc="Name of RA column",
    )
    dec_name = pexConfig.Field(
        dtype=str,
        doc="Name of Dec column",
    )
    ra_err_name = pexConfig.Field(
        dtype=str,
        doc="Name of RA error column",
        optional=True,
    )
    dec_err_name = pexConfig.Field(
        dtype=str,
        doc="Name of Dec error column",
        optional=True,
    )
    mag_column_list = pexConfig.ListField(
        dtype=str,
        doc="The values in the reference catalog are assumed to be in AB magnitudes. "
            "List of column names to use for photometric information.  At least one entry is required."
    )
    mag_err_column_map = pexConfig.DictField(
        keytype=str,
        itemtype=str,
        default={},
        doc="A map of magnitude column name (key) to magnitude error column (value)."
    )
    is_photometric_name = pexConfig.Field(
        dtype=str,
        optional=True,
        doc='Name of column stating if satisfactory for photometric calibration (optional).'
    )
    is_resolved_name = pexConfig.Field(
        dtype=str,
        optional=True,
        doc='Name of column stating if the object is resolved (optional).'
    )
    is_variable_name = pexConfig.Field(
        dtype=str,
        optional=True,
        doc='Name of column stating if the object is measured to be variable (optional).'
    )
    id_name = pexConfig.Field(
        dtype=str,
        optional=True,
        doc='Name of column to use as an identifier (optional).'
    )
    pm_ra_name = pexConfig.Field(
        dtype=str,
        doc="Name of proper motion RA column",
        optional=True,
    )
    pm_dec_name = pexConfig.Field(
        dtype=str,
        doc="Name of proper motion Dec column",
        optional=True,
    )
    pm_ra_err_name = pexConfig.Field(
        dtype=str,
        doc="Name of proper motion RA error column",
        optional=True,
    )
    pm_dec_err_name = pexConfig.Field(
        dtype=str,
        doc="Name of proper motion Dec error column",
        optional=True,
    )
    pm_scale = pexConfig.Field(
        dtype=float,
        doc="Scale factor by which to multiply proper motion values to obtain units of milliarcsec/year",
        default=1.0,
    )
    parallax_name = pexConfig.Field(
        dtype=str,
        doc="Name of parallax column",
        optional=True,
    )
    parallax_err_name = pexConfig.Field(
        dtype=str,
        doc="Name of parallax error column",
        optional=True,
    )
    parallax_scale = pexConfig.Field(
        dtype=float,
        doc="Scale factor by which to multiply parallax values to obtain units of milliarcsec",
        default=1.0,
    )
    epoch_name = pexConfig.Field(
        dtype=str,
        doc="Name of epoch column",
        optional=True,
    )
    epoch_format = pexConfig.Field(
        dtype=str,
        doc="Format of epoch column: any value accepted by astropy.time.Time, e.g. 'iso' or 'unix'",
        optional=True,
    )
    epoch_scale = pexConfig.Field(
        dtype=str,
        doc="Scale of epoch column: any value accepted by astropy.time.Time, e.g. 'utc'",
        optional=True,
    )
    extra_col_names = pexConfig.ListField(
        dtype=str,
        default=[],
        doc='Extra columns to add to the reference catalog.'
    )

    def setDefaults(self):
        # Newly ingested reference catalogs always have the latest format_version.
        self.dataset_config.format_version = LATEST_FORMAT_VERSION

    def validate(self):
        pexConfig.Config.validate(self)

        def assertAllOrNone(*names):
            """Raise ValueError unless all the named fields are set or are
            all none (or blank)
            """
            setNames = [name for name in names if bool(getattr(self, name))]
            if len(setNames) in (len(names), 0):
                return
            prefix = "Both or neither" if len(names) == 2 else "All or none"
            raise ValueError("{} of {} must be set, but only {} are set".format(
                prefix, ", ".join(names), ", ".join(setNames)))

        if not (self.ra_name and self.dec_name and self.mag_column_list):
            raise ValueError(
                "ra_name and dec_name and at least one entry in mag_column_list must be supplied.")
        if self.mag_err_column_map and set(self.mag_column_list) != set(self.mag_err_column_map.keys()):
            raise ValueError(
                "mag_err_column_map specified, but keys do not match mag_column_list: {} != {}".format(
                    sorted(self.mag_err_column_map.keys()), sorted(self.mag_column_list)))
        assertAllOrNone("ra_err_name", "dec_err_name")
        assertAllOrNone("epoch_name", "epoch_format", "epoch_scale")
        assertAllOrNone("pm_ra_name", "pm_dec_name")
        assertAllOrNone("pm_ra_err_name", "pm_dec_err_name")
        if self.pm_ra_err_name and not self.pm_ra_name:
            raise ValueError('"pm_ra/dec_name" must be specified if "pm_ra/dec_err_name" are specified')
        if (self.pm_ra_name or self.parallax_name) and not self.epoch_name:
            raise ValueError(
                '"epoch_name" must be specified if "pm_ra/dec_name" or "parallax_name" are specified')


class IngestIndexedReferenceTask(pipeBase.CmdLineTask):
    """Class for producing and loading indexed reference catalogs.

    This implements an indexing scheme based on hierarchical triangular
    mesh (HTM). The term index really means breaking the catalog into
    localized chunks called shards.  In this case each shard contains
    the entries from the catalog in a single HTM trixel

    For producing catalogs this task makes the following assumptions
    about the input catalogs:
    - RA, Dec, RA error and Dec error are all in decimal degrees.
    - Epoch is available in a column, in a format supported by astropy.time.Time.
    - There are no off-diagonal covariance terms, such as covariance
        between RA and Dec, or between PM RA and PM Dec. Gaia is a well
        known example of a catalog that has such terms, and thus should not
        be ingested with this task.

    Parameters
    ----------
    butler : `lsst.daf.persistence.Butler`
        Data butler for reading and writing catalogs
    """
    canMultiprocess = False
    ConfigClass = IngestIndexedReferenceConfig
    RunnerClass = IngestReferenceRunner
    _DefaultName = 'IngestIndexedReferenceTask'

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser.

        This returns a standard parser with an extra "files" argument.
        """
        parser = pipeBase.InputOnlyArgumentParser(name=cls._DefaultName)
        parser.add_argument("files", nargs="+", help="Names of files to index")
        return parser

    def __init__(self, *args, butler=None, **kwargs):
        self.butler = butler
        super().__init__(*args, **kwargs)
        self.indexer = IndexerRegistry[self.config.dataset_config.indexer.name](
            self.config.dataset_config.indexer.active)
        self.makeSubtask('file_reader')

    def createIndexedCatalog(self, inputFiles):
        """Index a set of files comprising a reference catalog.

        Outputs are persisted in the butler repository.

        Parameters
        ----------
        inputFiles : `list`
            A list of file paths to read.
        """
        schema, key_map = self._saveMasterSchema(inputFiles[0])
        # create an HTM we can interrogate about pixel ids
        htm = lsst.sphgeom.HtmPixelization(self.indexer.htm.get_depth())
        filenames = self._getButlerFilenames(htm)
        worker = IngestIndexManager(filenames,
                                    self.config,
                                    self.file_reader,
                                    self.indexer,
                                    schema,
                                    key_map,
                                    htm.universe()[0],
                                    addRefCatMetadata,
                                    self.log)
        worker.run(inputFiles)

        # write the config that was used to generate the refcat
        dataId = self.indexer.makeDataId(None, self.config.dataset_config.ref_dataset_name)
        self.butler.put(self.config.dataset_config, 'ref_cat_config', dataId=dataId)

    def _saveMasterSchema(self, filename):
        """Generate and save the master catalog schema.

        Parameters
        ----------
        filename : `str`
            An input file to read to get the input dtype.
        """
        arr = self.file_reader.run(filename)
        schema, key_map = self.makeSchema(arr.dtype)
        dataId = self.indexer.makeDataId('master_schema',
                                         self.config.dataset_config.ref_dataset_name)

        catalog = afwTable.SimpleCatalog(schema)
        addRefCatMetadata(catalog)
        self.butler.put(catalog, 'ref_cat', dataId=dataId)
        return schema, key_map

    def _getButlerFilenames(self, htm):
        """Get filenames from the butler for each output pixel."""
        filenames = {}
        start, end = htm.universe()[0]
        # path manipulation because butler.get() per pixel will take forever
        dataId = self.indexer.makeDataId(start, self.config.dataset_config.ref_dataset_name)
        path = self.butler.get('ref_cat_filename', dataId=dataId)[0]
        base = os.path.join(os.path.dirname(path), "%d"+os.path.splitext(path)[1])
        for pixelId in range(start, end):
            filenames[pixelId] = base % pixelId

        return filenames

    def makeSchema(self, dtype):
        """Make the schema to use in constructing the persisted catalogs.

        Parameters
        ----------
        dtype : `numpy.dtype`
            Data type describing each entry in ``config.extra_col_names``
            for the catalogs being ingested.

        Returns
        -------
        schemaAndKeyMap : `tuple` of (`lsst.afw.table.Schema`, `dict`)
            A tuple containing two items:
            - The schema for the output source catalog.
            - A map of catalog keys to use in filling the record
        """
        self.config.validate()  # just to be sure

        # make a schema with the standard fields
        schema = LoadReferenceObjectsTask.makeMinimalSchema(
            filterNameList=self.config.mag_column_list,
            addCentroid=False,
            addIsPhotometric=bool(self.config.is_photometric_name),
            addIsResolved=bool(self.config.is_resolved_name),
            addIsVariable=bool(self.config.is_variable_name),
            coordErrDim=2 if bool(self.config.ra_err_name) else 0,
            addProperMotion=2 if bool(self.config.pm_ra_name) else 0,
            properMotionErrDim=2 if bool(self.config.pm_ra_err_name) else 0,
            addParallax=bool(self.config.parallax_name),
            addParallaxErr=bool(self.config.parallax_err_name),
        )
        keysToSkip = set(("id", "centroid_x", "centroid_y", "hasCentroid"))
        key_map = {fieldName: schema[fieldName].asKey() for fieldName in schema.getOrderedNames()
                   if fieldName not in keysToSkip}

        def addField(name):
            if dtype[name].kind == 'U':
                # dealing with a string like thing.  Need to get type and size.
                at_size = dtype[name].itemsize
                return schema.addField(name, type=str, size=at_size)
            else:
                at_type = dtype[name].type
                return schema.addField(name, at_type)

        for col in self.config.extra_col_names:
            key_map[col] = addField(col)
        return schema, key_map
