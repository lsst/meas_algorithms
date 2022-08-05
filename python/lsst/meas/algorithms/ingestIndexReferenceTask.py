# This file is part of meas_algorithms.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


# TODO DM-31698: post-gen2 removal notes
# `DatasetConfig`, `ConvertReferenceCatalogBase`, and `ConvertReferenceCatalogConfig`
# should all be moved to to `convertReferenceCatalog.py` once gen2 butler
# has been removed.

__all__ = ["DatasetConfig", "ConvertReferenceCatalogBase", "ConvertReferenceCatalogConfig"]

import abc
import os.path

import astropy.units

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.geom
import lsst.sphgeom
import lsst.afw.table as afwTable
from lsst.daf.base import PropertyList
from .indexerRegistry import IndexerRegistry
from .readTextCatalogTask import ReadTextCatalogTask
from .loadReferenceObjects import ReferenceObjectLoader
from . import convertRefcatManager

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
    """Task runner for the reference catalog ingester (gen2 version).

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

        task.run(files)
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
        # TODO DM-31817: remove this default value.
        default='cal_ref_cat',
        doc="Name of this reference catalog to be used in the butler registry.",
    )
    indexer = IndexerRegistry.makeField(
        default='HTM',
        doc='Name of indexer algoritm to use.  Default is HTM',
    )


class ConvertReferenceCatalogConfig(pexConfig.Config):
    dataset_config = pexConfig.ConfigField(
        dtype=DatasetConfig,
        doc="Configuration for reading the ingested data",
    )
    n_processes = pexConfig.Field(
        dtype=int,
        doc=("Number of python processes to use when ingesting."),
        default=1
    )
    manager = pexConfig.ConfigurableField(
        target=convertRefcatManager.ConvertRefcatManager,
        doc="Multiprocessing manager to perform the actual conversion of values, file-by-file."
    )
    file_reader = pexConfig.ConfigurableField(
        target=ReadTextCatalogTask,
        doc='Task to use to read the files.  Default is to expect text files.'
    )
    ra_name = pexConfig.Field(
        dtype=str,
        doc="Name of RA column (values in decimal degrees)",
    )
    dec_name = pexConfig.Field(
        dtype=str,
        doc="Name of Dec column (values in decimal degrees)",
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
    coord_err_unit = pexConfig.Field(
        dtype=str,
        doc="Unit of RA/Dec error fields (astropy.unit.Unit compatible)",
        optional=True
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
        # gen3 refcats are all depth=7
        self.dataset_config.indexer['HTM'].depth = 7

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
        assertAllOrNone("ra_err_name", "dec_err_name", "coord_err_unit")
        if self.coord_err_unit is not None:
            result = astropy.units.Unit(self.coord_err_unit, parse_strict='silent')
            if isinstance(result, astropy.units.UnrecognizedUnit):
                msg = f"{self.coord_err_unit} is not a valid astropy unit string."
                raise pexConfig.FieldValidationError(ConvertReferenceCatalogConfig.coord_err_unit, self, msg)

        assertAllOrNone("epoch_name", "epoch_format", "epoch_scale")
        assertAllOrNone("pm_ra_name", "pm_dec_name")
        assertAllOrNone("pm_ra_err_name", "pm_dec_err_name")
        assertAllOrNone("parallax_name", "parallax_err_name")
        if self.pm_ra_err_name and not self.pm_ra_name:
            raise ValueError('"pm_ra/dec_name" must be specified if "pm_ra/dec_err_name" are specified')
        if (self.pm_ra_name or self.parallax_name) and not self.epoch_name:
            raise ValueError(
                '"epoch_name" must be specified if "pm_ra/dec_name" or "parallax_name" are specified')


class ConvertReferenceCatalogBase(pipeBase.Task, abc.ABC):
    """Base class for producing and loading indexed reference catalogs,
    shared between gen2 and gen3.

    This implements an indexing scheme based on hierarchical triangular
    mesh (HTM). The term index really means breaking the catalog into
    localized chunks called shards.  In this case each shard contains
    the entries from the catalog in a single HTM trixel

    For producing catalogs this task makes the following assumptions
    about the input catalogs:
    - RA, Dec are in decimal degrees.
    - Epoch is available in a column, in a format supported by astropy.time.Time.
    - There are no off-diagonal covariance terms, such as covariance
      between RA and Dec, or between PM RA and PM Dec. Support for such
     covariance would have to be added to to the config, including consideration
     of the units in the input catalog.
    """
    canMultiprocess = False
    ConfigClass = ConvertReferenceCatalogConfig

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.indexer = IndexerRegistry[self.config.dataset_config.indexer.name](
            self.config.dataset_config.indexer.active)
        self.makeSubtask('file_reader')

    def run(self, inputFiles):
        """Index a set of files comprising a reference catalog.

        Outputs are persisted in the butler repository.

        Parameters
        ----------
        inputFiles : `list`
            A list of file paths to read.
        """
        self._preRun()
        schema, key_map = self._saveMasterSchema(inputFiles[0])
        # create an HTM we can interrogate about pixel ids
        htm = lsst.sphgeom.HtmPixelization(self.indexer.htm.get_depth())
        filenames = self._getButlerFilenames(htm)
        worker = self.config.manager.target(filenames,
                                            self.config,
                                            self.file_reader,
                                            self.indexer,
                                            schema,
                                            key_map,
                                            htm.universe()[0],
                                            addRefCatMetadata,
                                            self.log)
        result = worker.run(inputFiles)

        self._persistConfig()
        self._postRun(result)

    def _preRun(self):
        """Any setup that has to be performed at the start of ``run``, but that
        cannot be performed during ``__init__`` (e.g. making directories).
        """
        pass

    def _postRun(self, result):
        """Any tasks that have to happen at the end of ``run``.

        Parameters
        ----------
        result
            The result returned from``worker.run()``.
        """
        pass

    def _getButlerFilenames(self, htm):
        """Get filenames from the butler for each output htm pixel.

        Parameters
        ----------
        htm : `lsst.sphgeom.HtmPixelization`
            The HTM pixelization scheme to be used to build filenames.

        Returns
        -------
        filenames : `list [str]`
            List of filenames to write each HTM pixel to.
        """
        filenames = {}
        start, end = htm.universe()[0]
        # path manipulation because butler.get() per pixel will take forever
        path = self._getOnePixelFilename(start)
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
        # make a schema with the standard fields
        schema = ReferenceObjectLoader.makeMinimalSchema(
            filterNameList=self.config.mag_column_list,
            addCentroid=False,
            addIsPhotometric=bool(self.config.is_photometric_name),
            addIsResolved=bool(self.config.is_resolved_name),
            addIsVariable=bool(self.config.is_variable_name),
            coordErrDim=2 if bool(self.config.ra_err_name) else 0,
            addProperMotion=2 if bool(self.config.pm_ra_name) else 0,
            properMotionErrDim=2 if bool(self.config.pm_ra_err_name) else 0,
            addParallax=bool(self.config.parallax_name),
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

    def _saveMasterSchema(self, filename):
        """Generate and save the master catalog schema.

        Parameters
        ----------
        filename : `str`
            An input file to read to get the input dtype.
        """
        arr = self.file_reader.run(filename)
        schema, key_map = self.makeSchema(arr.dtype)

        catalog = afwTable.SimpleCatalog(schema)
        addRefCatMetadata(catalog)
        self._writeMasterSchema(catalog)
        return schema, key_map

    @abc.abstractmethod
    def _getOnePixelFilename(self, start):
        """Return one example filename to help construct the rest of the
        per-htm pixel filenames.

        Parameters
        ----------
        start : `int`
            The first HTM index in this HTM pixelization.

        Returns
        -------
        filename : `str`
            Path to a single file that would be written to the output location.
        """
        pass

    @abc.abstractmethod
    def _persistConfig(self):
        """Write the config that was used to generate the refcat.
        """
        pass

    @abc.abstractmethod
    def _writeMasterSchema(self, catalog):
        """Butler put the master catalog schema.

        Parameters
        ----------
        catalog : `lsst.afw.table.SimpleCatalog`
            An empty catalog with a fully-defined schema that matches the
            schema used in each of the HTM pixel files.
        """
        pass
