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

__all__ = ["ConvertRefcatManager", "ConvertGaiaManager"]

from ctypes import c_int
import os.path
import itertools
import multiprocessing

import astropy.time
import astropy.units as u
import numpy as np

import lsst.sphgeom
import lsst.afw.table as afwTable
from lsst.afw.image import fluxErrFromABMagErr
import lsst.pex.config as pexConfig


# global shared counter to keep track of source ids
# (multiprocess sharing is most easily done with a global)
COUNTER = multiprocessing.Value(c_int, 0)
# global shared counter to keep track of number of files processed.
FILE_PROGRESS = multiprocessing.Value(c_int, 0)


class ConvertRefcatManagerConfig(pexConfig.Config):
    """Placeholder for ConfigurableField validation; refcat convert is
    configured by the parent convert Task.
    """
    pass


class ConvertRefcatManager:
    """
    Convert a reference catalog from external files into the LSST HTM sharded
    format, using a multiprocessing Pool to speed up the work.

    Parameters
    ----------
    filenames : `dict` [`int`, `str`]
        The HTM pixel id and filenames to convert the catalog into.
    config : `lsst.meas.algorithms.ConvertReferenceCatalogConfig`
        The Task configuration holding the field names.
    file_reader : `lsst.pipe.base.Task`
        The file reader to use to load the files.
    indexer : `lsst.meas.algorithms.HtmIndexer`
        The class used to compute the HTM pixel per coordinate.
    schema : `lsst.afw.table.Schema`
        The schema of the output catalog.
    key_map : `dict` [`str`, `lsst.afw.table.Key`]
        The mapping from output field names to keys in the Schema.
    htmRange : `tuple` [`int`]
        The start and end HTM pixel ids.
    addRefCatMetadata : callable
        A function called to add extra metadata to each output Catalog.
    log : `lsst.log.Log` or `logging.Logger`
        The log to send messages to.
    """
    _flags = ['photometric', 'resolved', 'variable']
    _DefaultName = 'convertRefcatManager'
    ConfigClass = ConvertRefcatManagerConfig

    def __init__(self, filenames, config, file_reader, indexer,
                 schema, key_map, htmRange, addRefCatMetadata, log):
        self.filenames = filenames
        self.config = config
        self.file_reader = file_reader
        self.indexer = indexer
        self.schema = schema
        self.key_map = key_map
        self.htmRange = htmRange
        self.addRefCatMetadata = addRefCatMetadata
        self.log = log
        if self.config.coord_err_unit is not None:
            # cache this to speed up coordinate conversions
            self.coord_err_unit = u.Unit(self.config.coord_err_unit)

    def run(self, inputFiles):
        """Index a set of input files from a reference catalog, and write the
        output to the appropriate filenames, in parallel.

        Parameters
        ----------
        inputFiles : `list`
            A list of file paths to read data from.

        Returns
        -------
        output : `dict` [`int`, `str`]
            The htm ids and the filenames that were written to.
        """
        global COUNTER, FILE_PROGRESS
        self.nInputFiles = len(inputFiles)

        with multiprocessing.Manager() as manager:
            COUNTER.value = 0
            FILE_PROGRESS.value = 0
            fileLocks = manager.dict()
            self.log.info("Creating %s file locks.", self.htmRange[1] - self.htmRange[0])
            for i in range(self.htmRange[0], self.htmRange[1]):
                fileLocks[i] = manager.Lock()
            self.log.info("File locks created.")
            with multiprocessing.Pool(self.config.n_processes) as pool:
                result = pool.starmap(self._convertOneFile, zip(inputFiles, itertools.repeat(fileLocks)))
            return {id: self.filenames[id] for item in result for id in item}

    def _convertOneFile(self, filename, fileLocks):
        """Read and process one file, and write its records to the correct
        indexed files, while handling exceptions in a useful way so that they
        don't get swallowed by the multiprocess pool.

        Parameters
        ----------
        filename : `str`
            The file to process.
        fileLocks : `dict` [`int`, `multiprocessing.Lock`]
            A Lock for each HTM pixel; each pixel gets one file written, and
            we need to block when one process is accessing that file.

        Returns
        -------
        pixels, files : `list` [`int`]
            The pixel ids that were written to.
        """
        global FILE_PROGRESS
        inputData = self.file_reader.run(filename)
        fluxes = self._getFluxes(inputData)
        coordErr = self._getCoordErr(inputData)
        matchedPixels = self.indexer.indexPoints(inputData[self.config.ra_name],
                                                 inputData[self.config.dec_name])
        pixel_ids = set(matchedPixels)
        for pixelId in pixel_ids:
            with fileLocks[pixelId]:
                self._doOnePixel(inputData, matchedPixels, pixelId, fluxes, coordErr)
        with FILE_PROGRESS.get_lock():
            oldPercent = 100 * FILE_PROGRESS.value / self.nInputFiles
            FILE_PROGRESS.value += 1
            percent = 100 * FILE_PROGRESS.value / self.nInputFiles
            # only log each "new percent"
            if np.floor(percent) - np.floor(oldPercent) >= 1:
                self.log.info("Completed %d / %d files: %d %% complete ",
                              FILE_PROGRESS.value,
                              self.nInputFiles,
                              percent)
        return pixel_ids

    def _doOnePixel(self, inputData, matchedPixels, pixelId, fluxes, coordErr):
        """Process one HTM pixel, appending to an existing catalog or creating
        a new catalog, as needed.

        Parameters
        ----------
        inputData : `numpy.ndarray`
            The data from one input file.
        matchedPixels : `numpy.ndarray`
            The row-matched pixel indexes corresponding to ``inputData``.
        pixelId : `int`
            The pixel index we are currently processing.
        fluxes : `dict` [`str`, `numpy.ndarray`]
            The values that will go into the flux and fluxErr fields in the
            output catalog.
        coordErr : `dict` [`str`, `numpy.ndarray`]
            The values that will go into the coord_raErr, coord_decErr, and
            coord_ra_dec_Cov fields in the output catalog (in radians).
        """
        idx = np.where(matchedPixels == pixelId)[0]
        catalog = self.getCatalog(pixelId, self.schema, len(idx))
        for outputRow, inputRow in zip(catalog[-len(idx):], inputData[idx]):
            self._fillRecord(outputRow, inputRow)

        global COUNTER
        with COUNTER.get_lock():
            self._setIds(inputData[idx], catalog)

        # set fluxes from the pre-computed array
        for name, array in fluxes.items():
            catalog[self.key_map[name]][-len(idx):] = array[idx]

        # set coordinate errors from the pre-computed array
        for name, array in coordErr.items():
            catalog[name][-len(idx):] = array[idx]

        catalog.writeFits(self.filenames[pixelId])

    def _setIds(self, inputData, catalog):
        """Fill the `id` field of catalog with a running index, filling the
        last values up to the length of ``inputData``.

        Fill with `self.config.id_name` if specified, otherwise use the
        global running counter value.

        Parameters
        ----------
        inputData : `numpy.ndarray`
            The input data that is being processed.
        catalog : `lsst.afw.table.SimpleCatalog`
            The output catalog to fill the ids.
        """
        global COUNTER
        size = len(inputData)
        if self.config.id_name:
            catalog['id'][-size:] = inputData[self.config.id_name]
        else:
            idEnd = COUNTER.value + size
            catalog['id'][-size:] = np.arange(COUNTER.value, idEnd)
            COUNTER.value = idEnd

    def getCatalog(self, pixelId, schema, nNewElements):
        """Get a catalog from disk or create it if it doesn't exist.

        Parameters
        ----------
        pixelId : `dict`
            Identifier for catalog to retrieve
        schema : `lsst.afw.table.Schema`
            Schema to use in catalog creation it does not exist.
        nNewElements : `int`
            The number of new elements that will be added to the catalog,
            so space can be preallocated.

        Returns
        -------
        catalog : `lsst.afw.table.SimpleCatalog`
            The new or read-and-resized catalog specified by `dataId`.
        """
        # This is safe, because we lock on this file before getCatalog is called.
        if os.path.isfile(self.filenames[pixelId]):
            catalog = afwTable.SimpleCatalog.readFits(self.filenames[pixelId])
            catalog.resize(len(catalog) + nNewElements)
            return catalog.copy(deep=True)  # ensure contiguity, so that column-assignment works
        catalog = afwTable.SimpleCatalog(schema)
        catalog.resize(nNewElements)
        self.addRefCatMetadata(catalog)
        return catalog

    @staticmethod
    def computeCoord(row, ra_name, dec_name):
        """Create an ICRS coord. from a row of a catalog being converted.

        Parameters
        ----------
        row : `numpy.ndarray`
            Row from catalog being converted.
        ra_name : `str`
            Name of RA key in catalog being converted.
        dec_name : `str`
            Name of Dec key in catalog being converted.

        Returns
        -------
        coord : `lsst.geom.SpherePoint`
            ICRS coordinate.
        """
        return lsst.geom.SpherePoint(row[ra_name], row[dec_name], lsst.geom.degrees)

    def _getCoordErr(self, inputData, ):
        """Compute the ra/dec error fields that will go into the output catalog.

        Parameters
        ----------
        inputData : `numpy.ndarray`
            The input data to compute fluxes for.

        Returns
        -------
        coordErr : `dict` [`str`, `numpy.ndarray`]
            The values that will go into the coord_raErr, coord_decErr, fields
            in the output catalog (in radians).

        Notes
        -----
        This does not currently handle the ra/dec covariance field,
        ``coord_ra_dec_Cov``. That field may require extra work, as its units
        may be more complicated in external catalogs.
        """
        result = {}
        if hasattr(self, "coord_err_unit"):
            result['coord_raErr'] = u.Quantity(inputData[self.config.ra_err_name],
                                               self.coord_err_unit).to_value(u.radian)
            result['coord_decErr'] = u.Quantity(inputData[self.config.dec_err_name],
                                                self.coord_err_unit).to_value(u.radian)
        return result

    def _setFlags(self, record, row):
        """Set flags in an output record.

        Parameters
        ----------
        record : `lsst.afw.table.SimpleRecord`
            Row from indexed catalog to modify.
        row : `numpy.ndarray`
            Row from catalog being converted.
        """
        names = record.schema.getNames()
        for flag in self._flags:
            if flag in names:
                attr_name = 'is_{}_name'.format(flag)
                record.set(self.key_map[flag], bool(row[getattr(self.config, attr_name)]))

    def _getFluxes(self, inputData):
        """Compute the flux fields that will go into the output catalog.

        Parameters
        ----------
        inputData : `numpy.ndarray`
            The input data to compute fluxes for.

        Returns
        -------
        fluxes : `dict` [`str`, `numpy.ndarray`]
            The values that will go into the flux and fluxErr fields in the
            output catalog.
        """
        result = {}
        for item in self.config.mag_column_list:
            result[item+'_flux'] = (inputData[item]*u.ABmag).to_value(u.nJy)
        if len(self.config.mag_err_column_map) > 0:
            for err_key in self.config.mag_err_column_map.keys():
                error_col_name = self.config.mag_err_column_map[err_key]
                # TODO: multiply by 1e9 here until we have a replacement (see DM-16903)
                # NOTE: copy the arrays because the numpy strides may not be useable by C++.
                fluxErr = fluxErrFromABMagErr(inputData[error_col_name].copy(),
                                              inputData[err_key].copy())*1e9
                result[err_key+'_fluxErr'] = fluxErr
        return result

    def _setProperMotion(self, record, row):
        """Set proper motion fields in a record of an indexed catalog.

        The proper motions are read from the specified columns,
        scaled appropriately, and installed in the appropriate
        columns of the output.

        Parameters
        ----------
        record : `lsst.afw.table.SimpleRecord`
            Row from indexed catalog to modify.
        row : structured `numpy.array`
            Row from catalog being converted.
        """
        if self.config.pm_ra_name is None:  # ConvertReferenceCatalogConfig.validate ensures all or none
            return
        radPerOriginal = np.radians(self.config.pm_scale)/(3600*1000)
        record.set(self.key_map["pm_ra"], row[self.config.pm_ra_name]*radPerOriginal*lsst.geom.radians)
        record.set(self.key_map["pm_dec"], row[self.config.pm_dec_name]*radPerOriginal*lsst.geom.radians)
        record.set(self.key_map["epoch"], self._epochToMjdTai(row[self.config.epoch_name]))
        if self.config.pm_ra_err_name is not None:  # pm_dec_err_name also, by validation
            record.set(self.key_map["pm_raErr"], row[self.config.pm_ra_err_name]*radPerOriginal)
            record.set(self.key_map["pm_decErr"], row[self.config.pm_dec_err_name]*radPerOriginal)

    def _setParallax(self, record, row):
        """Set the parallax fields in a record of a refcat.
        """
        if self.config.parallax_name is None:
            return
        scale = self.config.parallax_scale*lsst.geom.milliarcseconds
        record.set(self.key_map['parallax'], row[self.config.parallax_name]*scale)
        record.set(self.key_map['parallaxErr'], row[self.config.parallax_err_name]*scale)

    def _epochToMjdTai(self, nativeEpoch):
        """Convert an epoch in native format to TAI MJD (a float).
        """
        return astropy.time.Time(nativeEpoch, format=self.config.epoch_format,
                                 scale=self.config.epoch_scale).tai.mjd

    def _setExtra(self, record, row):
        """Set extra data fields in a record of an indexed catalog.

        Parameters
        ----------
        record : `lsst.afw.table.SimpleRecord`
            Row from indexed catalog to modify.
        row : structured `numpy.array`
            Row from catalog being converted.
        """
        for extra_col in self.config.extra_col_names:
            value = row[extra_col]
            # If data read from a text file contains string like entires,
            # numpy stores this as its own internal type, a numpy.str_
            # object. This seems to be a consequence of how numpy stores
            # string like objects in fixed column arrays. This checks
            # if any of the values to be added to the catalog are numpy
            # string types, and if they are, casts them to a python string
            # which is what the python c++ records expect
            if isinstance(value, np.str_):
                value = str(value)
            record.set(self.key_map[extra_col], value)

    def _fillRecord(self, record, row):
        """Fill a record in an indexed catalog to be persisted.

        Parameters
        ----------
        record : `lsst.afw.table.SimpleRecord`
            Row from indexed catalog to modify.
        row : structured `numpy.array`
            Row from catalog being converted.
        """
        record.setCoord(self.computeCoord(row, self.config.ra_name, self.config.dec_name))

        self._setFlags(record, row)
        self._setProperMotion(record, row)
        self._setParallax(record, row)
        self._setExtra(record, row)


class ConvertGaiaManager(ConvertRefcatManager):
    """Special-case convert manager to deal with Gaia fluxes.
    """
    def _getFluxes(self, input):
        result = {}

        def gaiaFluxToFlux(flux, zeroPoint):
            """Equations 5.19 and 5.30 from the Gaia calibration document define the
            conversion from Gaia electron/second fluxes to AB magnitudes.
            https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu5pho/sec_cu5pho_calibr/ssec_cu5pho_calibr_extern.html
            """
            result = ((zeroPoint + -2.5 * np.log10(flux))*u.ABmag).to_value(u.nJy)
            # set 0 instrumental fluxes to 0 (instead of NaN/inf from the math)
            result[flux == 0] = 0
            return result

        # Some fluxes are 0, so log10(flux) can give warnings. We handle the
        # zeros explicitly, so they warnings are irrelevant.
        with np.errstate(invalid='ignore', divide='ignore'):
            # The constants below come from table 5.3 in this document;
            # https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu5pho/sec_cu5pho_calibr/ssec_cu5pho_calibr_extern.html
            result['phot_g_mean_flux'] = gaiaFluxToFlux(input['phot_g_mean_flux'], 25.7934)
            result['phot_bp_mean_flux'] = gaiaFluxToFlux(input['phot_bp_mean_flux'], 25.3806)
            result['phot_rp_mean_flux'] = gaiaFluxToFlux(input['phot_rp_mean_flux'], 25.1161)

        result['phot_g_mean_fluxErr'] = result['phot_g_mean_flux'] / input['phot_g_mean_flux_over_error']
        result['phot_bp_mean_fluxErr'] = result['phot_bp_mean_flux'] / input['phot_bp_mean_flux_over_error']
        result['phot_rp_mean_fluxErr'] = result['phot_rp_mean_flux'] / input['phot_rp_mean_flux_over_error']

        return result
