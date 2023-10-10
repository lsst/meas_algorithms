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

__all__ = ["ReadFitsCatalogConfig", "ReadFitsCatalogTask"]

from astropy.table import Table

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase


class ReadFitsCatalogConfig(pexConfig.Config):
    hdu = pexConfig.Field(
        dtype=int,
        default=1,
        doc="HDU containing the desired binary table, 0-based but a binary table never occurs in HDU 0",
    )
    column_map = pexConfig.DictField(
        doc="Mapping of input column name: output column name; each specified column must exist, "
            "but additional columns in the input data are written using their original name. ",
        keytype=str,
        itemtype=str,
        default={},
    )


class ReadFitsCatalogTask(pipeBase.Task):
    """Read an object catalog from a FITS table

    The resulting `catalogArray` is a numpy structured array containing fields such as "name", "ra" and "dec"
    and a few rows of data. For more complicated cases config parameters allow you to rename columns
    and choose which HDU to read.
    """
    _DefaultName = 'readCatalog'
    ConfigClass = ReadFitsCatalogConfig

    def run(self, filename):
        """Read an object catalog from the specified FITS file

        Parameters
        ----------
        filename : `str`
            Path to specified FITS file

        Returns
        -------
        table : `np.array`
            a numpy structured array containing the specified columns
        """

        # Set character_as_bytes=False to ensure that all string columns are
        # set to python string types rather than byte arrays.
        table = Table.read(filename, hdu=self.config.hdu, character_as_bytes=False)
        if table is None:
            raise RuntimeError("No data found in %s HDU %s" % (filename, self.config.hdu))

        if not self.config.column_map:
            # take the data as it is
            return table.as_array()

        missingnames = set(self.config.column_map.keys()) - set(table.columns.keys())
        if missingnames:
            raise RuntimeError("Columns %s in column_map were not found in %s" % (missingnames, filename))

        for inname, outname in self.config.column_map.items():
            table.columns[inname].name = outname
        return table.as_array()
