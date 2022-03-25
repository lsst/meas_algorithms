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

__all__ = ["ReadTextCatalogConfig", "ReadTextCatalogTask"]

import numpy as np
from astropy.table import Table

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase


class ReadTextCatalogConfig(pexConfig.Config):
    header_lines = pexConfig.Field(
        dtype=int,
        default=0,
        doc='Number of lines to skip when reading the text reference file.'
    )
    colnames = pexConfig.ListField(
        dtype=str,
        default=[],
        doc="An ordered list of column names to use in ingesting the catalog. "
            "With an empty list, column names will be discovered from the first line "
            "after the skipped header lines."
    )
    delimiter = pexConfig.Field(
        dtype=str,
        default=',',
        doc='Delimiter to use when reading text reference files.  Comma is default.'
    )
    format = pexConfig.Field(
        dtype=str,
        default='csv',
        doc=("Format of files to read, from the astropy.table I/O list here:"
             "http://docs.astropy.org/en/stable/io/unified.html#built-in-table-readers-writers")
    )


class ReadTextCatalogTask(pipeBase.Task):
    """Read an object catalog from a text file
    """
    _DefaultName = 'readCatalog'
    ConfigClass = ReadTextCatalogConfig

    def run(self, filename):
        """Read an object catalog from the specified text file

        Parameters
        ----------
        filename : `string`
            Path to specified text file

        Returns
        -------
        A numpy structured array containing the specified columns
        """
        kwargs = {}
        if self.config.colnames:
            # Wrap in list() to avoid transferring a pex_config proxy object.
            kwargs['names'] = list(self.config.colnames)
            # if we specify the column names, then we need to just ignore the header lines.
            kwargs['data_start'] = self.config.header_lines
        else:
            # if we don't specify column names, start the header at this line.
            kwargs['header_start'] = self.config.header_lines

        # return a numpy array for backwards compatibility with other readers
        return np.array(Table.read(filename, format=self.config.format,
                                   delimiter=self.config.delimiter,
                                   **kwargs).as_array())
