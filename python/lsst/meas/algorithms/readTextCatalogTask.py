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
    fill_values = pexConfig.ListField(
        dtype=str,
        default=None,
        optional=True,
        doc=("A list giving [<match_string>, <fill_value>], which is used to mask"
             " the given values in the input file. '0' is suggested for the fill value in order to prevent"
             " changing the column datatype. The default behavior is to fill empty data with zeros. See "
             "https://docs.astropy.org/en/stable/io/ascii/read.html#bad-or-missing-values for more details."
             "Use `replace_missing_floats_with_nan` to change floats to NaN instead of <fill_value>.")
    )
    replace_missing_floats_with_nan = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="If True, replace missing data in float columns with NaN instead of zero. If `fill_values` is "
            "set, this parameter with replace the floats identified as missing by `fill_values`, and the fill"
            " value from `fill_values` will be overridden with NaN for floats."
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

        if self.config.fill_values:
            kwargs['fill_values'] = [list(self.config.fill_values)]

        table = Table.read(filename, format=self.config.format,
                           delimiter=self.config.delimiter,
                           **kwargs)

        # convert to a numpy array for backwards compatibility with other readers
        arr = np.array(table.as_array())

        if self.config.replace_missing_floats_with_nan:
            for column in table.columns:
                if (table.dtype[column] == np.float32) or (table.dtype[column] == np.float64):
                    arr[column][table.mask[column]] = np.nan

        return arr
