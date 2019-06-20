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

## @addtogroup LSST_task_documentation
## @{
## @page ReadTextCatalogTask
## @ref ReadTextCatalogTask_ "ReadTextCatalogTask"
## @copybrief ReadTextCatalogTask
## @}


class ReadTextCatalogTask(pipeBase.Task):
    r"""!Read an object catalog from a text file

    @anchor ReadTextCatalogTask_

    @section meas_algorithms_readTextCatalog_Contents  Contents

     - @ref meas_algorithms_readTextCatalog_Purpose
     - @ref meas_algorithms_readTextCatalog_Initialize
     - @ref meas_algorithms_readTextCatalog_Config
     - @ref meas_algorithms_readTextCatalog_Example

    @section meas_algorithms_readTextCatalog_Purpose  Description

    Read an object catalog from a text file. Designed to read foreign catalogs
    so they can be written out in a form suitable for IngestIndexedReferenceTask.

    The file is assumed to be encoded as UTF-8 (which is compatible with ASCII).

    @section meas_algorithms_readTextCatalog_Initialize  Task initialisation

    @copydoc \_\_init\_\_

    @section meas_algorithms_readTextCatalog_Config  Configuration parameters

    See @ref ReadTextCatalogConfig

    @section meas_algorithms_readTextCatalog_Example   A complete example of using ReadTextCatalogTask

    Given a file named `table.csv` containing the following:

        ra,     dec,    flux
        5.5,    -45.2,  12453
        19.6,   34.2,   32123

    you can read this file with the following code:

        from lsst.meas.algorithms.readTextCatalogTask import ReadTextCatalogTask
        task = ReadTextCatalogTask()
        catalogArray = task.run("table.csv")

    The resulting `catalogArray` is a numpy structured array containing three fields
    ("ra", "dec" and "flux") and two rows of data. For more complex cases,
    config parameters allow you to specify the names of the columns (instead of using automatic discovery)
    and set the number of rows to skip.
    """
    _DefaultName = 'readCatalog'
    ConfigClass = ReadTextCatalogConfig

    def run(self, filename):
        """Read an object catalog from the specified text file

        @param[in] filename  path to text file
        @return a numpy structured array containing the specified columns
        """
        kwargs = {}
        if self.config.colnames:
            kwargs['names'] = self.config.colnames
            # if we specify the column names, then we need to just ignore the header lines.
            kwargs['data_start'] = self.config.header_lines
        else:
            # if we don't specify column names, start the header at this line.
            kwargs['header_start'] = self.config.header_lines

        # return a numpy array for backwards compatibility with other readers
        return np.array(Table.read(filename, format=self.config.format,
                                   delimiter=self.config.delimiter,
                                   **kwargs).as_array())
