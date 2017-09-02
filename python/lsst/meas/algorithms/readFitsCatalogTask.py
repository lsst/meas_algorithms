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

__all__ = ["ReadFitsCatalogConfig", "ReadFitsCatalogTask"]

from astropy.io import fits

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

## \addtogroup LSST_task_documentation
## \{
## \page ReadFitsCatalogTask
## \ref ReadFitsCatalogTask_ "ReadFitsCatalogTask"
## \copybrief ReadFitsCatalogTask
## \}


class ReadFitsCatalogTask(pipeBase.Task):
    """!Read an object catalog from a FITS table

    @anchor ReadFitsCatalogTask_

    @section meas_algorithms_readFitsCatalog_Contents  Contents

     - @ref meas_algorithms_readFitsCatalog_Purpose
     - @ref meas_algorithms_readFitsCatalog_Initialize
     - @ref meas_algorithms_readFitsCatalog_Config
     - @ref meas_algorithms_readFitsCatalog_Example

    @section meas_algorithms_readFitsCatalog_Purpose  Description

    Read an object catalog from a FITS table. Designed to read foreign catalogs
    so they can be written out in a form suitable for IngestIndexedReferenceTask.

    @section meas_algorithms_readFitsCatalog_Initialize  Task initialisation

    @copydoc \_\_init\_\_

    @section meas_algorithms_readFitsCatalog_Config  Configuration parameters

    See @ref ReadFitsCatalogConfig

    @section meas_algorithms_readFitsCatalog_Example   A complete example of using ReadFitsCatalogTask

    Run the following code from the main directory of meas_algorithms:

        from lsst.meas.algorithms.readFitsCatalogTask import ReadFitsCatalogTask
        filePath = "tests/data/testReadFitsCatalog.fits"
        task = ReadFitsCatalogTask()
        catalogArray = task.run(filePath)

    The resulting `catalogArray` is a numpy structured array containing fields such as "name", "ra" and "dec"
    and a few rows of data. For more complicated cases config parameters allow you to rename columns
    and choose which HDU to read.
    """
    _DefaultName = 'readCatalog'
    ConfigClass = ReadFitsCatalogConfig

    def run(self, filename):
        """Read an object catalog from the specified FITS file

        @param[in] filename  path to FITS file
        @return a numpy structured array containing the specified columns
        """
        with fits.open(filename) as f:
            hdu = f[self.config.hdu]
            if hdu.data is None:
                raise RuntimeError("No data found in %s HDU %s" % (filename, self.config.hdu))
            if hdu.is_image:
                raise RuntimeError("%s HDU %s is an image" % (filename, self.config.hdu))

            if not self.config.column_map:
                # take the data as it is
                return hdu.data

            missingnames = set(self.config.column_map.keys()) - set(hdu.columns.names)
            if missingnames:
                raise RuntimeError("Columns %s in column_map were not found in %s" % (missingnames, filename))

            for inname, outname in self.config.column_map.items():
                hdu.columns[inname].name = outname
            return hdu.data
