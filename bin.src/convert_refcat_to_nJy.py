#!/usr/bin/env python
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

"""Convert old HTM reference catalogs to use nJy for fluxes, instead of Jy.

This will process all .fits files in the given directory, converting the units
of any `_flux`, `_fluxSigma`, and `_fluxErr` fields from Jy->nJy and making
their units field in the schema `nJy`, then overwriting the files with the
updated values. Also replaces `_fluxSigma` with `_fluxErr` in the schema,
per RFC-333. If all flux fields in the refcat schema have units of `'nJy'`,
the files are not modified.

Many of our old reference catalogs have no units for their fluxes: we assume
(as the algorithmic code did) that these are all in Jy units.

By default, only print the files and fields that will be modified, but do not
write the output.

If you are processing a large number of files (e.g. ps1_pv3), we recommend
capturing stdout to a log file, and using the -n8 option to parallelize it.
"""
import os.path
import glob

import concurrent.futures
import itertools

import lsst.afw.table
from lsst.meas.algorithms import DatasetConfig
from lsst.meas.algorithms.loadReferenceObjects import convertToNanojansky, hasNanojanskyFluxUnits
from lsst.meas.algorithms.ingestIndexReferenceTask import addRefCatMetadata
import lsst.log


def is_old_schema(config, filename):
    """Check whether this file's schema has "old-style" fluxes."""
    catalog = lsst.afw.table.SimpleCatalog.readFits(filename)
    return (config.format_version == 0) and (not hasNanojanskyFluxUnits(catalog.schema))


def process_one(filename, write=False, quiet=False):
    """Convert one file in-place from Jy (or no units) to nJy fluxes.

    Parameters
    ----------
    filename : `str`
        The file to convert.
    write : `bool`, optional
        Write the converted catalog out, overwriting the read in catalog?
    quiet : `bool`, optional
        Do not print messages about files read/written or fields found?
    """
    log = lsst.log.Log()
    if quiet:
        log.setLevel(lsst.log.WARN)

    log.info("Reading: %s", filename)
    catalog = lsst.afw.table.SimpleCatalog.readFits(filename)

    output = convertToNanojansky(catalog, log, doConvert=write)

    if write:
        addRefCatMetadata(output)
        output.writeFits(filename)
        log.info("Wrote: %s", filename)


def main():
    import argparse
    import sys

    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=CustomFormatter)
    parser.add_argument("path",
                        help="Directory (written by IngestIndexedReferenceTask) containing the"
                        " reference catalogs to overwrite."
                        " All files with a `.fits` extension in the directory will be processed,"
                        " including `master_schema.fits`, which must exist.")
    parser.add_argument('-n', '--nprocesses', default=1, type=int,
                        help="Number of processes to use when reading and writing files.")
    parser.add_argument('--write', action="store_true",
                        help="Write the corrected files (default just prints what would have changed).")
    parser.add_argument('--quiet', action="store_true",
                        help="Be less verbose about what files and fields are being converted.")
    args = parser.parse_args()

    schema_file = os.path.join(args.path, "master_schema.fits")
    if not os.path.isfile(schema_file):
        print("Error: Cannot find master_schema.fits in supplied path:", args.path)
        sys.exit(-1)
    configPath = os.path.join(args.path, 'config.py')
    config = DatasetConfig()
    config.load(configPath)
    if not is_old_schema(config, schema_file):
        print("Catalog does not contain old-style fluxes; nothing to convert.")
        sys.exit(0)

    files = glob.glob(os.path.join(args.path, "*.fits"))
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.nprocesses) as executor:
        futures = executor.map(process_one, files, itertools.repeat(args.write), itertools.repeat(args.quiet))
        # we have to at least loop over the futures, otherwise exceptions will be lost
        for future in futures:
            pass

    if args.write:
        config.format_version = 1
        # update the docstring to annotate the file.
        msg = "\nUpdated refcat from version 0->1 to have nJy flux units via convert_refcat_to_nJy.py"
        config._fields['format_version'].doc += msg
        config.save(configPath)
        if not args.quiet:
            print("Added `format_version=1` to config.py")


if __name__ == "__main__":
    main()
