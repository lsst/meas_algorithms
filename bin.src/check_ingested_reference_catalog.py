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

"""Test that an ingested reference catalog matches some subset of the original
input data.

Example usage (note both sets of quotes around the glob and resulting envvar!):
    FILEGLOB='/project/shared/data/gaia_dr2/gaia_source/csv/*.csv.gz'
    check_ingested_reference_catalog.py --config gaia_dr2_config.py --ref_name gaia-dr2 refcat/ "$FILEGLOB"
"""

import glob
import os.path
import random

import lsst.daf.persistence
import lsst.geom
from lsst.meas.algorithms import LoadIndexedReferenceObjectsConfig, LoadIndexedReferenceObjectsTask
from lsst.meas.algorithms.ingestIndexReferenceTask import IngestIndexedReferenceConfig


def do_one_file(filename, refObjLoader, reader, config, n):
    """
    Parameters
    ----------
    filename : `str`
        The file to read data from.
    refObjLoader : `lsst.meas.algorithms.LoadReferenceObjectsTask`
        Reference loader to use to search for the sources from ``filename``.
    reader : `lsst.meas.algorithms.ReadTextCatalogTask`
        File reader to use to load the data in ``filename``.
    config : `lsst.pex.config.Config`
        Configuration used to originally ingest the reference catalog, used
        to identify the mappings between input and output column names.
    n : `int`
        Number of sources from the file to check (randomly sampled).
    """
    data = reader.run(filename)
    filterName = config.mag_column_list[0]
    successCount = 0
    print("Checking:", os.path.basename(filename), end=" : ")
    for i in random.sample(range(len(data)), n):
        center = lsst.geom.SpherePoint(data[i][config.ra_name], data[i][config.dec_name], lsst.geom.degrees)
        radius = 1*lsst.geom.arcseconds  # don't need to search far to get an exact match!
        refCat = refObjLoader.loadSkyCircle(center, radius, filterName).refCat
        if any(refCat['id'] == data[i][config.id_name]):
            successCount += 1
    print(f"{successCount} / {n}")


def main():
    import argparse

    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=CustomFormatter)
    parser.add_argument("refCatPath",
                        help="Butler repo (written by IngestIndexedReferenceTask) containing the"
                        " reference catalogs to test.")
    parser.add_argument("inputGlob",
                        help="A glob pattern specifying the files (read by IngestIndexedReferenceTask) to "
                        " check against the refCatPath output. (e.g. '/datasets/foo/*.csv'")
    parser.add_argument("--nFiles", default=5, type=int,
                        help="Number of input files to test (randomly selected from inputGlob).")
    parser.add_argument("--nPerFile", default=100, type=int,
                        help="Number of objects to test per file (randomly selected).")
    parser.add_argument("--config",
                        help="A IngestIndexedReferenceConfig config file, for the field name mappings.")
    parser.add_argument("--ref_name",
                        help="The name of the reference catalog stored in refCatPath.")
    args = parser.parse_args()

    ingestConfig = IngestIndexedReferenceConfig()
    ingestConfig.load(args.config)

    butler = lsst.daf.persistence.Butler(args.refCatPath)
    refObjConfig = LoadIndexedReferenceObjectsConfig()
    refObjConfig.ref_dataset_name = args.ref_name
    refObjLoader = LoadIndexedReferenceObjectsTask(butler, config=refObjConfig)
    reader = ingestConfig.file_reader.target()

    files = glob.glob(args.inputGlob)
    files.sort()
    for filename in random.sample(files, args.nFiles):
        do_one_file(filename, refObjLoader, reader, ingestConfig, args.nPerFile)


if __name__ == "__main__":
    main()
