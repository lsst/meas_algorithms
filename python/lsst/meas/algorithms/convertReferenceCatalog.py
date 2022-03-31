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

"""
Convert an external reference catalog into the hierarchical triangular mesh
(HTM) sharded LSST-style format, to be ingested into the butler.
"""

__all__ = ["ConvertReferenceCatalogTask"]

import argparse
import glob
import os
import pathlib
import logging

import astropy

from . import ConvertReferenceCatalogBase


class ConvertReferenceCatalogTask(ConvertReferenceCatalogBase):
    """Class for producing HTM-indexed reference catalogs from external
    catalog data.

    Parameters
    ----------
    output_dir : `str`
        The path to write the output files to, in a subdirectory defined by
        ``DatasetConfig.ref_dataset_name``.
    """
    _DefaultName = 'ConvertReferenceCatalogTask'

    def __init__(self, *, output_dir=None, **kwargs):
        super().__init__(**kwargs)
        if output_dir is None:
            raise RuntimeError("Must specify output_dir.")
        self.base_dir = output_dir
        self.output_dir = os.path.join(output_dir, self.config.dataset_config.ref_dataset_name)
        self.ingest_table_file = os.path.join(self.base_dir, "filename_to_htm.ecsv")

    def _preRun(self):
        # Create the output path, if it doesn't exist; fail if the path exists:
        # we don't want to accidentally append to existing files.
        pathlib.Path(self.output_dir).mkdir(exist_ok=False)

    def _postRun(self, result):
        # Write the astropy table containing the htm->filename relationship
        dimension = f"htm{self.config.dataset_config.indexer.active.depth}"
        table = astropy.table.Table(names=("filename", dimension), dtype=('str', 'int'))
        for key in result:
            table.add_row((result[key], key))
        table.write(self.ingest_table_file)

    def _persistConfig(self):
        filename = os.path.join(self.output_dir, "config.py")
        with open(filename, 'w') as file:
            self.config.dataset_config.saveToStream(file)

    def _getOnePixelFilename(self, start):
        return os.path.join(self.output_dir, f"{self.indexer.htm}.fits")

    def _writeMasterSchema(self, catalog):
        filename = os.path.join(self.output_dir, "master_schema.fits")
        catalog.writeFits(filename)

    def _reduce_kwargs(self):
        # Need to be able to pickle this class to use the multiprocess manager.
        kwargs = super()._reduce_kwargs()
        kwargs['output_dir'] = self.base_dir
        return kwargs


def build_argparser():
    """Construct an argument parser for the ``convertReferenceCatalog`` script.

    Returns
    -------
    argparser : `argparse.ArgumentParser`
        The argument parser that defines the ``convertReferenceCatalog``
        command-line interface.
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='More information is available at https://pipelines.lsst.io.'
    )
    parser.add_argument("outputDir",
                        help="Path to write the output shard files, configs, and `ingest-files` table to.")
    parser.add_argument("configFile",
                        help="File containing the ConvertReferenceCatalogConfig fields.")
    # Use a "+"-list here, so we can produce a more useful error if the user
    # uses an unquoted glob that gets shell expanded.
    parser.add_argument("fileglob", nargs="+",
                        help="Quoted glob for the files to be read in and converted."
                             " Example (note required quotes to prevent shell expansion):"
                             ' "gaia_source/csv/GaiaSource*"')
    return parser


def run_convert(outputDir, configFile, fileglob):
    """Run `ConvertReferenceCatalogTask` on the input arguments.

    Parameters
    ----------
    outputDir : `str`
        Path to write the output files to.
    configFile : `str`
        File specifying the ``ConvertReferenceCatalogConfig`` fields.
    fileglob : `str`
        Quoted glob for the files to be read in and converted.
    """
    # We have to initialize the logger manually when running from the commandline.
    logging.basicConfig(level=logging.INFO, format="{name} {levelname}: {message}", style="{")

    config = ConvertReferenceCatalogTask.ConfigClass()
    config.load(configFile)
    config.validate()
    converter = ConvertReferenceCatalogTask(output_dir=outputDir, config=config)
    files = glob.glob(fileglob)
    converter.run(files)
    with open(os.path.join(outputDir, "convertReferenceCatalogConfig.py"), "w") as outfile:
        converter.config.saveToStream(outfile)
    msg = ("Completed refcat conversion.\n\n"
           "Ingest the resulting files with the following commands, substituting the path\n"
           "to your butler repo for `REPO`, and the ticket number you are tracking this\n"
           "ingest on for `DM-NNNNN`:\n"
           f"\n    butler register-dataset-type REPO {config.dataset_config.ref_dataset_name} "
           "SimpleCatalog htm7"
           "\n    butler ingest-files -t direct REPO gaia_dr2 refcats/DM-NNNNN "
           f"{converter.ingest_table_file}"
           "\n    butler collection-chain REPO --mode extend refcats refcats/DM-NNNNN")
    print(msg)


def main():
    args = build_argparser().parse_args()
    if len(args.fileglob) > 1:
        raise RuntimeError("Final argument must be a quoted file glob, not a shell-expanded list of files.")
    # Fileglob comes out as a length=1 list, so we can test it above.
    run_convert(args.outputDir, args.configFile, args.fileglob[0])
