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

import os.path
import sys
import unittest
import unittest.mock
import tempfile

from lsst.pex.config import FieldValidationError
from lsst.meas.algorithms import convertReferenceCatalog
import lsst.utils


class TestMain(lsst.utils.tests.TestCase):
    """Test mocking commandline arguments and calling
    ``convertReferenceCatalog.main()``.
    """
    def setUp(self):
        self.inpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/mockrefcat/")
        self.expected_files = [os.path.join(self.inpath, "123.fits"),
                               os.path.join(self.inpath, "124.fits"),
                               os.path.join(self.inpath, "125.fits")]

    def test_main_args(self):
        """Test that main configures the task and calls run() with the correct
        file list.
        """
        outpath = tempfile.mkdtemp()
        args = ["convertReferenceCatalog",
                outpath,
                os.path.join(self.inpath, "mock_config.py"),
                os.path.join(self.inpath, "*.fits")]
        with unittest.mock.patch.object(convertReferenceCatalog.ConvertReferenceCatalogTask, "run") as run, \
                unittest.mock.patch.object(sys, "argv", args):
            convertReferenceCatalog.main()
            # Test with sets because the glob can come out in any order.
            self.assertEqual(set(run.call_args.args[0]), set(self.expected_files))

    def test_main_args_bad_config(self):
        """Test that a bad config file produces a useful error, i.e. that
        main() validates the config.
        """
        outpath = tempfile.mkdtemp()
        args = ["convertReferenceCatalog",
                outpath,
                os.path.join(self.inpath, "bad_config.py"),
                os.path.join(self.inpath, "*.fits")]
        with self.assertRaisesRegex(FieldValidationError, "Field 'ra_name' failed validation"), \
                unittest.mock.patch.object(sys, "argv", args):
            convertReferenceCatalog.main()

    def test_main_args_expanded_glob(self):
        """Test that an un-quoted glob (i.e. list of files) fails with a
        useful error.
        """
        outpath = tempfile.mkdtemp()
        args = ["convertReferenceCatalog",
                outpath,
                os.path.join(self.inpath, "mock_config.py"),
                # an un-quoted glob will be shell-expanded to a list of files.
                "file1", "file2", "file3"]
        msg = "Final argument must be a quoted file glob, not a shell-expanded list of files."
        with self.assertRaisesRegex(RuntimeError, msg), \
                unittest.mock.patch.object(sys, "argv", args):
            convertReferenceCatalog.main()


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
