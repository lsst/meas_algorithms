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

from convertReferenceCatalogTestBase import makeConvertConfig


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
        outdir = tempfile.TemporaryDirectory()
        outpath = outdir.name
        args = ["convertReferenceCatalog",
                outpath,
                os.path.join(self.inpath, "mock_config.py"),
                os.path.join(self.inpath, "*.fits")]
        with unittest.mock.patch.object(convertReferenceCatalog.ConvertReferenceCatalogTask, "run") as run, \
                unittest.mock.patch.object(sys, "argv", args):
            convertReferenceCatalog.main()
            # Test with sets because the glob can come out in any order.
            self.assertEqual(set(run.call_args.args[0]), set(self.expected_files))
        # This is necessary to avoid a ResourceWarning.
        outdir.cleanup()

    def test_main_args_bad_config(self):
        """Test that a bad config file produces a useful error, i.e. that
        main() validates the config.
        """
        outdir = tempfile.TemporaryDirectory()
        outpath = outdir.name
        args = ["convertReferenceCatalog",
                outpath,
                os.path.join(self.inpath, "bad_config.py"),
                os.path.join(self.inpath, "*.fits")]
        with self.assertRaisesRegex(FieldValidationError, "Field 'ra_name' failed validation"), \
                unittest.mock.patch.object(sys, "argv", args):
            convertReferenceCatalog.main()
        # This is necessary to avoid a ResourceWarning.
        outdir.cleanup()

    def test_main_args_expanded_glob(self):
        """Test that an un-quoted glob (i.e. list of files) fails with a
        useful error.
        """
        outdir = tempfile.TemporaryDirectory()
        outpath = outdir.name
        args = ["convertReferenceCatalog",
                outpath,
                os.path.join(self.inpath, "mock_config.py"),
                # an un-quoted glob will be shell-expanded to a list of files.
                "file1", "file2", "file3"]
        msg = "Final argument must be a quoted file glob, not a shell-expanded list of files."
        with self.assertRaisesRegex(RuntimeError, msg), \
                unittest.mock.patch.object(sys, "argv", args):
            convertReferenceCatalog.main()
        # This is necessary to avoid a ResourceWarning.
        outdir.cleanup()


class ConvertReferenceCatalogConfigValidateTestCase(lsst.utils.tests.TestCase):
    """Test validation of ConvertReferenceCatalogConfig."""
    def testValidateRaDecMag(self):
        config = makeConvertConfig()
        config.validate()

        for name in ("ra_name", "dec_name", "mag_column_list"):
            with self.subTest(name=name):
                config = makeConvertConfig()
                setattr(config, name, None)
                with self.assertRaises(ValueError):
                    config.validate()

    def testValidateRaDecErr(self):
        # check that a basic config validates
        config = makeConvertConfig(withRaDecErr=True)
        config.validate()

        # check that a config with any of these fields missing does not validate
        for name in ("ra_err_name", "dec_err_name", "coord_err_unit"):
            with self.subTest(name=name):
                config = makeConvertConfig(withRaDecErr=True)
                setattr(config, name, None)
                with self.assertRaises(ValueError):
                    config.validate()

        # check that coord_err_unit must be an astropy unit
        config = makeConvertConfig(withRaDecErr=True)
        config.coord_err_unit = "nonsense unit"
        with self.assertRaisesRegex(ValueError, "is not a valid astropy unit string"):
            config.validate()

    def testValidateMagErr(self):
        config = makeConvertConfig(withMagErr=True)
        config.validate()

        # test for missing names
        for name in config.mag_column_list:
            with self.subTest(name=name):
                config = makeConvertConfig(withMagErr=True)
                del config.mag_err_column_map[name]
                with self.assertRaises(ValueError):
                    config.validate()

        # test for incorrect names
        for name in config.mag_column_list:
            with self.subTest(name=name):
                config = makeConvertConfig(withMagErr=True)
                config.mag_err_column_map["badName"] = config.mag_err_column_map[name]
                del config.mag_err_column_map[name]
                with self.assertRaises(ValueError):
                    config.validate()

    def testValidatePm(self):
        basicNames = ["pm_ra_name", "pm_dec_name", "epoch_name", "epoch_format", "epoch_scale"]

        for withPmErr in (False, True):
            config = makeConvertConfig(withPm=True, withPmErr=withPmErr)
            config.validate()
            del config

            if withPmErr:
                names = basicNames + ["pm_ra_err_name", "pm_dec_err_name"]
            else:
                names = basicNames
                for name in names:
                    with self.subTest(name=name, withPmErr=withPmErr):
                        config = makeConvertConfig(withPm=True, withPmErr=withPmErr)
                        setattr(config, name, None)
                        with self.assertRaises(ValueError):
                            config.validate()

    def testValidateParallax(self):
        """Validation should fail if any parallax-related fields are missing.
        """
        names = ["parallax_name", "epoch_name", "epoch_format", "epoch_scale", "parallax_err_name"]

        config = makeConvertConfig(withParallax=True)
        config.validate()
        del config

        for name in names:
            with self.subTest(name=name):
                config = makeConvertConfig(withParallax=True)
                setattr(config, name, None)
                with self.assertRaises(ValueError, msg=name):
                    config.validate()


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
