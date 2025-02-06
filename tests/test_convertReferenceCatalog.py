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

import astropy.units as u
import numpy as np
import os.path
import sys
import unittest
import unittest.mock
import tempfile
import itertools
import logging

from lsst.afw.table import SimpleCatalog
from lsst.pex.config import FieldValidationError
from lsst.meas.algorithms import (convertReferenceCatalog, ConvertReferenceCatalogTask, getRefFluxField)
from lsst.meas.algorithms.readTextCatalogTask import ReadTextCatalogTask
from lsst.meas.algorithms.htmIndexer import HtmIndexer
from lsst.meas.algorithms.convertRefcatManager import ConvertGaiaManager
from lsst.meas.algorithms.convertReferenceCatalog import addRefCatMetadata, _makeSchema

import lsst.utils

from convertReferenceCatalogTestBase import makeConvertConfig
import convertReferenceCatalogTestBase


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


class MakeSchemaTestCase(lsst.utils.tests.TestCase):
    """Test the function to make reference catalog schemas.
    """
    def testMakeSchema(self):
        """Make a schema and check it."""
        for filterNameList in (["r"], ["foo", "_bar"]):
            for (addIsPhotometric, addIsResolved, addIsVariable) in itertools.product((False, True),
                                                                                      (False, True),
                                                                                      (False, True)):
                argDict = dict(
                    filterNameList=filterNameList,
                    addIsPhotometric=addIsPhotometric,
                    addIsResolved=addIsResolved,
                    addIsVariable=addIsVariable,
                )
                refSchema = _makeSchema(**argDict)
                self.assertTrue("coord_ra" in refSchema)
                self.assertTrue("coord_dec" in refSchema)
                self.assertTrue("coord_raErr" in refSchema)
                self.assertTrue("coord_decErr" in refSchema)
                for filterName in filterNameList:
                    fluxField = filterName + "_flux"
                    self.assertIn(fluxField, refSchema)
                    self.assertNotIn("x" + fluxField, refSchema)
                    fluxErrField = fluxField + "Err"
                    self.assertIn(fluxErrField, refSchema)
                    self.assertEqual(getRefFluxField(refSchema, filterName), filterName + "_flux")
                self.assertEqual("resolved" in refSchema, addIsResolved)
                self.assertEqual("variable" in refSchema, addIsVariable)
                self.assertEqual("photometric" in refSchema, addIsPhotometric)
                self.assertEqual("photometric" in refSchema, addIsPhotometric)

                # The default for `fullPositionInformation` is False, so none
                # of the following should be included. We test setting these
                # all together below.
                self.assertNotIn("epoch", refSchema)
                self.assertNotIn("pm_ra", refSchema)
                self.assertNotIn("pm_dec", refSchema)
                self.assertNotIn("pm_flag", refSchema)
                self.assertNotIn("parallax", refSchema)
                self.assertNotIn("parallax_flag", refSchema)

    def testMakeSchema_fullCovariance(self):
        """Make a schema with full position information and coordinate
        covariance and test it."""
        refSchema = _makeSchema(filterNameList=["r"], fullPositionInformation=True)
        # Test that the epoch, proper motion and parallax terms are included in
        # the schema.
        self.assertIn("epoch", refSchema)
        self.assertIn("pm_ra", refSchema)
        self.assertIn("pm_dec", refSchema)
        self.assertIn("pm_flag", refSchema)
        self.assertIn("parallax", refSchema)
        self.assertIn("parallax_flag", refSchema)
        # Test that a sample of the 15 covariance terms are included in the schema.
        self.assertIn("coord_raErr", refSchema)
        self.assertIn("coord_decErr", refSchema)
        self.assertIn("coord_ra_coord_dec_Cov", refSchema)
        self.assertIn("pm_raErr", refSchema)
        self.assertIn("pm_ra_parallax_Cov", refSchema)
        self.assertIn("parallaxErr", refSchema)
        self.assertEqual(refSchema['coord_raErr'].asField().getUnits(), "rad")
        self.assertEqual(refSchema['coord_ra_coord_dec_Cov'].asField().getUnits(), "rad^2")
        self.assertEqual(refSchema['pm_raErr'].asField().getUnits(), "rad/year")
        self.assertEqual(refSchema['pm_dec_parallax_Cov'].asField().getUnits(), "rad^2/year")


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
        names = ["pm_ra_name", "pm_dec_name", "epoch_name", "epoch_format", "epoch_scale",
                 "pm_ra_err_name", "pm_dec_err_name"]

        config = makeConvertConfig(withPm=True)
        config.validate()
        del config

        for name in names:
            with self.subTest(name=name):
                config = makeConvertConfig(withPm=True)
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

    def testValidateCovariance(self):
        """Validation should fail if any position-related fields are empty if
        full_position_information is set.
        """
        names = ["ra_err_name", "dec_err_name", "coord_err_unit",
                 "parallax_name", "parallax_err_name",
                 "epoch_name", "epoch_format", "epoch_scale",
                 "pm_ra_name", "pm_dec_name", "pm_ra_err_name", "pm_dec_err_name"]

        for name in names:
            with self.subTest(name=name):
                config = makeConvertConfig(withRaDecErr=True, withParallax=True, withPm=True)
                config.full_position_information = True
                config.manager.retarget(ConvertGaiaManager)
                setattr(config, name, None)
                with self.assertRaises(ValueError, msg=name):
                    config.validate()


class ConvertGaiaManagerTests(convertReferenceCatalogTestBase.ConvertReferenceCatalogTestBase,
                              lsst.utils.tests.TestCase):
    """Unittests specific to the Gaia catalog.
    """
    def setUp(self):
        self.tempDir = tempfile.TemporaryDirectory()
        tempPath = self.tempDir.name
        self.log = logging.getLogger("lsst.TestConvertRefcatManager")
        self.config = convertReferenceCatalogTestBase.makeConvertConfig(withRaDecErr=True)
        self.config.id_name = 'id'
        self.config.full_position_information = True
        self.config.manager.retarget(ConvertGaiaManager)
        self.config.coord_err_unit = 'milliarcsecond'
        self.config.ra_err_name = 'ra_error'
        self.config.dec_err_name = 'dec_error'
        self.config.pm_ra_name = 'pmra'
        self.config.pm_dec_name = 'pmdec'
        self.config.pm_ra_err_name = 'pmra_error'
        self.config.pm_dec_err_name = 'pmdec_error'
        self.config.parallax_name = 'parallax'
        self.config.parallax_err_name = 'parallax_error'
        self.config.epoch_name = 'unixtime'
        self.config.epoch_format = 'unix'
        self.config.epoch_scale = 'tai'
        self.depth = 2  # very small depth, for as few pixels as possible.
        self.indexer = HtmIndexer(self.depth)
        self.htm = lsst.sphgeom.HtmPixelization(self.depth)
        converter = ConvertReferenceCatalogTask(output_dir=tempPath, config=self.config)
        dtype = [('id', '<f8'), ('ra', '<f8'), ('dec', '<f8'), ('ra_err', '<f8'), ('dec_err', '<f8'),
                 ('a', '<f8'), ('a_err', '<f8')]
        self.schema, self.key_map = converter.makeSchema(dtype)
        self.fileReader = ReadTextCatalogTask()

        self.fakeInput = self.makeSkyCatalog(outPath=None, size=5, idStart=6543)
        self.matchedPixels = np.array([1, 1, 2, 2, 3])
        self.tempDir2 = tempfile.TemporaryDirectory()
        tempPath = self.tempDir2.name
        self.filenames = {x: os.path.join(tempPath, "%d.fits" % x) for x in set(self.matchedPixels)}

        self.worker = ConvertGaiaManager(self.filenames,
                                         self.config,
                                         self.fileReader,
                                         self.indexer,
                                         self.schema,
                                         self.key_map,
                                         self.htm.universe()[0],
                                         addRefCatMetadata,
                                         self.log)

    def tearDown(self):
        self.tempDir.cleanup()
        self.tempDir2.cleanup()

    def test_positionSetting(self):
        """Test the _setProperMotion, _setParallax, and
        _setCoordinateCovariance methods.
        """
        outputCatalog = SimpleCatalog(self.worker.schema)
        outputCatalog.resize(len(self.fakeInput))

        # Set coordinate errors and covariances:
        coordErr = self.worker._getCoordErr(self.fakeInput)
        for name, array in coordErr.items():
            outputCatalog[name] = array

        for outputRow, inputRow in zip(outputCatalog, self.fakeInput):
            self.worker._setProperMotion(outputRow, inputRow)
            self.worker._setParallax(outputRow, inputRow)
            self.worker._setCoordinateCovariance(outputRow, inputRow)

        coordConvert = (self.worker.coord_err_unit).to(u.radian)
        pmConvert = (self.worker.config.pm_scale * u.milliarcsecond).to_value(u.radian)
        parallaxConvert = (self.worker.config.parallax_scale * u.milliarcsecond).to_value(u.radian)

        # Test a few combinations of coordinates, proper motion, and parallax.
        # Check that the covariance in the output catalog matches the
        # covariance calculated from the input, and also matches the covariance
        # calculated from the output catalog errors with the input correlation.
        ra_pmra_cov1 = (self.fakeInput['ra_error'] * self.fakeInput['pmra_error']
                        * self.fakeInput['ra_pmra_corr']) * coordConvert * pmConvert
        ra_pmra_cov2 = (outputCatalog['coord_raErr'] * outputCatalog['pm_raErr']
                        * self.fakeInput['ra_pmra_corr'])
        np.testing.assert_allclose(ra_pmra_cov1, outputCatalog['coord_ra_pm_ra_Cov'])
        np.testing.assert_allclose(ra_pmra_cov2, outputCatalog['coord_ra_pm_ra_Cov'])

        dec_parallax_cov1 = (self.fakeInput['dec_error'] * self.fakeInput['parallax_error']
                             * self.fakeInput['dec_parallax_corr']) * coordConvert * parallaxConvert
        dec_parallax_cov2 = (outputCatalog['coord_decErr'] * outputCatalog['parallaxErr']
                             * self.fakeInput['dec_parallax_corr'])
        np.testing.assert_allclose(dec_parallax_cov1, outputCatalog['coord_dec_parallax_Cov'])
        np.testing.assert_allclose(dec_parallax_cov2, outputCatalog['coord_dec_parallax_Cov'])

        pmdec_parallax_cov1 = (self.fakeInput['pmdec_error'] * self.fakeInput['parallax_error']
                               * self.fakeInput['parallax_pmdec_corr']) * pmConvert * parallaxConvert
        pmdec_parallax_cov2 = (outputCatalog['pm_decErr'] * outputCatalog['parallaxErr']
                               * self.fakeInput['parallax_pmdec_corr'])
        np.testing.assert_allclose(pmdec_parallax_cov1, outputCatalog['pm_dec_parallax_Cov'])
        np.testing.assert_allclose(pmdec_parallax_cov2, outputCatalog['pm_dec_parallax_Cov'])


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
