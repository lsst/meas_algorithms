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

__all__ = ["ConvertReferenceCatalogTestBase", "make_coord", "makeConvertConfig"]

import logging
import math
import string
import tempfile

import numpy as np
import astropy
import astropy.units as u

import lsst.daf.butler
from lsst.meas.algorithms import IndexerRegistry, ConvertRefcatManager
from lsst.meas.algorithms import ConvertReferenceCatalogConfig
import lsst.utils


def make_coord(ra, dec):
    """Make an ICRS coord given its RA, Dec in degrees."""
    return lsst.geom.SpherePoint(ra, dec, lsst.geom.degrees)


class ConvertReferenceCatalogCustomClass(ConvertRefcatManager):
    """Custom class to overload `ConvertRefcatManager._setCoordinateCovariance`
    """
    def _setCoordinateCovariance(self, record, row):
        """Coordinate covariance will not be used, so set to zero.
        """
        outputParams = ['coord_ra', 'coord_dec', 'pm_ra', 'pm_dec', 'parallax']
        for i in range(5):
            for j in range(i):
                record.set(self.key_map[f'{outputParams[j]}_{outputParams[i]}_Cov'], 0)


def makeConvertConfig(withMagErr=False, withRaDecErr=False, withPm=False,
                      withParallax=False, withFullPositionInformation=False):
    """Make a config for ConvertReferenceCatalogTask

    This is primarily intended to simplify tests of config validation,
    so fields that are not validated are not set.
    However, it can also be used to reduce boilerplate in other tests.
    """
    config = ConvertReferenceCatalogConfig()
    config.dataset_config.ref_dataset_name = "testRefCat"
    config.pm_scale = 1000.0
    config.parallax_scale = 1e3
    config.ra_name = 'ra'
    config.dec_name = 'dec'
    config.mag_column_list = ['a', 'b']

    if withMagErr:
        config.mag_err_column_map = {'a': 'a_err', 'b': 'b_err'}

    if withRaDecErr:
        config.ra_err_name = "ra_err"
        config.dec_err_name = "dec_err"
        config.coord_err_unit = "arcsecond"

    if withPm:
        config.pm_ra_name = "pm_ra"
        config.pm_dec_name = "pm_dec"
        config.pm_ra_err_name = "pm_ra_err"
        config.pm_dec_err_name = "pm_dec_err"

    if withParallax:
        config.parallax_name = "parallax"
        config.parallax_err_name = "parallax_err"

    if withPm or withParallax:
        config.epoch_name = "unixtime"
        config.epoch_format = "unix"
        config.epoch_scale = "utc"

    if withFullPositionInformation:
        config.full_position_information = True
        config.manager.retarget(ConvertReferenceCatalogCustomClass)

    return config


class ConvertReferenceCatalogTestBase:
    """Base class for tests involving ConvertReferenceCatalogTask
    """
    @classmethod
    def makeSkyCatalog(cls, outPath, size=1000, idStart=1, seed=123):
        """Make an on-sky catalog, and save it to a text file.

        The catalog columns mimic the columns from the native Gaia catalog.

        Parameters
        ----------
        outPath : `str` or None
            The directory to write the catalog to.
            Specify None to not write any output.
        size : `int`, (optional)
            Number of items to add to the catalog.
        idStart : `int`, (optional)
            First id number to put in the catalog.
        seed : `float`, (optional)
            Random seed for ``np.random``.

        Returns
        -------
        refCatPath : `str`
            Path to the created on-sky catalog.
        refCatOtherDelimiterPath : `str`
            Path to the created on-sky catalog with a different delimiter.
        refCatData : `np.ndarray`
            The data contained in the on-sky catalog files.
        """
        rng = np.random.Generator(np.random.MT19937(seed))
        ident = np.arange(idStart, size + idStart, dtype=int)
        ra = rng.random(size)*360.
        dec = np.degrees(np.arccos(2.*rng.random(size) - 1.))
        dec -= 90.
        ra_err = np.ones(size)*0.1  # arcsec
        dec_err = np.ones(size)*0.1  # arcsec
        a_mag = 16. + rng.random(size)*4.
        a_mag_err = 0.01 + rng.random(size)*0.2
        b_mag = 17. + rng.random(size)*5.
        b_mag_err = 0.02 + rng.random(size)*0.3
        is_photometric = rng.integers(2, size=size)
        is_resolved = rng.integers(2, size=size)
        is_variable = rng.integers(2, size=size)
        extra_col1 = rng.normal(size=size)
        extra_col2 = rng.normal(1000., 100., size=size)
        # compute proper motion and PM error in arcseconds/year
        # and let the convert task scale them to radians
        pm_amt_arcsec = cls.properMotionAmt.asArcseconds()
        pm_dir_rad = cls.properMotionDir.asRadians()
        pm_ra = np.ones(size)*pm_amt_arcsec*math.cos(pm_dir_rad)
        pm_dec = np.ones(size)*pm_amt_arcsec*math.sin(pm_dir_rad)
        pm_ra_err = np.ones(size)*cls.properMotionErr.asArcseconds()*abs(math.cos(pm_dir_rad))
        pm_dec_err = np.ones(size)*cls.properMotionErr.asArcseconds()*abs(math.sin(pm_dir_rad))
        parallax = np.ones(size)*0.1  # arcseconds
        parallax_error = np.ones(size)*0.003  # arcseconds
        ra_dec_corr = 2 * rng.random(size) - 1
        ra_parallax_corr = 2 * rng.random(size) - 1
        ra_pmra_corr = 2 * rng.random(size) - 1
        ra_pmdec_corr = 2 * rng.random(size) - 1
        dec_parallax_corr = 2 * rng.random(size) - 1
        dec_pmra_corr = 2 * rng.random(size) - 1
        dec_pmdec_corr = 2 * rng.random(size) - 1
        parallax_pmra_corr = 2 * rng.random(size) - 1
        parallax_pmdec_corr = 2 * rng.random(size) - 1
        pmra_pmdec_corr = 2 * rng.random(size) - 1
        unixtime = np.ones(size)*cls.epoch.unix

        def get_word(word_len):
            return "".join(rng.choice([s for s in string.ascii_letters], word_len))
        extra_col3 = np.array([get_word(num) for num in rng.integers(11, size=size)])

        dtype = np.dtype([('id', float), ('ra', float), ('dec', float),
                          ('ra_error', float), ('dec_error', float), ('a', float),
                          ('a_err', float), ('b', float), ('b_err', float), ('is_phot', int),
                          ('is_res', int), ('is_var', int), ('val1', float), ('val2', float),
                          ('val3', '|S11'), ('pmra', float), ('pmdec', float), ('pmra_error', float),
                          ('pmdec_error', float), ('parallax', float), ('parallax_error', float),
                          ('ra_dec_corr', float), ('ra_parallax_corr', float), ('ra_pmra_corr', float),
                          ('ra_pmdec_corr', float), ('dec_parallax_corr', float), ('dec_pmra_corr', float),
                          ('dec_pmdec_corr', float), ('parallax_pmra_corr', float),
                          ('parallax_pmdec_corr', float), ('pmra_pmdec_corr', float), ('unixtime', float)])

        arr = np.array(list(zip(ident, ra, dec, ra_err, dec_err, a_mag, a_mag_err, b_mag, b_mag_err,
                                is_photometric, is_resolved, is_variable, extra_col1, extra_col2, extra_col3,
                                pm_ra, pm_dec, pm_ra_err, pm_dec_err, parallax, parallax_error, ra_dec_corr,
                                ra_parallax_corr, ra_pmra_corr, ra_pmdec_corr, dec_parallax_corr,
                                dec_pmra_corr, dec_pmdec_corr, parallax_pmra_corr, parallax_pmdec_corr,
                                pmra_pmdec_corr, unixtime)),
                       dtype=dtype)
        if outPath is not None:
            # write the data with full precision; this is not realistic for
            # real catalogs, but simplifies tests based on round tripped data
            saveKwargs = dict(
                header="id,ra,dec,ra_err,dec_err,"
                       "a,a_err,b,b_err,is_phot,is_res,is_var,val1,val2,val3,"
                       "pm_ra,pm_dec,pm_ra_err,pm_dec_err,parallax,parallax_err,ra_dec_corr,"
                       "ra_parallax_corr,ra_pmra_corr,ra_pmdec_corr,dec_parallax_corr,"
                       "dec_pmra_corr,dec_pmdec_corr,parallax_pmra_corr,parallax_pmdec_corr,"
                       "pmra_pmdec_corr,unixtime",
                fmt=["%i", "%.15g", "%.15g", "%.15g", "%.15g",
                     "%.15g", "%.15g", "%.15g", "%.15g", "%i", "%i", "%i", "%.15g", "%.15g", "%s",
                     "%.15g", "%.15g", "%.15g", "%.15g", "%.15g", "%.15g", "%.15g", "%.15g", "%.15g", "%.15g",
                     "%.15g", "%.15g", "%.15g", "%.15g", "%.15g", "%.15g", "%.15g"]
            )

            np.savetxt(outPath+"/ref.txt", arr, delimiter=",", **saveKwargs)
            np.savetxt(outPath+"/ref_test_delim.txt", arr, delimiter="|", **saveKwargs)
            return outPath+"/ref.txt", outPath+"/ref_test_delim.txt", arr
        else:
            return arr

    @classmethod
    def tearDownClass(cls):
        cls.outDir.cleanup()
        del cls.outPath
        del cls.skyCatalogFile
        del cls.skyCatalogFileDelim
        del cls.skyCatalog
        del cls.testRas
        del cls.testDecs
        del cls.searchRadius
        del cls.compCats

    @classmethod
    def setUpClass(cls):
        cls.outDir = tempfile.TemporaryDirectory()
        cls.outPath = cls.outDir.name
        # arbitrary, but reasonable, amount of proper motion (angle/year)
        # and direction of proper motion
        cls.properMotionAmt = 3.0*lsst.geom.arcseconds
        cls.properMotionDir = 45*lsst.geom.degrees
        cls.properMotionErr = 1e-3*lsst.geom.arcseconds
        cls.epoch = astropy.time.Time(58206.861330339219, scale="tai", format="mjd")
        cls.skyCatalogFile, cls.skyCatalogFileDelim, cls.skyCatalog = cls.makeSkyCatalog(cls.outPath)
        cls.testRas = [210., 14.5, 93., 180., 286., 0.]
        cls.testDecs = [-90., -51., -30.1, 0., 27.3, 62., 90.]
        cls.searchRadius = 3. * lsst.geom.degrees
        cls.compCats = {}  # dict of center coord: list of IDs of stars within cls.searchRadius of center
        cls.depth = 4  # gives a mean area of 20 deg^2 per pixel, roughly matching a 3 deg search radius

        config = IndexerRegistry['HTM'].ConfigClass()
        # Match on disk comparison file
        config.depth = cls.depth
        cls.indexer = IndexerRegistry['HTM'](config)
        for ra in cls.testRas:
            for dec in cls.testDecs:
                tupl = (ra, dec)
                cent = make_coord(*tupl)
                cls.compCats[tupl] = []
                for rec in cls.skyCatalog:
                    if make_coord(rec['ra'], rec['dec']).separation(cent) < cls.searchRadius:
                        cls.compCats[tupl].append(rec['id'])

        cls.testRepoPath = cls.outPath+"/test_repo"

    def setUp(self):
        self.repoPath = tempfile.TemporaryDirectory()  # cleaned up automatically when test ends
        self.butler = self.makeTemporaryRepo(self.repoPath.name, self.depth)
        self.logger = logging.getLogger('lsst.ReferenceObjectLoader')

    def tearDown(self):
        self.repoPath.cleanup()

    @staticmethod
    def makeTemporaryRepo(rootPath, depth):
        """Create a temporary butler repository, configured to support a given
        htm pixel depth, to use for a single test.

        Parameters
        ----------
        rootPath : `str`
            Root path for butler.
        depth : `int`
            HTM pixel depth to be used in this test.

        Returns
        -------
        butler : `lsst.daf.butler.Butler`
            The newly created and instantiated butler.
        """
        dimensionConfig = lsst.daf.butler.DimensionConfig()
        dimensionConfig['skypix']['common'] = f'htm{depth}'
        lsst.daf.butler.Butler.makeRepo(rootPath, dimensionConfig=dimensionConfig)
        return lsst.daf.butler.Butler(rootPath, writeable=True)

    def checkAllRowsInRefcat(self, refObjLoader, skyCatalog, config):
        """Check that every item in ``skyCatalog`` is in the converted catalog,
        and check that fields are correct in it.

        Parameters
        ----------
        refObjLoader : `lsst.meas.algorithms.ReferenceObjectLoader`
            A reference object loader to use to search for rows from
            ``skyCatalog``.
        skyCatalog : `np.ndarray`
            The original data to compare with.
        config : `lsst.meas.algorithms.LoadReferenceObjectsConfig`
            The Config that was used to generate the refcat.
        """
        for row in skyCatalog:
            center = lsst.geom.SpherePoint(row['ra'], row['dec'], lsst.geom.degrees)
            with self.assertLogs(self.logger.name, level="INFO") as cm:
                cat = refObjLoader.loadSkyCircle(center, 2*lsst.geom.arcseconds, filterName='a').refCat
            self.assertIn("Loading reference objects from testRefCat in region", cm.output[0])
            self.assertGreater(len(cat), 0, "No objects found in loaded catalog.")
            msg = f"input row not found in loaded catalog:\nrow:\n{row}\n{row.dtype}\n\ncatalog:\n{cat[0]}"
            self.assertEqual(row['id'], cat[0]['id'], msg)
            # coordinates won't match perfectly due to rounding in radian/degree conversions
            self.assertFloatsAlmostEqual(row['ra'], cat[0]['coord_ra'].asDegrees(),
                                         rtol=1e-14, msg=msg)
            self.assertFloatsAlmostEqual(row['dec'], cat[0]['coord_dec'].asDegrees(),
                                         rtol=1e-14, msg=msg)
            if config.coord_err_unit is not None:
                # coordinate errors are not lsst.geom.Angle, so we have to use the
                # `units` field to convert them, and they are float32, so the tolerance is wider.
                raErr = cat[0]['coord_raErr']*u.Unit(cat.schema['coord_raErr'].asField().getUnits())
                decErr = cat[0]['coord_decErr']*u.Unit(cat.schema['coord_decErr'].asField().getUnits())
                self.assertFloatsAlmostEqual(row['ra_error'], raErr.to_value(config.coord_err_unit),
                                             rtol=1e-7, msg=msg)
                self.assertFloatsAlmostEqual(row['dec_error'], decErr.to_value(config.coord_err_unit),
                                             rtol=1e-7, msg=msg)

            if config.parallax_name is not None:
                self.assertFloatsAlmostEqual(row['parallax'], cat[0]['parallax'].asArcseconds())
                parallaxErr = cat[0]['parallaxErr'].asArcseconds()
                # larger tolerance: input data is float32
                self.assertFloatsAlmostEqual(row['parallax_error'], parallaxErr, rtol=3e-8)
