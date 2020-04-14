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

__all__ = ["IngestIndexCatalogTestBase", "make_coord", "makeIngestIndexConfig"]

import math
import os.path
import shutil
import string
import tempfile

import numpy as np
import astropy
import astropy.units as u

import lsst.daf.persistence as dafPersist
from lsst.meas.algorithms import IndexerRegistry
from lsst.meas.algorithms import IngestIndexedReferenceTask
import lsst.utils


def make_coord(ra, dec):
    """Make an ICRS coord given its RA, Dec in degrees."""
    return lsst.geom.SpherePoint(ra, dec, lsst.geom.degrees)


def makeIngestIndexConfig(withMagErr=False, withRaDecErr=False, withPm=False, withPmErr=False,
                          withParallax=False):
    """Make a config for IngestIndexedReferenceTask

    This is primarily intended to simplify tests of config validation,
    so fields that are not validated are not set.
    However, it can calso be used to reduce boilerplate in other tests.
    """
    config = IngestIndexedReferenceTask.ConfigClass()
    config.pm_scale = 1000.0
    config.parallax_scale = 1e3
    config.ra_name = 'ra_icrs'
    config.dec_name = 'dec_icrs'
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

    if withPmErr:
        config.pm_ra_err_name = "pm_ra_err"
        config.pm_dec_err_name = "pm_dec_err"

    if withParallax:
        config.parallax_name = "parallax"
        config.parallax_err_name = "parallax_err"

    if withPm or withParallax:
        config.epoch_name = "unixtime"
        config.epoch_format = "unix"
        config.epoch_scale = "utc"

    return config


class IngestIndexCatalogTestBase:
    """Base class for tests involving IngestIndexedReferenceTask
    """
    @classmethod
    def makeSkyCatalog(cls, outPath, size=1000, idStart=1, seed=123):
        """Make an on-sky catalog, and save it to a text file.

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
        np.random.seed(seed)
        ident = np.arange(idStart, size + idStart, dtype=int)
        ra = np.random.random(size)*360.
        dec = np.degrees(np.arccos(2.*np.random.random(size) - 1.))
        dec -= 90.
        ra_err = np.ones(size)*0.1  # arcsec
        dec_err = np.ones(size)*0.1  # arcsec
        a_mag = 16. + np.random.random(size)*4.
        a_mag_err = 0.01 + np.random.random(size)*0.2
        b_mag = 17. + np.random.random(size)*5.
        b_mag_err = 0.02 + np.random.random(size)*0.3
        is_photometric = np.random.randint(2, size=size)
        is_resolved = np.random.randint(2, size=size)
        is_variable = np.random.randint(2, size=size)
        extra_col1 = np.random.normal(size=size)
        extra_col2 = np.random.normal(1000., 100., size=size)
        # compute proper motion and PM error in arcseconds/year
        # and let the ingest task scale them to radians
        pm_amt_arcsec = cls.properMotionAmt.asArcseconds()
        pm_dir_rad = cls.properMotionDir.asRadians()
        pm_ra = np.ones(size)*pm_amt_arcsec*math.cos(pm_dir_rad)
        pm_dec = np.ones(size)*pm_amt_arcsec*math.sin(pm_dir_rad)
        pm_ra_err = np.ones(size)*cls.properMotionErr.asArcseconds()*abs(math.cos(pm_dir_rad))
        pm_dec_err = np.ones(size)*cls.properMotionErr.asArcseconds()*abs(math.sin(pm_dir_rad))
        parallax = np.ones(size)*0.1  # arcseconds
        parallax_error = np.ones(size)*0.003  # arcseconds
        unixtime = np.ones(size)*cls.epoch.unix

        def get_word(word_len):
            return "".join(np.random.choice([s for s in string.ascii_letters], word_len))
        extra_col3 = np.array([get_word(num) for num in np.random.randint(11, size=size)])

        dtype = np.dtype([('id', float), ('ra_icrs', float), ('dec_icrs', float),
                          ('ra_err', float), ('dec_err', float), ('a', float),
                          ('a_err', float), ('b', float), ('b_err', float), ('is_phot', int),
                          ('is_res', int), ('is_var', int), ('val1', float), ('val2', float),
                          ('val3', '|S11'), ('pm_ra', float), ('pm_dec', float), ('pm_ra_err', float),
                          ('pm_dec_err', float), ('parallax', float), ('parallax_error', float),
                          ('unixtime', float)])

        arr = np.array(list(zip(ident, ra, dec, ra_err, dec_err, a_mag, a_mag_err, b_mag, b_mag_err,
                                is_photometric, is_resolved, is_variable, extra_col1, extra_col2, extra_col3,
                                pm_ra, pm_dec, pm_ra_err, pm_dec_err, parallax, parallax_error, unixtime)),
                       dtype=dtype)
        if outPath is not None:
            # write the data with full precision; this is not realistic for
            # real catalogs, but simplifies tests based on round tripped data
            saveKwargs = dict(
                header="id,ra_icrs,dec_icrs,ra_err,dec_err,"
                       "a,a_err,b,b_err,is_phot,is_res,is_var,val1,val2,val3,"
                       "pm_ra,pm_dec,pm_ra_err,pm_dec_err,parallax,parallax_err,unixtime",
                fmt=["%i", "%.15g", "%.15g", "%.15g", "%.15g",
                     "%.15g", "%.15g", "%.15g", "%.15g", "%i", "%i", "%i", "%.15g", "%.15g", "%s",
                     "%.15g", "%.15g", "%.15g", "%.15g", "%.15g", "%.15g", "%.15g"]
            )

            np.savetxt(outPath+"/ref.txt", arr, delimiter=",", **saveKwargs)
            np.savetxt(outPath+"/ref_test_delim.txt", arr, delimiter="|", **saveKwargs)
            return outPath+"/ref.txt", outPath+"/ref_test_delim.txt", arr
        else:
            return arr

    @classmethod
    def tearDownClass(cls):
        try:
            shutil.rmtree(cls.outPath)
        except Exception:
            print("WARNING: failed to remove temporary dir %r" % (cls.outPath,))
        del cls.outPath
        del cls.skyCatalogFile
        del cls.skyCatalogFileDelim
        del cls.skyCatalog
        del cls.testRas
        del cls.testDecs
        del cls.searchRadius
        del cls.compCats
        del cls.testButler

    @classmethod
    def setUpClass(cls):
        cls.obs_test_dir = lsst.utils.getPackageDir('obs_test')
        cls.input_dir = os.path.join(cls.obs_test_dir, "data", "input")

        cls.outPath = tempfile.mkdtemp()
        cls.testCatPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data",
                                       "testHtmIndex.fits")
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
                    if make_coord(rec['ra_icrs'], rec['dec_icrs']).separation(cent) < cls.searchRadius:
                        cls.compCats[tupl].append(rec['id'])

        cls.testRepoPath = cls.outPath+"/test_repo"
        config = makeIngestIndexConfig(withMagErr=True, withRaDecErr=True, withPm=True, withPmErr=True,
                                       withParallax=True)
        # To match on disk test data
        config.dataset_config.indexer.active.depth = cls.depth
        config.id_name = 'id'
        config.pm_scale = 1000.0  # arcsec/yr --> mas/yr
        config.parallax_scale = 1e3  # arcsec -> milliarcsec
        # np.savetxt prepends '# ' to the header lines, so use a reader that understands that
        config.file_reader.format = 'ascii.commented_header'
        # run the intest once to create a butler repo we can compare to
        IngestIndexedReferenceTask.parseAndRun(args=[cls.input_dir, "--output", cls.testRepoPath,
                                                     cls.skyCatalogFile], config=config)
        cls.defaultDatasetName = config.dataset_config.ref_dataset_name
        cls.testDatasetName = 'diff_ref_name'
        cls.testButler = dafPersist.Butler(cls.testRepoPath)
        os.symlink(os.path.join(cls.testRepoPath, 'ref_cats', cls.defaultDatasetName),
                   os.path.join(cls.testRepoPath, 'ref_cats', cls.testDatasetName))

    def checkAllRowsInRefcat(self, refObjLoader, skyCatalog, config):
        """Check that every item in ``skyCatalog`` is in the ingested catalog,
        and check that fields are correct in it.

        Parameters
        ----------
        refObjLoader : `lsst.meas.algorithms.LoadIndexedReferenceObjectsTask`
            A reference object loader to use to search for rows from
            ``skyCatalog``.
        skyCatalog : `np.ndarray`
            The original data to compare with.
        config : `lsst.meas.algorithms.LoadIndexedReferenceObjectsConfig`
            The Config that was used to generate the refcat.
        """
        for row in skyCatalog:
            center = lsst.geom.SpherePoint(row['ra_icrs'], row['dec_icrs'], lsst.geom.degrees)
            cat = refObjLoader.loadSkyCircle(center, 2*lsst.geom.arcseconds, filterName='a').refCat
            self.assertGreater(len(cat), 0, "No objects found in loaded catalog.")
            msg = f"input row not found in loaded catalog:\nrow:\n{row}\n{row.dtype}\n\ncatalog:\n{cat[0]}"
            self.assertEqual(row['id'], cat[0]['id'], msg)
            # coordinates won't match perfectly due to rounding in radian/degree conversions
            self.assertFloatsAlmostEqual(row['ra_icrs'], cat[0]['coord_ra'].asDegrees(),
                                         rtol=1e-14, msg=msg)
            self.assertFloatsAlmostEqual(row['dec_icrs'], cat[0]['coord_dec'].asDegrees(),
                                         rtol=1e-14, msg=msg)
            # coordinate errors are not lsst.geom.Angle, so we have to use the
            # `units` field to convert them, and they are float32, so the tolerance is wider.
            raErr = cat[0]['coord_raErr']*u.Unit(cat.schema['coord_raErr'].asField().getUnits())
            decErr = cat[0]['coord_decErr']*u.Unit(cat.schema['coord_decErr'].asField().getUnits())
            self.assertFloatsAlmostEqual(row['ra_err'], raErr.to_value(config.coord_err_unit),
                                         rtol=1e-7, msg=msg)
            self.assertFloatsAlmostEqual(row['dec_err'], decErr.to_value(config.coord_err_unit),
                                         rtol=1e-7, msg=msg)

            if config.parallax_name is not None:
                self.assertFloatsAlmostEqual(row['parallax'], cat[0]['parallax'].asArcseconds())
                self.assertFloatsAlmostEqual(row['parallax_error'], cat[0]['parallaxErr'].asArcseconds())
