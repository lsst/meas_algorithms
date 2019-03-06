#!/usr/bin/env python
"""Write a small trivial reference catalog, to test reading catalog versions.

Takes no arguments, writes a `versionX/` directory containing the new test
repo to the current directory. Requires that obs_test be setup.

This code was used to write the version0/ and version1/ test repos in
`meas_algorithms/tests/data/` (after suitably configuring the code to produce
a "version 0" catalog that looks like most of our existing catalogs (e.g.
fluxSigma instead of fluxErr, no units for fluxes).
"""
import os.path

import numpy as np
import astropy

from lsst.meas.algorithms import IngestIndexedReferenceTask
from lsst.meas.algorithms.ingestIndexReferenceTask import LATEST_FORMAT_VERSION
import lsst.utils


def make_skyCatalog(outPath, size=100):
    """Write a text file containing an on-sky catalog to be ingested."""
    np.random.seed(123)
    epoch = astropy.time.Time(58206.861330339219, scale="tai", format="mjd")
    ident = np.arange(1, size+1, dtype=int)
    # ra/dec in degrees, centered on (10,20)
    ra = 10 + np.random.random(size)
    dec = 20 + np.random.random(size)
    ra_err = np.ones(size)*0.1  # arcsec
    dec_err = np.ones(size)*0.1  # arcsec
    a_mag = 16. + np.random.random(size)*4.
    a_mag_err = 0.01 + np.random.random(size)*0.2
    b_mag = 17. + np.random.random(size)*5.
    b_mag_err = 0.02 + np.random.random(size)*0.3
    is_photometric = np.random.randint(2, size=size)
    is_resolved = np.random.randint(2, size=size)
    is_variable = np.random.randint(2, size=size)
    # compute proper motion and PM error in arcseconds/year
    # and let the ingest task scale them to radians
    pm_amt_arcsec = (3.0*lsst.geom.arcseconds).asArcseconds()
    pm_dir_rad = (45*lsst.geom.degrees).asRadians()
    pm_ra = np.ones(size)*pm_amt_arcsec*np.cos(pm_dir_rad)
    pm_dec = np.ones(size)*pm_amt_arcsec*np.sin(pm_dir_rad)
    properMotionErr = 1e-3*lsst.geom.arcseconds
    pm_ra_err = np.ones(size)*properMotionErr.asArcseconds()*abs(np.cos(pm_dir_rad))
    pm_dec_err = np.ones(size)*properMotionErr.asArcseconds()*abs(np.sin(pm_dir_rad))
    unixtime = np.ones(size)*epoch.unix

    dtype = np.dtype([('id', float), ('ra_icrs', float), ('dec_icrs', float),
                     ('ra_err', float), ('dec_err', float), ('a', float),
                      ('a_err', float), ('b', float), ('b_err', float), ('is_phot', int),
                      ('is_res', int), ('is_var', int), ('pm_ra', float),
                      ('pm_dec', float), ('pm_ra_err', float),
                      ('pm_dec_err', float), ('unixtime', float)])

    arr = np.array(list(zip(ident, ra, dec, ra_err, dec_err, a_mag, a_mag_err, b_mag, b_mag_err,
                            is_photometric, is_resolved, is_variable,
                            pm_ra, pm_dec, pm_ra_err, pm_dec_err, unixtime)), dtype=dtype)
    # write the data with full precision; this is not realistic for
    # real catalogs, but simplifies tests based on round tripped data
    saveKwargs = dict(
        header="id,ra_icrs,dec_icrs,ra_err,dec_err,"
               "a,a_err,b,b_err,is_phot,is_res,is_var,"
               "pm_ra,pm_dec,pm_ra_err,pm_dec_err,unixtime",
        fmt=["%i", "%.15g", "%.15g", "%.15g", "%.15g",
             "%.15g", "%.15g", "%.15g", "%.15g", "%i", "%i", "%i",
             "%.15g", "%.15g", "%.15g", "%.15g", "%.15g"]
    )

    np.savetxt(outPath+"/ref.txt", arr, delimiter=",", **saveKwargs)
    return outPath+"/ref.txt"


inPath = os.path.join(lsst.utils.getPackageDir('meas_algorithms'), 'tests/data')
outPath = os.path.join(inPath, "version%s" % LATEST_FORMAT_VERSION)

skyCatalogFile = make_skyCatalog(inPath)

config = IngestIndexedReferenceTask.ConfigClass()
config.ra_name = 'ra_icrs'
config.dec_name = 'dec_icrs'
config.mag_column_list = ['a', 'b']
config.mag_err_column_map = {'a': 'a_err', 'b': 'b_err'}
config.dataset_config.indexer.active.depth = 4

OBS_TEST_DIR = lsst.utils.getPackageDir('obs_test')
obsTestPath = os.path.join(OBS_TEST_DIR, "data", "input")
IngestIndexedReferenceTask.parseAndRun(
    args=[obsTestPath, "--output", outPath, skyCatalogFile],
    config=config)

# cleanup files and make the mapper independent of obs_test's path
os.remove(skyCatalogFile)
os.remove(os.path.join(outPath, 'repositoryCfg.yaml'))
with open(os.path.join(outPath, '_mapper'), 'w') as mapperFile:
    mapperFile.write('lsst.obs.test.TestMapper\n')

print("Wrote new catalog to:", outPath)
