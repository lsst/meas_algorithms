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

__all__ = ["plantSources", "makeRandomTransmissionCurve", "makeDefectList",
           "MockReferenceObjectLoaderFromFiles", "MockRefcatDataId"]

import numpy as np
import esutil

import lsst.geom
import lsst.afw.image as afwImage
from lsst.pipe.base import InMemoryDatasetHandle
from lsst import sphgeom
from . import SingleGaussianPsf
from . import Defect

from . import ReferenceObjectLoader
import lsst.afw.table as afwTable


def plantSources(bbox, kwid, sky, coordList, addPoissonNoise=True):
    """Make an exposure with stars (modelled as Gaussians)

    Parameters
    ----------
    bbox : `lsst.geom.Box2I`
        Parent bbox of exposure
    kwid : `int`
        Kernal width (and height; kernal is square)
    sky : `float`
        Amount of sky background (counts)
    coordList : `list [tuple]`
        A list of [x, y, counts, sigma] where:
            * x,y are relative to exposure origin
            * counts is the integrated counts for the star
            * sigma is the Gaussian sigma in pixels
    addPoissonNoise : `bool`
        If True: add Poisson noise to the exposure
    """
    # make an image with sources
    img = afwImage.ImageD(bbox)
    meanSigma = 0.0
    for coord in coordList:
        x, y, counts, sigma = coord
        meanSigma += sigma

        # make a single gaussian psf
        psf = SingleGaussianPsf(kwid, kwid, sigma)

        # make an image of it and scale to the desired number of counts
        thisPsfImg = psf.computeImage(lsst.geom.PointD(x, y))
        thisPsfImg *= counts

        # bbox a window in our image and add the fake star image
        psfBox = thisPsfImg.getBBox()
        psfBox.clip(bbox)
        if psfBox != thisPsfImg.getBBox():
            thisPsfImg = thisPsfImg[psfBox, afwImage.PARENT]
        imgSeg = img[psfBox, afwImage.PARENT]
        imgSeg += thisPsfImg
    meanSigma /= len(coordList)

    img += sky

    # add Poisson noise
    if (addPoissonNoise):
        np.random.seed(seed=1)  # make results reproducible
        imgArr = img.getArray()
        imgArr[:] = np.random.poisson(imgArr)

    # bundle into a maskedimage and an exposure
    mask = afwImage.Mask(bbox)
    var = img.convertFloat()
    img -= sky
    mimg = afwImage.MaskedImageF(img.convertFloat(), mask, var)
    exposure = afwImage.makeExposure(mimg)

    # insert an approximate psf
    psf = SingleGaussianPsf(kwid, kwid, meanSigma)
    exposure.setPsf(psf)

    return exposure


def makeRandomTransmissionCurve(rng, minWavelength=4000.0, maxWavelength=7000.0, nWavelengths=200,
                                maxRadius=80.0, nRadii=30, perturb=0.05):
    """Create a random TransmissionCurve with nontrivial spatial and
    wavelength variation.

    Parameters
    ----------
    rng : numpy.random.RandomState
        Random number generator.
    minWavelength : float
        Average minimum wavelength for generated TransmissionCurves (will be
        randomly perturbed).
    maxWavelength : float
        Average maximum wavelength for generated TransmissionCurves (will be
        randomly perturbed).
    nWavelengths : int
        Number of samples in the wavelength dimension.
    maxRadius : float
        Average maximum radius for spatial variation (will be perturbed).
    nRadii : int
        Number of samples in the radial dimension.
    perturb: float
        Fraction by which wavelength and radius bounds should be randomly
        perturbed.
    """
    dWavelength = maxWavelength - minWavelength

    def perturbed(x, s=perturb*dWavelength):
        return x + 2.0*s*(rng.rand() - 0.5)

    wavelengths = np.linspace(perturbed(minWavelength), perturbed(maxWavelength), nWavelengths)
    radii = np.linspace(0.0, perturbed(maxRadius, perturb*maxRadius), nRadii)
    throughput = np.zeros(wavelengths.shape + radii.shape, dtype=float)
    # throughput will be a rectangle in wavelength, shifting to higher wavelengths and shrinking
    # in height with radius, going to zero at all bounds.
    peak0 = perturbed(0.9, 0.05)
    start0 = perturbed(minWavelength + 0.25*dWavelength)
    stop0 = perturbed(minWavelength + 0.75*dWavelength)
    for i, r in enumerate(radii):
        mask = np.logical_and(wavelengths >= start0 + r, wavelengths <= stop0 + r)
        throughput[mask, i] = peak0*(1.0 - r/1000.0)
    return afwImage.TransmissionCurve.makeRadial(throughput, wavelengths, radii)


def makeDefectList():
    """Create a list of defects that can be used for testing.

    Returns
    -------
    defectList = `list` [`lsst.meas.algorithms.Defect`]
        The list of defects.
    """
    defectList = [Defect(lsst.geom.Box2I(lsst.geom.Point2I(962, 0),
                                         lsst.geom.Extent2I(2, 4611))),
                  Defect(lsst.geom.Box2I(lsst.geom.Point2I(1316, 0),
                                         lsst.geom.Extent2I(2, 4611))),
                  Defect(lsst.geom.Box2I(lsst.geom.Point2I(1576, 0),
                                         lsst.geom.Extent2I(4, 4611))),
                  Defect(lsst.geom.Box2I(lsst.geom.Point2I(1626, 0),
                                         lsst.geom.Extent2I(2, 4611))),
                  Defect(lsst.geom.Box2I(lsst.geom.Point2I(1994, 252),
                                         lsst.geom.Extent2I(2, 4359))),
                  Defect(lsst.geom.Box2I(lsst.geom.Point2I(1426, 702),
                                         lsst.geom.Extent2I(2, 3909))),
                  Defect(lsst.geom.Box2I(lsst.geom.Point2I(1526, 1140),
                                         lsst.geom.Extent2I(2, 3471))),
                  Defect(lsst.geom.Box2I(lsst.geom.Point2I(856, 2300),
                                         lsst.geom.Extent2I(2, 2311))),
                  Defect(lsst.geom.Box2I(lsst.geom.Point2I(858, 2328),
                                         lsst.geom.Extent2I(2, 65))),
                  Defect(lsst.geom.Box2I(lsst.geom.Point2I(859, 2328),
                                         lsst.geom.Extent2I(1, 56))),
                  Defect(lsst.geom.Box2I(lsst.geom.Point2I(844, 2796),
                                         lsst.geom.Extent2I(4, 1814))),
                  Defect(lsst.geom.Box2I(lsst.geom.Point2I(1366, 2804),
                                         lsst.geom.Extent2I(2, 1806))),
                  Defect(lsst.geom.Box2I(lsst.geom.Point2I(1766, 3844),
                                         lsst.geom.Extent2I(2, 766))),
                  Defect(lsst.geom.Box2I(lsst.geom.Point2I(1872, 4228),
                                         lsst.geom.Extent2I(2, 382))),
                  ]

    return defectList


class MockRefcatDataId:
    """Mock reference catalog dataId.

    The reference catalog dataId is only used to retrieve a region property.

    Parameters
    ----------
    region : `lsst.sphgeom.Region`
        The region associated with this mock dataId.
    """
    def __init__(self, region):
        self._region = region

    @property
    def region(self):
        return self._region


class MockReferenceObjectLoaderFromFiles(ReferenceObjectLoader):
    """A simple mock of ReferenceObjectLoader.

    This mock ReferenceObjectLoader uses a set of files on disk to create
    mock dataIds and data reference handles that can be accessed
    without a butler. The files must be afw catalog files in the reference
    catalog format, sharded with HTM pixelization.

    Parameters
    ----------
    filenames : `list` [`str`]
        Names of files to use.
    config : `lsst.meas.astrom.LoadReferenceObjectsConfig`, optional
        Configuration object if necessary to override defaults.
    htmLevel : `int`, optional
        HTM level to use for the loader.
    """
    def __init__(self, filenames, name='cal_ref_cat', config=None, htmLevel=4):
        dataIds, refCats = self._createDataIdsAndRefcats(filenames, htmLevel, name)

        super().__init__(dataIds, refCats, name=name, config=config)

    def _createDataIdsAndRefcats(self, filenames, htmLevel, name):
        """Create mock dataIds and refcat handles.

        Parameters
        ----------
        filenames : `list` [`str`]
            Names of files to use.
        htmLevel : `int`
            HTM level to use for the loader.
        name : `str`
            Name of reference catalog (for logging).

        Returns
        -------
        dataIds : `list` [`MockRefcatDataId`]
            List of mock dataIds.
        refCats : `list` [`lsst.pipe.base.InMemoryDatasetHandle`]
            List of mock deferred dataset handles.

        Raises
        ------
        RuntimeError if any file contains sources that cover more than one HTM
            pixel at level ``htmLevel``.
        """
        pixelization = sphgeom.HtmPixelization(htmLevel)
        htm = esutil.htm.HTM(htmLevel)

        dataIds = []
        refCats = []

        for filename in filenames:
            cat = afwTable.BaseCatalog.readFits(filename)

            ids = htm.lookup_id(np.rad2deg(cat['coord_ra']), np.rad2deg(cat['coord_dec']))

            if len(np.unique(ids)) != 1:
                raise RuntimeError(f"File {filename} contains more than one pixel at level {htmLevel}")

            dataIds.append(MockRefcatDataId(pixelization.pixel(ids[0])))
            refCats.append(InMemoryDatasetHandle(cat, name=name))

        return dataIds, refCats
