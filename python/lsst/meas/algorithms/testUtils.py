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

import numpy as np

import lsst.geom
import lsst.afw.image as afwImage
from . import SingleGaussianPsf
from . import Defect


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
