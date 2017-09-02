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
import numpy

import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
from . import SingleGaussianPsf


def plantSources(bbox, kwid, sky, coordList, addPoissonNoise=True):
    """Make an exposure with stars (modelled as Gaussians)

    @param bbox: parent bbox of exposure
    @param kwid: kernel width (and height; kernel is square)
    @param sky: amount of sky background (counts)
    @param coordList: a list of [x, y, counts, sigma], where:
        * x,y are relative to exposure origin
        * counts is the integrated counts for the star
        * sigma is the Gaussian sigma in pixels
    @param addPoissonNoise: add Poisson noise to the exposure?
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
        thisPsfImg = psf.computeImage(afwGeom.PointD(int(x), int(y)))
        thisPsfImg *= counts

        # bbox a window in our image and add the fake star image
        imgSeg = img.Factory(img, thisPsfImg.getBBox())
        imgSeg += thisPsfImg
    meanSigma /= len(coordList)

    img += sky

    # add Poisson noise
    if (addPoissonNoise):
        numpy.random.seed(seed=1)  # make results reproducible
        imgArr = img.getArray()
        imgArr[:] = numpy.random.poisson(imgArr)

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
