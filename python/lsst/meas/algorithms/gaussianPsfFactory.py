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
from __future__ import absolute_import, division

__all__ = ["GaussianPsfFactory", "SigmaPerFwhm"]

import math

from lsst.pex.config import Config, Field, ConfigurableField
from .singleGaussianPsf import SingleGaussianPsf
from .doubleGaussianPsf import DoubleGaussianPsf

SigmaPerFwhm = 1.0 / (2.0 * math.sqrt(2.0 * math.log(2.0)))


def isPositive(x):
    return x > 0


class GaussianPsfFactory(Config):
    """Factory for simple Gaussian PSF models

    Provides a high-level interface to DoubleGaussianPsf and SingleGaussianPsf
    by specifying Gaussian PSF model width in FWHM instead of sigma,
    and supporting computing kernel size as a multiple of PSF width.
    This makes it suitable for tasks where PSF width is not known in advance.
    """
    size = Field(
        doc="Kernel size (width and height) (pixels); if None then sizeFactor is used",
        dtype=int,
        optional=True,
        default=None,
        check=isPositive,
    )
    sizeFactor = Field(
        doc = "Kernel size as a factor of fwhm (dimensionless); " +
        "size = sizeFactor * fwhm; ignored if size is not None",
        dtype=float,
        optional=False,
        default=3.0,
        check=isPositive,
    )
    minSize = Field(
        doc="Minimum kernel size if using sizeFactor (pixels); ignored if size is not None",
        dtype=int,
        optional=True,
        default=5,
        check=isPositive,
    )
    maxSize = Field(
        doc="Maximum kernel size if using sizeFactor (pixels); ignored if size is not None",
        dtype=int,
        optional=True,
        default=None,
        check=isPositive,
    )
    defaultFwhm = Field(
        doc="Default FWHM of Gaussian model of core of star (pixels)",
        dtype=float,
        default=3.0,
        check=isPositive,
    )
    addWing = Field(
        doc="Add a Gaussian to represent wings?",
        dtype=bool,
        optional=False,
        default=True,
    )
    wingFwhmFactor = Field(
        doc="wing width, as a multiple of core width (dimensionless); ignored if addWing false",
        dtype=float,
        optional=False,
        default=2.5,
        check=isPositive,
    )
    wingAmplitude = Field(
        doc="wing amplitude, as a multiple of core amplitude (dimensionless); ignored if addWing false",
        dtype=float,
        optional=False,
        default=0.1,
        check=isPositive,
    )

    def computeSizeAndSigma(self, fwhm=None):
        """Compute kernel size and star width as sigma

        kernel size will be odd unless minSize or maxSize is used and that value is even.

        @param[in] fwhm: FWHM of core star (pixels); if None then defaultFwhm is used
        @return two values:
        - kernel size (width == height) in pixels
        - sigma equivalent to supplied fwhm, assuming a Gaussian (pixels)

        @warning assumes a valid config
        """
        if fwhm is None:
            fwhm = self.defaultFwhm

        if self.size is not None:
            size = self.size
        else:
            desSize = (int(self.sizeFactor * fwhm) // 2) * 2 + 1  # make result odd
            if self.minSize and self.minSize > desSize:
                size = self.minSize
            elif self.maxSize and self.maxSize < desSize:
                size = self.maxSize
            else:
                size = desSize

        return size, fwhm * SigmaPerFwhm

    def validate(self):
        Config.validate(self)
        if self.minSize and self.maxSize and self.minSize > self.maxSize:
            raise RuntimeError("minSize=%s > maxSize=%s" % (self.minSize, self.maxSize))

    def apply(self, fwhm=None):
        """Construct a GaussianPsf

        @param[in] self: an instance of ConfigClass
        @param[in] fwhm: FWHM of core of star (pixels); if None then self.defaultFwhm is used
        @return a DoubleGaussianPsf if self.addWing is True, else a SingleGaussianPsf
        """
        kernelSize, sigma = self.computeSizeAndSigma(fwhm)
        if self.addWing:
            wingsSigma = sigma * self.wingFwhmFactor
            return DoubleGaussianPsf(kernelSize, kernelSize, sigma, wingsSigma, self.wingAmplitude)
        else:
            return SingleGaussianPsf(kernelSize, kernelSize, sigma)

    @classmethod
    def makeField(cls, doc):
        """Make an lsst.pex.config.ConfigurableField
        """
        def applyWrapper(config, **kwargs):
            """Construct a Gaussian PSF

            @param[in] config: an instance of GaussianPsfFactory
            """
            return config.apply(**kwargs)
        return ConfigurableField(
            doc=doc,
            target=applyWrapper,
            ConfigClass=cls
        )
