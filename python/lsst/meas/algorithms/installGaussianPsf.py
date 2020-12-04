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
__all__ = ["InstallGaussianPsfConfig", "InstallGaussianPsfTask"]

import math

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.meas.algorithms import SingleGaussianPsf

FwhmPerSigma = 2.0*math.sqrt(2.0*math.log(2.0))


class InstallGaussianPsfConfig(pexConfig.Config):
    """Config for InstallGaussianPsfTask
    """
    fwhm = pexConfig.Field(
        dtype=float,
        default=1.5 * FwhmPerSigma,
        doc="Estimated FWHM of simple Gaussian PSF model, in pixels. "
        "Ignored if input exposure has a PSF model."
    )
    width = pexConfig.RangeField(
        dtype=int,
        doc="Width and height of PSF model, in pixels. Must be odd.",
        default=11,
        min=1,
    )

    def validate(self):
        if self.width % 2 == 0:
            raise RuntimeError("width=%s must be odd" % (self.width,))


class InstallGaussianPsfTask(pipeBase.Task):
    r"""Install a Gaussian PSF model in an exposure.
    If the exposure already has a PSF model then the new model
    has the same sigma and size (width and height in pixels) of the existing model.
    """
    ConfigClass = InstallGaussianPsfConfig
    _DefaultName = "installSimplePsfModel"

    def run(self, exposure):
        """Set exposure's PSF to a simple PSF model

        The sigma and width of the new simple PSF model matches the sigma and width of the current model,
        if any, else the config parameters are used.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure in which to replace or add the PSF model.
        """
        if exposure.hasPsf():
            psfModel = exposure.getPsf()
            psfSigma = psfModel.computeShape().getDeterminantRadius()
            width, height = psfModel.computeImage().getDimensions()
        else:
            psfSigma = self.config.fwhm / FwhmPerSigma
            width = height = self.config.width

        if psfSigma <= 0:
            raise RuntimeError("psfSigma = %s <= 0" % (psfSigma,))

        self.log.debug("installing a simple Gaussian PSF model with width=%s, height=%s, FWHM=%0.3f",
                       width, height, psfSigma*FwhmPerSigma)
        psfModel = SingleGaussianPsf(width, height, psfSigma)
        exposure.setPsf(psfModel)
