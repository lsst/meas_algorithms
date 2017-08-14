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
    """!Config for InstallGaussianPsfTask"""
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


## \addtogroup LSST_task_documentation
## \{
## \page InstallGaussianPsfTask
## \ref InstallGaussianPsfTask_ "InstallGaussianPsfTask"
## \copybrief InstallGaussianPsfTask
## \}

class InstallGaussianPsfTask(pipeBase.Task):
    """!Install a Gaussian PSF model in an exposure

    @anchor InstallGaussianPsfTask_

    @section pipe_tasks_installGaussianPsf_Contents  Contents

     - @ref pipe_tasks_installGaussianPsf_Purpose
     - @ref pipe_tasks_installGaussianPsf_Initialize
     - @ref pipe_tasks_installGaussianPsf_IO
     - @ref pipe_tasks_installGaussianPsf_Config
     - @ref pipe_tasks_installGaussianPsf_Metadata
     - @ref pipe_tasks_installGaussianPsf_Debug
     - @ref pipe_tasks_installGaussianPsf_Example

    @section pipe_tasks_installGaussianPsf_Purpose  Description

    Install a Gaussian PSF model in an exposure.
    If the exposure already has a PSF model then the new model
    has the same sigma and size (width and height in pixels) of the existing model.
    If the exposure does not have a PSF model then the PSF sigma and size
    are taken from the config.

    At present the produced model is always circularly symmetric, but it is planned
    to change this to an elliptical PSF model (only for the case that the exposure
    already has a PSF model), once the necessary PSF object is available.

    A variant of this task may someday exist to estimate the PSF
    from the pixel data if no PSF model is present.

    @section pipe_tasks_installGaussianPsf_Initialize  Task initialisation

    @copydoc \_\_init\_\_

    @section pipe_tasks_installGaussianPsf_IO  Invoking the Task

    The main method is `run`.

    @section pipe_tasks_installGaussianPsf_Config  Configuration parameters

    See @ref InstallGaussianPsfConfig

    @section pipe_tasks_installGaussianPsf_Debug  Debug variables

    This task has no debug display

    @section pipe_tasks_installGaussianPsf_Example  A complete example of using InstallGaussianPsfTask

        from lsst.afw.image import ExposureF
        from lsst.meas.algorithms.installGaussianPsf import InstallGaussianPsfTask, FwhmPerSigma

        exposure = ExposureF(100, 100)
        task = InstallGaussianPsfTask()
        task.run(exposure=exposure)

        # This particular exposure had no PSF model to begin with, so the new PSF model
        # uses the config's FWHM. However, measured FWHM is based on the truncated
        # PSF image, so it does not exactly match the input
        measFwhm = exposure.getPsf().computeShape().getDeterminantRadius() * FwhmPerSigma
        assert abs(measFwhm - task.config.fwhm) < 1e-3
    """
    ConfigClass = InstallGaussianPsfConfig
    _DefaultName = "installSimplePsfModel"

    def run(self, exposure):
        """!Set exposure's PSF to a simple PSF model

        The sigma and width of the new simple PSF model matches the sigma and width of the current model,
        if any, else the config parameters are used.

        @param[in,out] exposure  exposure to which to replace or add the PSF model
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
