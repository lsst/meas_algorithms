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

__all__ = ("ReinterpolatePixelsConfig", "ReinterpolatePixelsTask")


import math

import lsst.afw.detection as afwDetection
import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.afw.image import MaskedImageF


class ReinterpolatePixelsConfig(pexConfig.Config):
    """Config for ReinterpolatePixelsTask"""

    kernelFwhm = pexConfig.Field[float](
        doc="FWHM of double Gaussian smoothing kernel.",
        default=1.0,
    )

    growFootprints = pexConfig.Field[int](
        doc="Number of pixels to grow footprints for pixels.",
        default=1,
    )

    maskNameList = pexConfig.ListField[str](
        doc="Mask plane name.",
        default=["SAT"],
    )

    fallbackValue = pexConfig.Field[float](
        doc="Value of last resort for interpolation.",
        default=None,
        optional=True,
    )


class ReinterpolatePixelsTask(pipeBase.Task):
    """Task to reinterpolate pixels"""

    ConfigClass = ReinterpolatePixelsConfig
    _DefaultName = "reinterpolatePixels"

    def run(self, exposure):
        """Reinterpolate pixels over ``exposure``.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to interpolate over.

        Notes
        -----
        ``exposure`` is modified in place and will become the reinterpolated exposure.
        """
        mask = exposure.mask

        if self.config.growFootprints > 0 and "SAT" in self.config.maskNameList:
            # If we are interpolating over an area larger than the original
            # masked region, we need to expand the original mask bit to the
            # full area to explain why we interpolated there.
            from lsst.ip.isr import growMasks

            growMasks(mask, radius=self.config.growFootprints, maskNameList=["SAT"], maskValue="SAT")

        from lsst.ip.isr import Defects

        thresh = afwDetection.Threshold(
            mask.getPlaneBitMask(self.config.maskNameList), afwDetection.Threshold.BITMASK
        )
        fpSet = afwDetection.FootprintSet(mask, thresh)
        defectList = Defects.fromFootprintList(fpSet.getFootprints())
        self._interpolateDefectList(exposure, defectList)

    def _interpolateDefectList(self, exposure, defectList):
        """Interpolate over defects specified in a defect list.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to process.
        defectList : `lsst.meas.algorithms.Defects`
            List of defects to interpolate over.

        Notes
        -----
        ``exposure`` is modified in place and will become the reinterpolated exposure.
        """

        psf = self._createPsf(self.config.kernelFwhm)
        if self.config.fallbackValue is None:
            fallbackValue = afwMath.makeStatistics(exposure.image, afwMath.MEANCLIP).getValue()
        else:
            fallbackValue = self.config.fallbackValue
        if "INTRP" not in exposure.mask.getMaskPlaneDict():
            exposure.mask.addMaskPlane("INTRP")
        maskedImage = MaskedImageF(image=exposure.image, mask=exposure.mask, variance=exposure.variance)
        measAlg.interpolateOverDefects(maskedImage, psf, defectList, fallbackValue, True)
        exposure.image = maskedImage.image
        exposure.mask = maskedImage.mask
        exposure.variance = maskedImage.variance

    @staticmethod
    def _createPsf(fwhm):
        """Make a double Gaussian PSF.

        Parameters
        ----------
        fwhm : float
            FWHM of double Gaussian smoothing kernel.

        Returns
        -------
        psf : `lsst.meas.algorithms.DoubleGaussianPsf`
            The created smoothing kernel.
        """
        ksize = 4 * int(fwhm) + 1
        return measAlg.DoubleGaussianPsf(ksize, ksize, fwhm / (2 * math.sqrt(2 * math.log(2))))
