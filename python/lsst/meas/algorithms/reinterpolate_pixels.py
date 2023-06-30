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

import itertools

import lsst.afw.detection as afwDetection
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.meas.algorithms import Defect, interpolateOverDefects


class ReinterpolatePixelsConfig(pexConfig.Config):
    """Config for ReinterpolatePixelsTask"""

    growSaturatedFootprints = pexConfig.Field[int](
        doc="Number of pixels to grow footprints for saturated pixels. Be aware that when reinterpolating "
        "after the interpolation performed in ip_isr, where the masks are intentionally left grown as a side "
        "effect of interpolation, passing a value greater than zero here will further grow them. In the case "
        "of this additional growth, we will reset the mask to its original state after reinterpolation.",
        default=0,
        optional=True,
    )

    maskNameList = pexConfig.ListField[str](
        doc="Names of mask planes whose image-plane values should be interpolated.",
        default=["SAT"],
    )

    fallbackValue = pexConfig.Field[float](
        doc="Value of last resort for interpolation. If not provided (or set to None), it is assigned the "
        "clipped mean value of the exposure image.",
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

        if self.config.growSaturatedFootprints > 0 and "SAT" in self.config.maskNameList:
            # If we are interpolating over an area larger than the original
            # masked region, we need to expand the original mask bit to the
            # full area to explain why we interpolated there.
            self._growMasks(
                mask, radius=self.config.growSaturatedFootprints, maskNameList=["SAT"], maskValue="SAT"
            )

        thresh = afwDetection.Threshold(
            mask.getPlaneBitMask(self.config.maskNameList), afwDetection.Threshold.BITMASK
        )
        fpSet = afwDetection.FootprintSet(mask, thresh)
        defectList = self._fromFootprintList(fpSet.getFootprints())
        self._interpolateDefectList(exposure, defectList)

        if self.config.growSaturatedFootprints > 0 and "SAT" in self.config.maskNameList:
            # Reset the mask to its original state prior to expanding the SAT mask.
            self.originalFpSet.setMask(mask, "SAT")

    def _fromFootprintList(self, fpList):
        """Compute a defect list from a footprint list.

        Parameters
        ----------
        fpList : `list` of `lsst.afw.detection.Footprint`
            Footprint list to process.

        Returns
        -------
        defects : `list` of `lsst.meas.algorithms.Defect`
            List of defects.

        Notes
        -----
        By using `itertools.chain.from_iterable()`, the nested iterables are merged into a flat iterable,
        allowing convenient iteration over all `lsst.meas.algorithms.Defect` objects.
        """
        return list(
            itertools.chain.from_iterable(map(Defect, afwDetection.footprintToBBoxList(fp)) for fp in fpList)
        )

    def _growMasks(self, mask, radius=0, maskNameList=["BAD"], maskValue="BAD"):
        """Grow a mask by an amount and add to the requested plane.

        Parameters
        ----------
        mask : `lsst.afw.image.Mask`
            Mask image to process.
        radius : scalar
            Amount to grow the mask.
        maskNameList : `str` or `list` [`str`]
            Mask names that should be grown.
        maskValue : `str`
            Mask plane to assign the newly masked pixels to.
        """
        if radius > 0:
            thresh = afwDetection.Threshold(
                mask.getPlaneBitMask(maskNameList), afwDetection.Threshold.BITMASK
            )
            fpSet = afwDetection.FootprintSet(mask, thresh)
            self.originalFpSet = fpSet
            fpSet = afwDetection.FootprintSet(fpSet, rGrow=radius, isotropic=False)
            fpSet.setMask(mask, maskValue)

    def _interpolateDefectList(self, exposure, defectList):
        """Interpolate over defects specified in a defect list.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to process.
        defectList : `list` of `lsst.meas.algorithms.Defect`
            List of defects to interpolate over.

        Notes
        -----
        ``exposure`` is modified in place and will become the reinterpolated exposure.
        """

        # Although `interpolateOverDefects` requires the
        # `lsst.afw.detection.Psf` argument to be passed, it does not actually
        # utilize it. A dummy `lsst.afw.detection.Psf` object must still be
        # provided to prevent a TypeError.
        dummyPsf = afwDetection.Psf()
        if self.config.fallbackValue is None:
            fallbackValue = afwMath.makeStatistics(exposure.image, afwMath.MEANCLIP).getValue()
        else:
            fallbackValue = self.config.fallbackValue
        if "INTRP" not in exposure.mask.getMaskPlaneDict():
            exposure.mask.addMaskPlane("INTRP")
        interpolateOverDefects(exposure.maskedImage, dummyPsf, defectList, fallbackValue, True)
