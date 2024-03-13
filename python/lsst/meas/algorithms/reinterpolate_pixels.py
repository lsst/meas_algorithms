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
# from lsst.meas.algorithms import Defect
from . import Defect, interpolateOverDefects


class ReinterpolatePixelsConfig(pexConfig.Config):
    """Config for ReinterpolatePixelsTask"""

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

        thresh = afwDetection.Threshold(
            mask.getPlaneBitMask(self.config.maskNameList), afwDetection.Threshold.BITMASK
        )
        fpSet = afwDetection.FootprintSet(mask, thresh)
        defectList = self._fromFootprintList(fpSet.getFootprints())
        self._interpolateDefectList(exposure, defectList)

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
        When reinterpolating following the interpolation in ip_isr, be aware that the masks are intentionally
        left grown as a side-effect of that interpolation.
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
