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

__all__ = ["ComputePsfChiSqTask", "ComputePsfChiSqConfig"]

import numpy as np

from lsst.afw.detection import Psf
from lsst.afw.image import ImageF, MaskedImageF
from lsst.geom import Point2D
from lsst.pex.config import Config, Field, ListField
from lsst.pipe.base import Struct, Task


class ComputePsfChiSqConfig(Config):
    padding = Field(
        "Number of pixels to grow the PSF model bbox by on all sides.",
        dtype=int,
        default=10,
    )
    bad_mask_planes = ListField(
        "Mask planes to identify pixels to drop from the calculation.",
        dtype=str,
        default=["SAT", "SUSPECT"],
    )


class ComputePsfChiSqTask(Task):
    ConfigClass = ComputePsfChiSqConfig
    _DefaultName = "computePsfChiSq"
    config: ComputePsfChiSqConfig

    def run(
        self, *, masked_image: MaskedImageF, psf: Psf, x: np.ndarray, y: np.ndarray
    ) -> Struct:
        bitmask = masked_image.mask.getPlaneBitMask(self.config.bad_mask_planes)
        result = Struct(
            chisq=np.zeros(x.shape, dtype=np.float32),
            snr=np.zeros(x.shape, dtype=np.float32),
            n_missing_pixels=np.zeros(x.shape, dtype=np.uint16),
        )
        for i, (sx, sy) in enumerate(zip(x, y)):
            inner_model = psf.computeImage(Point2D(sx, sy))
            inner_bbox = inner_model.getBBox()
            outer_bbox = inner_bbox.dilatedBy(self.config.padding)
            n_clipped_pixels = 0
            if not masked_image.getBBox().contains(outer_bbox):
                original_bbox_area = outer_bbox.getArea()
                outer_bbox.clip(masked_image.getBBox())
                n_clipped_pixels = original_bbox_area - outer_bbox.getArea()
            outer_model = ImageF(outer_bbox)
            outer_model.array[:, :] = 0
            outer_model[inner_bbox] = inner_model.convertF()
            outer_data = masked_image.image[outer_bbox]
            outer_variance = masked_image.variance[outer_bbox]
            mask = np.logical_not(masked_image.mask[outer_bbox].array & bitmask)
            flat_model = outer_model.array[mask]
            flat_data = outer_data.array[mask]
            flat_variance = outer_variance.array[mask]
            sigma = flat_variance**0.5
            weighted_model = flat_model / sigma
            weighted_data = flat_data / sigma
            inv_var_alpha = np.dot(weighted_model, weighted_model)
            alpha = np.dot(weighted_model, weighted_data) / inv_var_alpha
            weighted_residuals = weighted_model * alpha
            weighted_residuals -= weighted_data
            result.chisq[i] = sum(weighted_residuals**2) / np.sum(mask)
            result.snr[i] = alpha * inv_var_alpha**0.5
            result.n_missing_pixels[i] = (
                outer_bbox.getArea() - np.sum(mask) + n_clipped_pixels
            )
        return result
