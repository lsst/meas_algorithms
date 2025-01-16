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

__all__ = ["MaskBrightStarsTask", "MaskBrightStarsConfig"]

import math

import numpy as np
from astropy.modeling import models, fitting

import lsst.afw.math
import lsst.pex.config as pexConfig
import lsst.pipe.base


class MaskBrightStarsConfig(pexConfig.Config):
    scale = pexConfig.Field(
        doc="Signal-to-noise to grow each mask to; zero means do not mask.",
        dtype=float,
        default=2
    )
    name = pexConfig.Field(
        doc="Name of mask plane to set for pixels affected by bright sources.",
        dtype=str,
        default="BRIGHT"
    )

import lsst.afw.display
display = lsst.afw.display.Display()
display.frame = 1


class MaskBrightStarsTask(lsst.pipe.base.Task):
    ConfigClass = MaskBrightStarsConfig
    _DefaultName = "maskBrightStars"

    def run(self, exposure, stars, flux_field):
        mask = 2**exposure.mask.addMaskPlane(self.config.name)

        footprints = []
        stdev = lsst.afw.math.makeStatistics(exposure.maskedImage, lsst.afw.math.STDEVCLIP).getValue()

        for star in stars:
            result = self._mask_star(exposure, stdev, mask, star, flux_field)
            if result is not None:
                footprints.append(result)

        return lsst.pipe.base.Struct(masked_footprints=footprints)

    def _mask_star(self, exposure, stdev, mask, star, flux_field):
        centroid = lsst.geom.Point2D(star["slot_Centroid_x"], star["slot_Centroid_y"])
        instflux = exposure.photoCalib.nanojanskyToInstFlux(star[flux_field], centroid)
        sigma = exposure.psf.computeShape(centroid).getDeterminantRadius()

        psf = exposure.psf.computeImage(centroid)
        y, x = np.mgrid[:psf.getWidth(), :psf.getHeight()]
        model_init = models.Gaussian2D(1, psf.getWidth()/2, psf.getHeight()/2, sigma, sigma) + \
            models.Gaussian2D(.5, psf.getWidth()/2, psf.getHeight()/2, 2*sigma, 2*sigma)
        fitter = fitting.SLSQPLSQFitter()
        model_fit = fitter(model_init, x, y, psf.array)
        sigma_2 = max(model_fit.x_stddev_0, model_fit.y_stddev_0, model_fit.x_stddev_1, model_fit.y_stddev_1)
        print(model_fit, sigma_2)
        # import os; print(os.getpid()); import ipdb; ipdb.set_trace();

        # TODO: is it 2pi or 4pi here? Pretty sure it's 2pi.
        amplitude = instflux / (sigma**2 * 2*math.pi)
        try:
            radius = math.sqrt(-math.log(self.config.scale / amplitude) * 2.0*2*sigma_2**2)
        except ValueError:
            self.log.debug("Skipping source at (%s) due to flux greater than configured S/N scale (%s > %s)",
                           centroid, amplitude, self.config.scale)
            return
        print(centroid, instflux, sigma, sigma_2.value, amplitude, radius)

        # ceil, so that we're always over-estimating the affected mask.
        spans = lsst.afw.geom.SpanSet.fromShape(math.ceil(radius), offset=lsst.geom.Point2I(centroid))
        footprint = lsst.afw.detection.Footprint(spans, exposure.getBBox())
        footprint.clipTo(exposure.getBBox())
        # TODO: check if there's anything left first?
        footprint.spans.setMask(exposure.mask, mask)
        # import os; print(os.getpid()); import ipdb; ipdb.set_trace();

        return footprint
