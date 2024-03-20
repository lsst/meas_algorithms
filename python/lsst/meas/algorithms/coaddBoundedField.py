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

import numpy as np

from lsst.geom import Point2D
from lsst.utils import continueClass

from ._algorithmsLib import CoaddBoundedField

__all__ = ["CoaddBoundedField"]


@continueClass
class CoaddBoundedField:  # noqa: F811
    def evaluate(self, x, y=None):
        """Evaluate the CoaddBoundedField.

        This accepts either a Point2D or an array of x and y
        positions.

        When arrays are passed, this uses a vectorized version
        of CoaddBoundedField::evaluate(). If the coadd bounded
        field has throwOnMissing then this will return NaN
        for missing values; otherwise it will return 0.0.

        Parameters
        ----------
        x : `lsst.geom.Point2D` or `np.ndarray`
            Array of x values.
        y : `np.ndarray`, optional
            Array of y values.

        Returns
        -------
        values : `float` or `np.ndarray`
            Evaluated value or array of values.
        """
        if isinstance(x, Point2D):
            return self._evaluate(x)

        _x = np.atleast_1d(x).ravel()
        _y = np.atleast_1d(y).ravel()

        if len(_x) != len(_y):
            raise ValueError("x and y arrays must be the same length.")

        ra, dec = self.getCoaddWcs().pixelToSkyArray(_x, _y)

        sums = np.zeros(len(_x))
        wts = np.zeros_like(sums)

        for elt in self.getElements():
            ix, iy = elt.wcs.skyToPixelArray(ra, dec)
            in_box = elt.field.getBBox().contains(ix, iy)
            if elt.validPolygon:
                in_box &= elt.validPolygon.contains(ix, iy)
            sums[in_box] += elt.weight*elt.field.evaluate(ix[in_box], iy[in_box])
            wts[in_box] += elt.weight

        good = (wts > 0)
        values = np.zeros(len(_x))
        values[good] = sums[good]/wts[good]

        if self.getThrowOnMissing():
            values[~good] = np.nan

        return values
