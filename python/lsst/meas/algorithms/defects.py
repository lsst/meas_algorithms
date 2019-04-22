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
"""Support for image defects"""

__all__ = ("Defects",)

import collections.abc
from deprecated.sphinx import deprecated

import lsst.geom
import lsst.pex.policy as policy
from . import Defect


@deprecated(reason="Policy defect files no longer supported (will be removed after v18)",
            category=FutureWarning)
def policyToBadRegionList(policyFile):
    """Given a Policy file describing a CCD's bad pixels, return a vector of BadRegion::Ptr"""

    badPixelsPolicy = policy.Policy.createPolicy(policyFile)
    badPixels = []

    if badPixelsPolicy.exists("Defects"):
        d = badPixelsPolicy.getArray("Defects")
        for reg in d:
            x0 = reg.get("x0")
            width = reg.get("width")
            if not width:
                x1 = reg.get("x1")
                width = x1 - x0 - 1

            y0 = reg.get("y0")
            if reg.exists("height"):
                height = reg.get("height")
            else:
                y1 = reg.get("y1")
                height = y1 - y0 - 1

            bbox = lsst.geom.BoxI(lsst.geom.PointI(x0, y0), lsst.geom.ExtentI(width, height))
            badPixels.append(Defect(bbox))

    return badPixels


class Defects(collections.abc.MutableSequence):
    """Collection of `lsst.meas.algorithms.Defect`.

    Parameters
    ----------
    defectList : iterable of `lsst.meas.algorithms.Defect`
                 or `lsst.geom.BoxI`, optional
        Collections of defects to apply to the image.
    """

    def __init__(self, defectList=None):
        self._defects = []

        if defectList is None:
            return

        # Ensure that type checking
        for d in defectList:
            self.append(d)

    def _check_value(self, value):
        """Check that the supplied value is a `~lsst.meas.algorithms.Defect`
        or can be converted to one.

        Parameters
        ----------
        value : `object`
            Value to check.

        Returns
        -------
        new : `~lsst.meas.algorithms.Defect`
            Either the supplied value or a new object derived from it.

        Raises
        ------
        ValueError
            Raised if the supplied value can not be converted to
            `~lsst.meas.algorithms.Defect`
        """
        if isinstance(value, Defect):
            pass
        elif isinstance(value, lsst.geom.BoxI):
            value = Defect(value)
        else:
            raise ValueError("Defects must be of type Defect or BoxI, not '{value}'")
        return value

    def __len__(self):
        return len(self._defects)

    def __getitem__(self, index):
        return self._defects[index]

    def __setitem__(self, index, value):
        """Can be given a `~lsst.meas.algorithms.Defect` or a `lsst.geom.BoxI`
        """
        self._defects[index] = self._check_value(value)

    def __iter__(self):
        return iter(self._defects)

    def __delitem__(self, index):
        del self._defects[index]

    def insert(self, index, value):
        self._defects.insert(index, self._check_value(value))
