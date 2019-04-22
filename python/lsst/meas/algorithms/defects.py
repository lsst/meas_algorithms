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

import logging
import collections.abc
from deprecated.sphinx import deprecated
import numpy as np

import lsst.geom
import lsst.pex.policy as policy
import lsst.afw.table

from . import Defect

log = logging.getLogger(__name__)


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

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False

        # Assume equal if bounding boxes are equal
        for d1, d2 in zip(self, other):
            if d1.getBBox() != d2.getBBox():
                return False

        return True

    def insert(self, index, value):
        self._defects.insert(index, self._check_value(value))

    def toTable(self):
        """Convert defect list to `~lsst.afw.table.BaseCatalog`

        Returns
        -------
        table : `lsst.afw.table.BaseCatalog`
            Defects in tabular form.
        """
        schema = lsst.afw.table.Schema()
        x = schema.addField("X", type="D", units="pix")
        y = schema.addField("Y", type="D", units="pix")
        shape = schema.addField("SHAPE", type="String", size=16)
        r = schema.addField("R", type="ArrayD", size=2, units="pix")
        rotang = schema.addField("ROTANG", type="D", units="deg")
        component = schema.addField("COMPONENT", type="I")
        table = lsst.afw.table.BaseCatalog(schema)

        for i, defect in enumerate(self._defects):
            box = defect.getBBox()
            record = table.addNew()
            record.set(x, box.getCenterX())
            record.set(y, box.getCenterY())
            record.set(shape, "BOX")
            record.set(r, np.array([box.getWidth(), box.getHeight()], dtype=np.float64))
            record.set(rotang, 0.0)
            record.set(component, i)

        return table

    def writeFits(self, *args):
        """Write defect list to FITS.

        Parameters
        ----------
        *args
            Arguments to be forwarded to
            `lsst.afw.table.BaseCatalog.writeFits`.
        """
        table = self.toTable()
        table.writeFits(*args)

    @classmethod
    def fromTable(cls, table):
        """Construct a `Defects` from the contents of a
        `~lsst.afw.table.BaseCatalog`.

        Parameters
        ----------
        table : `lsst.afw.table.BaseCatalog`
            Table with one row per defect.

        Returns
        -------
        defects : `Defects`
            A `Defects` list.
        """

        defectList = []

        schema = table.getSchema()

        # Check schema to see which definitions we have
        if "X" in schema and "Y" in schema and "R" in schema and "SHAPE" in schema:
            # This is a FITS region style table
            isFitsRegion = True

        elif "x0" in schema and "y0" in schema and "width" in schema and "height" in schema:
            # This is a classic LSST-style defect table
            isFitsRegion = False

        else:
            raise ValueError("Unsupported schema for defects extraction")

        for r in table:
            record = r.extract("*")

            if isFitsRegion:
                if record["SHAPE"] != "BOX":
                    log.warning("Defect lists can only be defined using BOX not %s",
                                record["SHAPE"])
                box = lsst.geom.Box2I.makeCenteredBox(lsst.geom.Point2D(record["X"], record["Y"]),
                                                      lsst.geom.Extent2I(record["R"]))

            elif "x0" in record and "y0" in record and "width" in record and "height" in record:
                # This is a classic LSST-style defect table
                box = lsst.geom.Box2I(lsst.geom.Point2I(record["x0"], record["y0"]),
                                      lsst.geom.Extent2I(record["width"], record["height"]))

            defectList.append(box)

        return Defects(defectList)

    @classmethod
    def readFits(cls, *args):
        """Read defect list from FITS table.

        Parameters
        ----------
        *args
            Arguments to be forwarded to
            `lsst.afw.table.BaseCatalog.writeFits`.

        Returns
        -------
        defects : `Defects`
            Defects read from a FITS table.
        """
        table = lsst.afw.table.BaseCatalog.readFits(*args)
        return cls.fromTable(table)
