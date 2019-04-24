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
import itertools
import collections.abc
from deprecated.sphinx import deprecated
import numpy as np
import copy
import datetime

import lsst.geom
import lsst.pex.policy as policy
import lsst.afw.table
import lsst.afw.detection
import lsst.afw.image
from lsst.daf.base import PropertyList

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

    _OBSTYPE = "defects"
    """The calibration type used for ingest."""

    def __init__(self, defectList=None, metadata=None):
        self._defects = []

        if metadata is not None:
            self._metadata = metadata
        else:
            self.setMetadata()

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
        elif isinstance(value, lsst.geom.PointI):
            value = Defect(lsst.geom.Box2I(value, lsst.geom.Extent2I(1, 1)))
        elif isinstance(value, lsst.afw.image.DefectBase):
            value = Defect(value.getBBox())
        else:
            raise ValueError(f"Defects must be of type Defect, BoxI, or PointI, not '{value!r}'")
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
        """Compare if two `Defects` are equal.

        Two `Defects` are equal if their bounding boxes are equal and in
        the same order.  Metadata content is ignored.
        """
        if not isinstance(other, self.__class__):
            return False

        # Assume equal if bounding boxes are equal
        for d1, d2 in zip(self, other):
            if d1.getBBox() != d2.getBBox():
                return False

        return True

    def __str__(self):
        return "Defects(" + ",".join(str(d.getBBox()) for d in self) + ")"

    def insert(self, index, value):
        self._defects.insert(index, self._check_value(value))

    def getMetadata(self):
        """Retrieve metadata associated with these `Defects`.

        Returns
        -------
        meta : `lsst.daf.base.PropertyList`
            Metadata. The returned `~lsst.daf.base.PropertyList` can be
            modified by the caller and the changes will be written to
            external files.
        """
        return self._metadata

    def setMetadata(self, metadata=None):
        """Store a copy of the supplied metadata with the defects.

        Parameters
        ----------
        metadata : `lsst.daf.base.PropertyList`, optional
            Metadata to associate with the defects.  Will be copied and
            overwrite existing metadata.  If not supplied the existing
            metadata will be reset.
        """
        if metadata is None:
            self._metadata = PropertyList()
        else:
            self._metadata = copy.copy(metadata)

        # Ensure that we have the obs type required by calibration ingest
        self._metadata["OBSTYPE"] = self._OBSTYPE

    def copy(self):
        """Copy the defects to a new list, creating new defects from the
        bounding boxes.

        Returns
        -------
        new : `Defects`
            New list with new `Defect` entries.

        Notes
        -----
        This is not a shallow copy in that new `Defect` instances are
        created from the original bounding boxes.  It's also not a deep
        copy since the bounding boxes are not recreated.
        """
        return self.__class__(d.getBBox() for d in self)

    def transpose(self):
        """Make a transposed copy of this defect list.

        Returns
        -------
        retDefectList : `Defects`
            Transposed list of defects.
        """
        retDefectList = self.__class__()
        for defect in self:
            bbox = defect.getBBox()
            dimensions = bbox.getDimensions()
            nbbox = lsst.geom.Box2I(lsst.geom.Point2I(bbox.getMinY(), bbox.getMinX()),
                                    lsst.geom.Extent2I(dimensions[1], dimensions[0]))
            retDefectList.append(nbbox)
        return retDefectList

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

        # Set some metadata in the table (force OBSTYPE to exist)
        metadata = copy.copy(self.getMetadata())
        metadata["OBSTYPE"] = self._OBSTYPE
        table.setMetadata(metadata)

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

        # Add some additional headers useful for tracking purposes
        metadata = table.getMetadata()
        now = datetime.datetime.utcnow()
        metadata["DATE"] = now.isoformat()
        metadata["CALIB_CREATION_DATE"] = now.strftime("%Y-%m-%d")
        metadata["CALIB_CREATION_TIME"] = now.strftime("%T %Z").strip()

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

        return cls(defectList)

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
        defects = cls.fromTable(table)
        print(f"Meta: {table.getMetadata()!r}")
        defects.setMetadata(table.getMetadata())
        return defects

    @classmethod
    def readLsstDefectsFile(cls, filename):
        """Read defects information from an LSST format text file.

        Parameters
        ----------
        filename : `str`
            Name of text file containing the defect information.

        Returns
        -------
        defects : `Defects`
            The defects.

        Notes
        -----
        These defect text files are used as the human readable definitions
        of defects in calibration data definition repositories.  The format
        is to use four columns defined as follows:

        x0 : `int`
            X coordinate of bottom left corner of box.
        y0 : `int`
            Y coordinate of bottom left corner of box.
        width : `int`
            X extent of the box.
        height : `int`
            Y extent of the box.
        """
        # Use loadtxt so that ValueError is thrown if the file contains a
        # non-integer value. genfromtxt converts bad values to -1.
        defect_array = np.loadtxt(filename,
                                  dtype=[("x0", "int"), ("y0", "int"),
                                         ("x_extent", "int"), ("y_extent", "int")])

        return cls(lsst.geom.Box2I(lsst.geom.Point2I(row["x0"], row["y0"]),
                                   lsst.geom.Extent2I(row["x_extent"], row["y_extent"]))
                   for row in defect_array)

    @classmethod
    def fromFootprintList(cls, fpList):
        """Compute a defect list from a footprint list, optionally growing
        the footprints.

        Parameters
        ----------
        fpList : `list` of `lsst.afw.detection.Footprint`
            Footprint list to process.

        Returns
        -------
        defects : `Defects`
            List of defects.
        """
        return cls(itertools.chain.from_iterable(lsst.afw.detection.footprintToBBoxList(fp)
                                                 for fp in fpList))
