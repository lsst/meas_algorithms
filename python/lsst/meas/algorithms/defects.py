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
import contextlib
import numpy as np
import copy
import datetime
import math
import numbers
import os.path
import warnings
import astropy.table

import lsst.geom
import lsst.afw.table
import lsst.afw.detection
import lsst.afw.image
import lsst.afw.geom
from lsst.daf.base import PropertyList

from . import Defect

log = logging.getLogger(__name__)

SCHEMA_NAME_KEY = "DEFECTS_SCHEMA"
SCHEMA_VERSION_KEY = "DEFECTS_SCHEMA_VERSION"


class Defects(collections.abc.MutableSequence):
    """Collection of `lsst.meas.algorithms.Defect`.

    Parameters
    ----------
    defectList : iterable of `lsst.meas.algorithms.Defect`
                 or `lsst.geom.BoxI`, optional
        Collections of defects to apply to the image.
    metadata : `lsst.daf.base.PropertyList`, optional
        Metadata to associate with the defects.  Will be copied and
        overwrite existing metadata.  If not supplied the existing
        metadata will be reset.
    normalize : `bool`
        If True, normalization is applied to the defects in ``defectList`` to
        remove duplicates, eliminate overlaps, etc.

    Notes
    -----
    Defects are stored within this collection in a "reduced" or "normalized"
    form: rather than simply storing the bounding boxes which are added to the
    collection, we eliminate overlaps and duplicates. This normalization
    procedure may introduce overhead when adding many new defects; it may be
    temporarily disabled using the `Defects.bulk_update` context manager if
    necessary.
    """

    _OBSTYPE = "defects"
    """The calibration type used for ingest."""

    def __init__(self, defectList=None, metadata=None, *, normalize=True):
        self._defects = []

        if defectList is not None:
            self._bulk_update = True
            for d in defectList:
                self.append(d)
        self._bulk_update = False

        if normalize:
            self._normalize()

        if metadata is not None:
            self._metadata = metadata
        else:
            self.setMetadata()

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
        self._normalize()

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

        # checking the bboxes with zip() only works if same length
        if len(self) != len(other):
            return False

        # Assume equal if bounding boxes are equal
        for d1, d2 in zip(self, other):
            if d1.getBBox() != d2.getBBox():
                return False

        return True

    def __str__(self):
        return "Defects(" + ",".join(str(d.getBBox()) for d in self) + ")"

    def _normalize(self):
        """Recalculate defect bounding boxes for efficiency.

        Notes
        -----
        Ideally, this would generate the provably-minimal set of bounding
        boxes necessary to represent the defects. At present, however, that
        doesn't happen: see DM-24781. In the cases of substantial overlaps or
        duplication, though, this will produce a much reduced set.
        """
        # In bulk-update mode, normalization is a no-op.
        if self._bulk_update:
            return

        # work out the minimum and maximum bounds from all defect regions.
        minX, minY, maxX, maxY = float('inf'), float('inf'), float('-inf'), float('-inf')
        for defect in self:
            bbox = defect.getBBox()
            minX = min(minX, bbox.getMinX())
            minY = min(minY, bbox.getMinY())
            maxX = max(maxX, bbox.getMaxX())
            maxY = max(maxY, bbox.getMaxY())

        region = lsst.geom.Box2I(lsst.geom.Point2I(minX, minY),
                                 lsst.geom.Point2I(maxX, maxY))

        mi = lsst.afw.image.MaskedImageF(region)
        self.maskPixels(mi, maskName="BAD")
        self._defects = Defects.fromMask(mi, "BAD")._defects

    @contextlib.contextmanager
    def bulk_update(self):
        """Temporarily suspend normalization of the defect list.
        """
        self._bulk_update = True
        try:
            yield
        finally:
            self._bulk_update = False
            self._normalize()

    def insert(self, index, value):
        self._defects.insert(index, self._check_value(value))
        self._normalize()

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

    def maskPixels(self, maskedImage, maskName="BAD"):
        """Set mask plane based on these defects.

        Parameters
        ----------
        maskedImage : `lsst.afw.image.MaskedImage`
            Image to process.  Only the mask plane is updated.
        maskName : str, optional
            Mask plane name to use.
        """
        # mask bad pixels
        mask = maskedImage.getMask()
        bitmask = mask.getPlaneBitMask(maskName)
        for defect in self:
            bbox = defect.getBBox()
            lsst.afw.geom.SpanSet(bbox).clippedTo(mask.getBBox()).setMask(mask, bitmask)

    def toFitsRegionTable(self):
        """Convert defect list to `~lsst.afw.table.BaseCatalog` using the
        FITS region standard.

        Returns
        -------
        table : `lsst.afw.table.BaseCatalog`
            Defects in tabular form.

        Notes
        -----
        The table created uses the
        `FITS regions <https://fits.gsfc.nasa.gov/registry/region.html>`_
        definition tabular format.  The ``X`` and ``Y`` coordinates are
        converted to FITS Physical coordinates that have origin pixel (1, 1)
        rather than the (0, 0) used in LSST software.
        """

        nrows = len(self._defects)

        schema = lsst.afw.table.Schema()
        x = schema.addField("X", type="D", units="pix", doc="X coordinate of center of shape")
        y = schema.addField("Y", type="D", units="pix", doc="Y coordinate of center of shape")
        shape = schema.addField("SHAPE", type="String", size=16, doc="Shape defined by these values")
        r = schema.addField("R", type="ArrayD", size=2, units="pix", doc="Extents")
        rotang = schema.addField("ROTANG", type="D", units="deg", doc="Rotation angle")
        component = schema.addField("COMPONENT", type="I", doc="Index of this region")
        table = lsst.afw.table.BaseCatalog(schema)
        table.resize(nrows)

        if nrows:
            # Adding entire columns is more efficient than adding
            # each element separately
            xCol = []
            yCol = []
            rCol = []

            for i, defect in enumerate(self._defects):
                box = defect.getBBox()
                center = box.getCenter()
                # Correct for the FITS 1-based offset
                xCol.append(center.getX() + 1.0)
                yCol.append(center.getY() + 1.0)

                width = box.width
                height = box.height

                if width == 1 and height == 1:
                    # Call this a point
                    shapeType = "POINT"
                else:
                    shapeType = "BOX"

                # Strings have to be added per row
                table[i][shape] = shapeType

                rCol.append(np.array([width, height], dtype=np.float64))

            # Assign the columns
            table[x] = np.array(xCol, dtype=np.float64)
            table[y] = np.array(yCol, dtype=np.float64)

            table[r] = np.array(rCol)
            table[rotang] = np.zeros(nrows, dtype=np.float64)
            table[component] = np.arange(nrows)

        # Set some metadata in the table (force OBSTYPE to exist)
        metadata = copy.copy(self.getMetadata())
        metadata["OBSTYPE"] = self._OBSTYPE
        metadata[SCHEMA_NAME_KEY] = "FITS Region"
        metadata[SCHEMA_VERSION_KEY] = 1
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
        table = self.toFitsRegionTable()

        # Add some additional headers useful for tracking purposes
        metadata = table.getMetadata()
        now = datetime.datetime.utcnow()
        metadata["DATE"] = now.isoformat()
        metadata["CALIB_CREATION_DATE"] = now.strftime("%Y-%m-%d")
        metadata["CALIB_CREATION_TIME"] = now.strftime("%T %Z").strip()

        table.writeFits(*args)

    def toSimpleTable(self):
        """Convert defects to a simple table form that we use to write
        to text files.

        Returns
        -------
        table : `lsst.afw.table.BaseCatalog`
            Defects in simple tabular form.

        Notes
        -----
        These defect tables are used as the human readable definitions
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
        schema = lsst.afw.table.Schema()
        x = schema.addField("x0", type="I", units="pix",
                            doc="X coordinate of bottom left corner of box")
        y = schema.addField("y0", type="I", units="pix",
                            doc="Y coordinate of bottom left corner of box")
        width = schema.addField("width", type="I", units="pix",
                                doc="X extent of box")
        height = schema.addField("height", type="I", units="pix",
                                 doc="Y extent of box")
        table = lsst.afw.table.BaseCatalog(schema)

        nrows = len(self._defects)
        table.resize(nrows)

        if nrows:

            xCol = []
            yCol = []
            widthCol = []
            heightCol = []

            for defect in self._defects:
                box = defect.getBBox()
                xCol.append(box.getBeginX())
                yCol.append(box.getBeginY())
                widthCol.append(box.getWidth())
                heightCol.append(box.getHeight())

            table[x] = np.array(xCol, dtype=np.int64)
            table[y] = np.array(yCol, dtype=np.int64)
            table[width] = np.array(widthCol, dtype=np.int64)
            table[height] = np.array(heightCol, dtype=np.int64)

        # Set some metadata in the table (force OBSTYPE to exist)
        metadata = copy.copy(self.getMetadata())
        metadata["OBSTYPE"] = self._OBSTYPE
        metadata[SCHEMA_NAME_KEY] = "Simple"
        metadata[SCHEMA_VERSION_KEY] = 1
        table.setMetadata(metadata)

        return table

    def writeText(self, filename):
        """Write the defects out to a text file with the specified name.

        Parameters
        ----------
        filename : `str`
            Name of the file to write.  The file extension ".ecsv" will
            always be used.

        Returns
        -------
        used : `str`
            The name of the file used to write the data (which may be
            different from the supplied name given the change to file
            extension).

        Notes
        -----
        The file is written to ECSV format and will include any metadata
        associated with the `Defects`.
        """

        # Using astropy table is the easiest way to serialize to ecsv
        afwTable = self.toSimpleTable()
        table = afwTable.asAstropy()

        metadata = afwTable.getMetadata()
        now = datetime.datetime.utcnow()
        metadata["DATE"] = now.isoformat()
        metadata["CALIB_CREATION_DATE"] = now.strftime("%Y-%m-%d")
        metadata["CALIB_CREATION_TIME"] = now.strftime("%T %Z").strip()

        table.meta = metadata.toDict()

        # Force file extension to .ecsv
        path, ext = os.path.splitext(filename)
        filename = path + ".ecsv"
        table.write(filename, format="ascii.ecsv")
        return filename

    @staticmethod
    def _get_values(values, n=1):
        """Retrieve N values from the supplied values.

        Parameters
        ----------
        values : `numbers.Number` or `list` or `np.array`
            Input values.
        n : `int`
            Number of values to retrieve.

        Returns
        -------
        vals : `list` or `np.array` or `numbers.Number`
            Single value from supplied list if ``n`` is 1, or `list`
            containing first ``n`` values from supplied values.

        Notes
        -----
        Some supplied tables have vectors in some columns that can also
        be scalars.  This method can be used to get the first number as
        a scalar or the first N items from a vector as a vector.
        """
        if n == 1:
            if isinstance(values, numbers.Number):
                return values
            else:
                return values[0]

        return values[:n]

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

        Notes
        -----
        Two table formats are recognized.  The first is the
        `FITS regions <https://fits.gsfc.nasa.gov/registry/region.html>`_
        definition tabular format written by `toFitsRegionTable` where the
        pixel origin is corrected from FITS 1-based to a 0-based origin.
        The second is the legacy defects format using columns ``x0``, ``y0``
        (bottom left hand pixel of box in 0-based coordinates), ``width``
        and ``height``.

        The FITS standard regions can only read BOX, POINT, or ROTBOX with
        a zero degree rotation.
        """

        defectList = []

        schema = table.getSchema()

        # Check schema to see which definitions we have
        if "X" in schema and "Y" in schema and "R" in schema and "SHAPE" in schema:
            # This is a FITS region style table
            isFitsRegion = True

            # Preselect the keys
            xKey = schema["X"].asKey()
            yKey = schema["Y"].asKey()
            shapeKey = schema["SHAPE"].asKey()
            rKey = schema["R"].asKey()
            rotangKey = schema["ROTANG"].asKey()

        elif "x0" in schema and "y0" in schema and "width" in schema and "height" in schema:
            # This is a classic LSST-style defect table
            isFitsRegion = False

            # Preselect the keys
            xKey = schema["x0"].asKey()
            yKey = schema["y0"].asKey()
            widthKey = schema["width"].asKey()
            heightKey = schema["height"].asKey()

        else:
            raise ValueError("Unsupported schema for defects extraction")

        for record in table:

            if isFitsRegion:
                # Coordinates can be arrays (some shapes in the standard
                # require this)
                # Correct for FITS 1-based origin
                xcen = cls._get_values(record[xKey]) - 1.0
                ycen = cls._get_values(record[yKey]) - 1.0
                shape = record[shapeKey].upper()
                if shape == "BOX":
                    box = lsst.geom.Box2I.makeCenteredBox(lsst.geom.Point2D(xcen, ycen),
                                                          lsst.geom.Extent2I(cls._get_values(record[rKey],
                                                                                             n=2)))
                elif shape == "POINT":
                    # Handle the case where we have an externally created
                    # FITS file.
                    box = lsst.geom.Point2I(xcen, ycen)
                elif shape == "ROTBOX":
                    # Astropy regions always writes ROTBOX
                    rotang = cls._get_values(record[rotangKey])
                    # We can support 0 or 90 deg
                    if math.isclose(rotang % 90.0, 0.0):
                        # Two values required
                        r = cls._get_values(record[rKey], n=2)
                        if math.isclose(rotang % 180.0, 0.0):
                            width = r[0]
                            height = r[1]
                        else:
                            width = r[1]
                            height = r[0]
                        box = lsst.geom.Box2I.makeCenteredBox(lsst.geom.Point2D(xcen, ycen),
                                                              lsst.geom.Extent2I(width, height))
                    else:
                        log.warning("Defect can not be defined using ROTBOX with non-aligned rotation angle")
                        continue
                else:
                    log.warning("Defect lists can only be defined using BOX or POINT not %s", shape)
                    continue

            else:
                # This is a classic LSST-style defect table
                box = lsst.geom.Box2I(lsst.geom.Point2I(record[xKey], record[yKey]),
                                      lsst.geom.Extent2I(record[widthKey], record[heightKey]))

            defectList.append(box)

        defects = cls(defectList)
        defects.setMetadata(table.getMetadata())

        # Once read, the schema headers are irrelevant
        metadata = defects.getMetadata()
        for k in (SCHEMA_NAME_KEY, SCHEMA_VERSION_KEY):
            if k in metadata:
                del metadata[k]

        return defects

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

    @classmethod
    def readText(cls, filename):
        """Read defect list from standard format text table file.

        Parameters
        ----------
        filename : `str`
            Name of the file containing the defects definitions.

        Returns
        -------
        defects : `Defects`
            Defects read from a FITS table.
        """
        with warnings.catch_warnings():
            # Squash warnings due to astropy failure to close files; we think
            # this is a real problem, but the warnings are even worse.
            # https://github.com/astropy/astropy/issues/8673
            warnings.filterwarnings("ignore", category=ResourceWarning, module="astropy.io.ascii")
            table = astropy.table.Table.read(filename)

        # Need to convert the Astropy table to afw table
        schema = lsst.afw.table.Schema()
        for colName in table.columns:
            schema.addField(colName, units=str(table[colName].unit),
                            type=table[colName].dtype.type)

        # Create AFW table that is required by fromTable()
        afwTable = lsst.afw.table.BaseCatalog(schema)

        afwTable.resize(len(table))
        for colName in table.columns:
            # String columns will fail -- currently we do not expect any
            afwTable[colName] = table[colName]

        # Copy in the metadata from the astropy table
        metadata = PropertyList()
        for k, v in table.meta.items():
            metadata[k] = v
        afwTable.setMetadata(metadata)

        # Extract defect information from the table itself
        return cls.fromTable(afwTable)

    @classmethod
    def readLsstDefectsFile(cls, filename):
        """Read defects information from a legacy LSST format text file.

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

        Files of this format were used historically to represent defects
        in simple text form.  Use `Defects.readText` and `Defects.writeText`
        to use the more modern format.
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
        # normalize is set to False to avoid recursively calling
        # fromMask/fromFootprintList in Defects.__init__.
        return cls(itertools.chain.from_iterable(lsst.afw.detection.footprintToBBoxList(fp)
                                                 for fp in fpList), normalize=False)

    @classmethod
    def fromMask(cls, maskedImage, maskName):
        """Compute a defect list from a specified mask plane.

        Parameters
        ----------
        maskedImage : `lsst.afw.image.MaskedImage`
            Image to process.
        maskName : `str` or `list`
            Mask plane name, or list of names to convert.

        Returns
        -------
        defects : `Defects`
            Defect list constructed from masked pixels.
        """
        mask = maskedImage.getMask()
        thresh = lsst.afw.detection.Threshold(mask.getPlaneBitMask(maskName),
                                              lsst.afw.detection.Threshold.BITMASK)
        fpList = lsst.afw.detection.FootprintSet(mask, thresh).getFootprints()
        return cls.fromFootprintList(fpList)
