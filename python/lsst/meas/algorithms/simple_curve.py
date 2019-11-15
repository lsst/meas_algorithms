#
# LSST Data Management System
#
# Copyright 2019  AURA/LSST.
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

__all__ = ["Curve", "AmpCurve", "DetectorCurve", "ImageCurve"]

from scipy.interpolate import interp1d
from astropy.table import Table
import astropy.units as u
from abc import ABC, abstractmethod
import datetime
import os
import numpy

import lsst.afw.cameraGeom.utils as cgUtils
from lsst.geom import Point2I


class Curve(ABC):
    """ An abstract class to represent an arbitrary curve with
    interpolation.
    """
    mode = ''
    subclasses = dict()

    def __init__(self, wavelength, efficiency, metadata):
        if not (isinstance(wavelength, u.Quantity) and wavelength.unit.physical_type == 'length'):
            raise ValueError('The wavelength must be a quantity with a length sense.')
        if not isinstance(efficiency, u.Quantity) or efficiency.unit != u.percent:
            raise ValueError('The efficiency must be a quantity with units of percent.')
        self.wavelength = wavelength
        self.efficiency = efficiency
        # make sure needed metadata is set if built directly from ctor.
        metadata.update({'MODE': self.mode, 'TYPE': 'QE'})
        self.metadata = metadata

    @classmethod
    @abstractmethod
    def fromTable(cls, table):
        """Class method for constructing a `Curve` object.

        Parameters
        ----------
        table : `astropy.table.Table`
            Table containing metadata and columns necessary
            for constructing a `Curve` object.

        Returns
        -------
        curve : `Curve`
            A `Curve` subclass of the appropriate type according
            to the table metadata
        """
        pass

    @abstractmethod
    def toTable(self):
        """Convert this `Curve` object to an `astropy.table.Table`.

        Parameters
        ----------
        Paramters vary by subclass implementation.

        Returns
        -------
        table : `astropy.table.Table`
            A table object containing the data from this `Curve`.
        """
        pass

    @abstractmethod
    def evaluate(self, detector, position, wavelength, kind='linear'):
        """Interpolate the curve at the specified position and wavelength.

        Parameters
        ----------
        detector : `lsst.afw.cameraGeom.Detector`
            Is used to find the appropriate curve given the position for
            curves that vary over the detector.  Ignored in the case where
            there is only a single curve per detector.
        position : `lsst.geom.Point2D`
            The position on the detector at which to evaluate the curve.
        wavelength : `astropy.units.Quantity`
            The wavelength(s) at which to make the interpolation.
        kind : `str`
            The type of interpolation to do.  See documentation for
            `scipy.interpolate.interp1d` for accepted values.

        Returns
        -------
        value : `astropy.units.Quantity`
            Interpolated value(s).  Number of values returned will match the
            length of `wavelength`.
        """
        pass

    @classmethod
    def __init_subclass__(cls, **kwargs):
        """Register subclasses with the abstract base class"""
        super().__init_subclass__(**kwargs)
        if cls.mode in Curve.subclasses:
            raise ValueError(f'Class for mode, {cls.mode}, already defined')
        Curve.subclasses[cls.mode] = cls

    @abstractmethod
    def __eq__(self, other):
        """Define equality for this class"""
        pass

    def compare_metadata(self, other,
                         keys_to_compare=['MODE', 'TYPE', 'CALIBDATE', 'INSTRUME', 'OBSTYPE', 'DETECTOR']):
        """Compare metadata in this object to another.

        Parameters
        ----------
        other : `Curve`
            The object with which to compare metadata.
        keys_to_compare : `list`
            List of metadata keys to compare.

        Returns
        -------
        same : `bool`
            Are the metadata the same?
        """
        for k in keys_to_compare:
            ret = ret and (self.metadata[k] == other.metadata[k])
        return ret

    def interpolate(self, wavelengths, values, wavelength, kind):
        """Interplate the curve at the specified wavelength(s).

        Parameters
        ----------
        wavelengths : `astropy.units.Quantity`
            The wavelength values for the curve.
        values : `astropy.units.Quantity`
            The y-values for the curve.
        wavelength : `astropy.units.Quantity`
            The wavelength(s) at which to make the interpolation.
        kind : `str`
            The type of interpolation to do.  See documentation for
            `scipy.interpolate.interp1d` for accepted values.

        Returns
        -------
        value : `astropy.units.Quantity`
            Interpolated value(s)
        """
        if not hasattr(wavelength, 'unit'):
            raise ValueError("Wavelengths at which to interpolate must be astropy quantities")
        interp_wavelength = wavelength.to(wavelengths.unit)
        f = interp1d(wavelengths, values, kind=kind)
        return f(interp_wavelength.value)*values.unit

    def getMetadata(self):
        """Return metadata

        Returns
        -------
        metadata : `dict`
            Dictionary of metadata for this curve.
        """
        # Needed to duck type as an object that can be ingested
        return self.metadata

    @classmethod
    def readText(cls, filename):
        """Class method for constructing a `Curve` object from
        the standardized text format.

        Parameters
        ----------
        filename : `str`
            Path to the text file to read.

        Returns
        -------
        curve : `Curve`
            A `Curve` subclass of the appropriate type according
            to the table metadata
        """
        table = Table.read(filename, format='ascii.ecsv')
        return cls.subclasses[table.meta['MODE']].fromTable(table)

    @classmethod
    def readFits(cls, filename):
        """Class method for constructing a `Curve` object from
        the standardized FITS format.

        Parameters
        ----------
        filename : `str`
            Path to the FITS file to read.

        Returns
        -------
        curve : `Curve`
            A `Curve` subclass of the appropriate type according
            to the table metadata
        """
        table = Table.read(filename, format='fits')
        return cls.subclasses[table.meta['MODE']].fromTable(table)

    @staticmethod
    def _check_cols(cols, table):
        """Check that the columns are in the table"""
        for col in cols:
            if col not in table.columns.keys():
                raise ValueError(f'The table must include a column named "{col}".')

    def _prep_write(self):
        """Compute standard metadata before writing file out"""
        now = datetime.datetime.utcnow()
        table = self.toTable()
        metadata = table.meta
        metadata["DATE"] = now.isoformat()
        metadata["CALIB_CREATION_DATE"] = now.strftime("%Y-%m-%d")
        metadata["CALIB_CREATION_TIME"] = now.strftime("%T %Z").strip()
        return table

    def writeText(self, filename):
        """ Write the `Curve` out to a text file.

        Parameters
        ----------
        filename : `str`
            Path to the text file to write.

        Returns
        -------
        filename : `str`
            Because this method forces a particular extension return
            the name of the file actually written.
        """
        table = self._prep_write()
        # Force file extension to .ecsv
        path, ext = os.path.splitext(filename)
        filename = path + ".ecsv"
        table.write(filename, format="ascii.ecsv")
        return filename

    def writeFits(self, filename):
        """ Write the `Curve` out to a FITS file.

        Parameters
        ----------
        filename : `str`
            Path to the FITS file to write.

        Returns
        -------
        filename : `str`
            Because this method forces a particular extension return
            the name of the file actually written.
        """
        table = self._prep_write()
        # Force file extension to .ecsv
        path, ext = os.path.splitext(filename)
        filename = path + ".fits"
        table.write(filename, format="fits")
        return filename


class DetectorCurve(Curve):
    """Subclass of `Curve` that represents a single curve per detector.

    Parameters
    ----------
    wavelength : `astropy.units.Quantity`
        Wavelength values for this curve
    efficiency : `astropy.units.Quantity`
        Quantum efficiency values for this curve
    metadata : `dict`
        Dictionary of metadata for this curve
    """
    mode = 'DETECTOR'

    def __eq__(self, other):
        return (self.compare_metadata(other) and
                numpy.array_equal(self.wavelength, other.wavelength) and
                numpy.array_equal(self.wavelength, other.wavelength))

    @classmethod
    def fromTable(cls, table):
        # Docstring inherited from base classs
        cls._check_cols(['wavelength', 'efficiency'], table)
        return cls(table['wavelength'].quantity, table['efficiency'].quantity, table.meta)

    def toTable(self):
        # Docstring inherited from base classs
        return Table({'wavelength': self.wavelength, 'efficiency': self.efficiency}, meta=self.metadata)

    def evaluate(self, detector, position, wavelength, kind='linear'):
        # Docstring inherited from base classs
        return self.interpolate(self.wavelength, self.efficiency, wavelength, kind=kind)


class AmpCurve(Curve):
    """Subclass of `Curve` that represents a curve per amp.

    Parameters
    ----------
    amp_name_list : iterable of `str`
        The name of the amp for each entry
    wavelength : `astropy.units.Quantity`
        Wavelength values for this curve
    efficiency : `astropy.units.Quantity`
        Quantum efficiency values for this curve
    metadata : `dict`
        Dictionary of metadata for this curve
    """
    mode = 'AMP'

    def __init__(self, amp_name_list, wavelength, efficiency, metadata):
        super().__init__(wavelength, efficiency, metadata)
        amp_names = set(amp_name_list)
        self.data = {}
        for amp_name in amp_names:
            idx = numpy.where(amp_name_list == amp_name)[0]
            # Deal with the case where the keys are bytes from FITS
            name = amp_name
            if isinstance(name, bytes):
                name = name.decode()
            self.data[name] = (wavelength[idx], efficiency[idx])

    def __eq__(self, other):
        ret = self.compare_metadata(other)
        for k in self.data.keys():
            ret = ret and numpy.array_equal(self.data[k][0], other.data[k][0])
            ret = ret and numpy.array_equal(self.data[k][1], other.data[k][1])
        return ret

    @classmethod
    def fromTable(cls, table):
        # Docstring inherited from base classs
        cls._check_cols(['amp_name', 'wavelength', 'efficiency'], table)
        return cls(table['amp_name'], table['wavelength'].quantity,
                   table['efficiency'].quantity, table.meta)

    def toTable(self):
        # Docstring inherited from base classs
        wavelength = None
        efficiency = None
        names = numpy.array([])
        # Loop over the amps and concatenate into three same length columns to feed
        # to the Table constructor.
        for amp_name, val in self.data.items():
            # This will preserve the quantity
            if wavelength is None:
                wunit = val[0].unit
                wavelength = val[0].value
            else:
                wavelength = numpy.concatenate([wavelength, val[0].value])
            if efficiency is None:
                eunit = val[1].unit
                efficiency = val[1].value
            else:
                efficiency = numpy.concatenate([efficiency, val[1].value])
            names = numpy.concatenate([names, numpy.full(val[0].shape, amp_name)])
        names = numpy.array(names)
        # Note that in future, the astropy.unit should make it through concatenation
        return Table({'amp_name': names, 'wavelength': wavelength*wunit, 'efficiency': efficiency*eunit},
                     meta=self.metadata)

    def evaluate(self, detector, position, wavelength, kind='linear'):
        # Docstring inherited from base classs
        amp = cgUtils.findAmp(detector, Point2I(position))  # cast to Point2I if Point2D passed
        w, e = self.data[amp.getName()]
        return self.interpolate(w, e, wavelength, kind=kind)


class ImageCurve(Curve):
    mode = 'IMAGE'

    def fromTable(self, table):
        # Docstring inherited from base classs
        raise NotImplementedError()

    def toTable(self):
        # Docstring inherited from base classs
        raise NotImplementedError()

    def evaluate(self, detector, position, wavelength, kind='linear'):
        # Docstring inherited from base classs
        raise NotImplementedError()
