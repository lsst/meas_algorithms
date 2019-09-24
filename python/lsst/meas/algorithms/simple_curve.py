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
from abc import ABC, abstractmethod
import datetime
import os
import numpy

import lsst.afw.cameraGeom.utils as cgUtils
from lsst.geom import Point2I


class Curve(ABC):
    mode = ''
    subclasses = dict()

    @abstractmethod
    def fromTable(self):
        pass

    @abstractmethod
    def toTable(self):
        pass

    @abstractmethod
    def evaluate(self, detector, position, wavelength, kind='linear'):
        pass

    @classmethod
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        if cls.mode in Curve.subclasses:
            raise ValueError(f'Class for mode, {cls.mode}, already defined')
        Curve.subclasses[cls.mode] = cls

    def interpolate(self, wavelengths, values, wavelength, kind):
        return interp1d(wavelengths, values, wavelength, kind=kind)

    @classmethod
    def readText(cls, filename):
        table = Table.read(filename, format='ascii.ecsv')
        return cls.subclasses[table.meta['MODE']].fromTable(table)

    @classmethod
    def readFits(cls, filename):
        table = Table.read(filename, format='fits')
        return cls.subclasses[table.meta['MODE']].fromTable(table)

    @staticmethod
    def _check_cols(cols, table):
        for col in cols:
            if col not in table.columns.keys():
                raise ValueError(f'The table must include a column named "{col}".')

    def _prep_write(self):
        now = datetime.datetime.utcnow()
        table = self.toTable()
        metadata = table.meta
        metadata["DATE"] = now.isoformat()
        metadata["CALIB_CREATION_DATE"] = now.strftime("%Y-%m-%d")
        metadata["CALIB_CREATION_TIME"] = now.strftime("%T %Z").strip()
        return table

    def writeText(self, filename):
        table = self._prep_write()
        # Force file extension to .ecsv
        path, ext = os.path.splitext(filename)
        filename = path + ".ecsv"
        table.write(filename, format="ascii.ecsv")
        return filename

    def writeFits(self, filename):
        table = self._prep_write()
        # Force file extension to .ecsv
        path, ext = os.path.splitext(filename)
        filename = path + ".fits"
        table.write(filename, format="fits")
        return filename


class DetectorCurve(Curve):
    mode = 'DETECTOR'

    def __init__(self, wavelength, efficiency, metadata):
        super().__init__()
        self.wavelength = numpy.array(wavelength)
        self.efficiency = numpy.array(efficiency)
        # make sure needed metadata is set if built directly from ctor.
        metadata.update({'MODE': 'DETECTOR', 'TYPE': 'QE'})
        self.metadata = metadata

    @classmethod
    def fromTable(cls, table):
        cls._check_cols(['wavelength', 'efficiency'], table)
        return cls(table['wavelength'].quantity, table['efficiency'].quantity, table.meta)

    def toTable(self):
        return Table({'wavelength': self.wavelength, 'efficiency': self.efficiency}, meta=self.metadata)

    def evaluate(self, detector, position, wavelength, kind='linear'):
        self.interpolate(self.wavelength, self.efficiency, wavelength, kind=kind)


class AmpCurve(Curve):
    mode = 'AMP'

    def __init__(self, amp_name_list, wavelength, efficiency, metadata):
        super().__init__()
        amp_names = set(amp_name_list.flat)
        self.data = {}
        for amp_name in amp_names:
            idx = numpy.where(amp_name_list == amp_name)[0]
            self.data[amp_name] = (wavelength[idx], efficiency[idx])
        # make sure needed metadata is set if built directly from ctor.
        metadata.update({'MODE': 'AMP', 'TYPE': 'QE'})
        self.metadata = metadata

    @classmethod
    def fromTable(cls, table):
        cls._check_cols(['amp_name', 'wavelength', 'efficiency'], table)
        return cls(table['amp_name'], table['wavelength'].quantity,
                   table['efficiency'].quantity, table.meta)

    def toTable(self):
        wavelength = None
        efficiency = None
        names = numpy.array([])
        for amp_name, val in self.data.items():
            # This will be preserve quantity
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
        return Table({'amp_name': names, 'wavelength': wavelength*wunit, 'efficiency': efficiency*eunit},
                     meta=self.metadata)

    def evaluate(self, detector, position, wavelength, kind='linear'):
        amp = cgUtils.findAmp(detector, Point2I(position))  # cast to Point2I if Point2D passed
        w, e = self.data[amp.getName()]
        return self.interpolate(w, e, wavelength, kind=kind)


class ImageCurve(Curve):
    mode = 'IMAGE'

    def fromTable(self, table):
        raise NotImplementedError()

    def toTable(self):
        raise NotImplementedError()

    def evaluate(self, detector, position, wavelength, kind='linear'):
        raise NotImplementedError()
