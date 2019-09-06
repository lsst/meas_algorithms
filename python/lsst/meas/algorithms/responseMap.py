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

__all__ = ["ResponseMap",]

from scipy.interpolate import interp1d
from astropy.table import Table

import lsst.afw.cameraGeom.utils as cg_utils



class ResponseMap(dict):
    def __init__(self, mode, *args, **kwargs):
        super().__inti__(*args, **kwargs)
        modes = ['DETECTOR', 'AMP', 'IMAGE']
        if mode not in modes:
            raise ValueError(f"Mode {mode} not supported.  Must be on of {modes}.")
        else:
            self.mode = mode

    def evaluate(self, detector, position, wavelength, kind='linear'):
        if self.mode == "DETECTOR":
            curve = self[detector.getName()]
            return curve.interpolate(wavelength, kind=kind)
        elif self.mode == "AMP":
            curves = self[detector.getName()]
            amp = cg_utils.findAmp(detector, Point2I(position))  # cast to Point2I if Point2D passed
            curve = curves[amp.getName()]
            return curve.interpolate(wavelength, kind=kind)
        elif self.mode == "IMAGE":
            raise NotImplementedError
        else:
            raise ValueError(f"Mode {mode} not supported.")
        




class Curve:
    def __init__(self, wavelength, response, metadata):
        self.wavelength = numpy.array(wavelength)
        self.response = numpy.array(response)
        self.metadata = metadata

    def interpolate(self, val, kind):
        return interp1d(self.wavelength, self.response, val, kind=kind)

    @classmethod
    def fromText(cls, filename):
        table = Table.read(filename, format='ascii.ecsv')
        return cls(table['wavelength'], table['response'], table.metadata)

    @classmethod
    def fromFits(self, filename):
        raise NotImplementedError
