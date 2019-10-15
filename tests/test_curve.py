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

import os
import unittest
import numpy as np
from scipy import signal
import astropy.units as u

import lsst.meas.algorithms as algorithms
import lsst.utils.tests
from lsst.geom import Point2I, Point2D, Box2I, Extent2I
import lsst.afw.cameraGeom.utils as cgUtils

TESTDIR = os.path.abspath(os.path.dirname(__file__))


class MockAmp:
    def __init__(self, name, bbox):
        self.name = name
        self.box = bbox

    def getName(self):
        return self.name

    def getBBox(self):
        return self.box


class CurveTestCase(lsst.utils.tests.TestCase):
    """Tests for the Curve class"""

    def curve_tester(self, curve_class, args):
        curve = curve_class(*args)

        # Serialization round trip
        table = curve.toTable()
        curve2 = curve_class.fromTable(table)
        self.assertEqual(curve, curve2)

        # via FITS
        with lsst.utils.tests.getTempFilePath(".fits") as tmpFile:
            curve.writeFits(tmpFile)
            curve2 = algorithms.Curve.readFits(tmpFile)

        self.assertEqual(curve2, curve)

        # via text file
        with lsst.utils.tests.getTempFilePath(".ecsv") as tmpFile:
            curve.writeText(tmpFile)
            curve2 = algorithms.Curve.readText(tmpFile)

        self.assertEqual(curve2, curve)

        # Check bad values
        with self.assertRaises(ValueError):
            # test that raised when non quantities are passed
            nargs = []
            for arg in args:
                if hasattr(arg, 'unit'):
                    nargs.append(arg.value)
                else:
                    nargs.append(arg)
            _ = curve_class(*nargs)

            nargs = []
            for arg in args:
                if hasattr(arg, 'unit'):
                    nargs.append(arg._set_unit(None))
                else:
                    nargs.append(arg)
            _ = curve_class(*nargs)

    def interp_tester(self, curve_class, args, detector):
        curve = curve_class(*args)
        w = 3500*u.angstrom
        xs = np.linspace(0, 1023, 33)
        ys = np.linspace(0, 1023, 33)
        val_map = {'A': 0.9329662, 'B': 0.7463730}
        # Does interpolation work
        for x, y in zip(xs, ys):
            point = Point2D(x, y)
            if detector:
                amp = cgUtils.findAmp(detector, Point2I(point))
                value = val_map[amp.getName()]
            else:
                value = 0.9329662
            interp_val = curve.evaluate(detector, point, w)
            self.assertAlmostEqual(interp_val.value, value, places=5)
            self.assertEqual(interp_val.unit, u.percent)
        # Does interpolation work with arrays
        w_arr = np.linspace(320, 430, 70)*u.nm
        out_arr = curve.evaluate(detector, point, w_arr)
        self.assertEqual(len(w_arr), len(out_arr))
        # Does interpolation with different units work as expected
        point = Point2D(500., 500.)
        val1 = curve.evaluate(detector, point, w)
        new_w = w.to(u.mm)
        val2 = curve.evaluate(detector, point, new_w)
        self.assertEqual(val1.value, val2.value)
        # interpolation with non-quantity should raise
        with self.assertRaises(ValueError):
            interp_val = curve.evaluate(detector, point, w.value)
        # Does out of band interpolation do something reasonable
        with self.assertRaises(ValueError):
            w = 0.*u.angstrom
            interp_val = curve.evaluate(detector, point, w)

    def test_curve(self):
        wavelength = np.linspace(3000, 5000, 150)*u.angstrom
        efficiency = signal.gaussian(len(wavelength), std=100)*u.percent
        # Future versions of astropy will pass unit through concatenation
        amp_wavelength = np.concatenate([wavelength.value, wavelength.value])*u.angstrom  # Two amps
        amp_efficiency = np.concatenate([efficiency.value, efficiency.value*0.8])*u.percent  # Two amps
        amp_name = np.concatenate([['A' for el in wavelength], ['B' for el in wavelength]])
        metadata = dict([('MODE', 'AMP'), ('TYPE', 'QE'), ('CALIBDATE', '1970-01-01T00:00:00'),
                         ('INSTRUME', 'ts8'), ('OBSTYPE', 'qe_curve'), ('DETECTOR', 99),
                         ('DATE', '2019-09-27T22:15:13.518320'), ('CALIB_CREATION_DATE', '2019-09-27'),
                         ('CALIB_CREATION_TIME', '22:15:13')])
        amp_input = (amp_name, amp_wavelength, amp_efficiency, metadata)
        detector_input = (wavelength, efficiency, metadata)
        amplist = [MockAmp('A', Box2I(Point2I(0, 0), Extent2I(512, 1025))),
                   MockAmp('B', Box2I(Point2I(512, 10), Extent2I(512, 1024)))]
        for curve, args, detector in zip((algorithms.AmpCurve, algorithms.DetectorCurve),
                                         (amp_input, detector_input),
                                         (amplist, None)):
            self.curve_tester(curve, args)
            self.interp_tester(curve, args, detector)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
