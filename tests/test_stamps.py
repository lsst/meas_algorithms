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

import unittest
import numpy as np
import tempfile

from lsst.meas.algorithms import stamps
from lsst.afw import image as afwImage
from lsst.daf.base import PropertyList
import lsst.geom as geom
import lsst.utils.tests

np.random.seed(90210)


def make_stamps(n_stamps=3):
    stampSize = 25
    # create dummy stamp stamps
    stampImages = [afwImage.maskedImage.MaskedImageF(stampSize, stampSize)
                   for _ in range(n_stamps)]
    for stampIm in stampImages:
        stampImArray = stampIm.image.array
        stampImArray += np.random.rand(stampSize, stampSize)
        stampMaskArray = stampIm.mask.array
        stampMaskArray += 10
        stampVarArray = stampIm.variance.array
        stampVarArray += 1000.
    ras = np.random.rand(n_stamps)*360.
    decs = np.random.rand(n_stamps)*180 - 90
    stamp_list = [stamps.Stamp(stamp_im=stampIm,
                               position=geom.SpherePoint(geom.Angle(ra, geom.degrees),
                                                         geom.Angle(dec, geom.degrees)))
                  for stampIm, ra, dec in zip(stampImages, ras, decs)]
    metadata = PropertyList()
    metadata['RA_DEG'] = ras
    metadata['DEC_DEG'] = decs

    return stamps.Stamps(stamp_list, metadata=metadata)


class StampsTestCase(lsst.utils.tests.TestCase):
    """Test Stamps.
    """
    def testAppend(self):
        """Test ability to append to a Stamps object
        """
        ss = make_stamps()
        s = ss[-1]
        ss.append(s)
        self.roundtrip(ss)
        # check if appending something other than a Stamp raises
        with self.assertRaises(ValueError):
            ss.append('hello world')

    def testExtend(self):
        ss = make_stamps()
        ss2 = make_stamps()
        ss.extend([s for s in ss2])
        # check if extending with something other than a Stamps
        # object raises
        with self.assertRaises(ValueError):
            ss.extend(['hello', 'world'])

    def testIO(self):
        """Test the class' write and readFits methods.
        """
        self.roundtrip(make_stamps())

    def testIOone(self):
        """Test the class' write and readFits methods.
        """
        self.roundtrip(make_stamps(1))

    def roundtrip(self, ss):
        """Round trip a Stamps object to disk and check values
        """
        with tempfile.NamedTemporaryFile() as f:
            ss.writeFits(f.name)
            options = PropertyList()
            ss2 = stamps.Stamps.readFitsWithOptions(f.name, options)
            self.assertEqual(len(ss), len(ss2))
            for s1, s2 in zip(ss, ss2):
                self.assertMaskedImagesAlmostEqual(s1.stamp_im, s2.stamp_im)
                self.assertAlmostEqual(s1.position.getRa().asDegrees(),
                                       s2.position.getRa().asDegrees())
                self.assertAlmostEqual(s1.position.getDec().asDegrees(),
                                       s2.position.getDec().asDegrees())


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
