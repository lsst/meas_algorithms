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
from lsst.daf.base import PropertySet
import lsst.geom as geom
import lsst.utils.tests


class StampsTestCase(lsst.utils.tests.TestCase):
    """Test Stamps.
    """
    def setUp(self):
        np.random.seed(90210)
        stampSize = 25
        # create dummy stamp stamps
        stampImages = [afwImage.maskedImage.MaskedImageF(stampSize, stampSize)
                       for _ in range(3)]
        for stampIm in stampImages:
            stampImArray = stampIm.image.array
            stampImArray += np.random.rand(stampSize, stampSize)
        ras = np.random.rand(3)*360.
        decs = np.random.rand(3)*180 - 90
        sizes = [stampSize for _ in range(3)]
        self.stamps = [stamps.Stamp(stamp_im=stampIm,
                                    position=geom.SpherePoint(geom.Angle(ra, geom.degrees),
                                                              geom.Angle(dec, geom.degrees)),
                                    size=size)
                       for stampIm, ra, dec, size in zip(stampImages, ras, decs, sizes)]
        self.ss = stamps.Stamps(self.stamps)
        self.ss._metadata['RA_DEG'] = ras
        self.ss._metadata['DEC_DEG'] = decs
        self.ss._metadata['SIZE'] = sizes

    def tearDown(self):
        del self.ss
        del self.stamps

    def testIO(self):
        """Test the class' write and readFits methods.

        The ``options`` argument to the read method is only used to read
        sub-BBoxes, which is handled by the Butler. Tests of this are done in
        afw.
        """
        with tempfile.NamedTemporaryFile() as f:
            self.ss.writeFits(f.name)
            options = PropertySet()
            ss2 = stamps.Stamps.readFitsWithOptions(f.name, options)
            self.assertEqual(len(self.ss), len(ss2))
            for s1, s2 in zip(self.ss, ss2):
                self.assertMaskedImagesAlmostEqual(s1.stamp_im, s2.stamp_im)
                self.assertAlmostEqual(s1.position.getRa().asDegrees(),
                                       s2.position.getRa().asDegrees())
                self.assertAlmostEqual(s1.position.getDec().asDegrees(),
                                       s2.position.getDec().asDegrees())
                self.assertEqual(s1.size, s2.size)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
