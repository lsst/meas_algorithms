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
from tempfile import NamedTemporaryFile

import lsst.utils.tests
import numpy as np
from lsst.afw.image import MaskedImageF
from lsst.daf.base import PropertyList
from lsst.meas.algorithms import BrightStarStamp, BrightStarStamps
from lsst.meas.algorithms.stamps import StampsBase


class BrightStarStampsTestCase(lsst.utils.tests.TestCase):
    """Test BrightStarStamps."""

    def setUp(self):
        rng = np.random.Generator(np.random.MT19937(seed=5))
        stamp_size = (25, 25)

        # Generate simulated bright star stamps
        brightStarStamps = []
        self.metadatas = []
        for i in range(3):
            stamp = MaskedImageF(*stamp_size)
            stamp_array = stamp.image.array
            stamp_array += rng.random(stamp_size)
            psf = None
            wcs = None
            metadata = PropertyList()
            metadata.set("VISIT", i)
            metadata.set("DETECTOR", i + 1)
            metadata.set("REF_ID", f"ref{i}")
            metadata.set("REF_MAG", float(i * 5))
            metadata.set("POSITION_X", rng.random())
            metadata.set("POSITION_Y", rng.random())
            metadata.set("FOCAL_PLANE_RADIUS", rng.random())
            metadata.set("FOCAL_PLANE_ANGLE_DEGREES", rng.random())
            metadata.set("SCALE", rng.random())
            metadata.set("SCALE_ERR", rng.random())
            metadata.set("PEDESTAL", rng.random())
            metadata.set("PEDESTAL_ERR", rng.random())
            metadata.set("PEDESTAL_SCALE_COV", rng.random())
            metadata.set("GRADIENT_X", rng.random())
            metadata.set("GRADIENT_Y", rng.random())
            metadata.set("GLOBAL_REDUCED_CHI_SQUARED", rng.random())
            metadata.set("GLOBAL_DEGREES_OF_FREEDOM", rng.integers(1, 100))
            metadata.set("PSF_REDUCED_CHI_SQUARED", rng.random())
            metadata.set("PSF_DEGREES_OF_FREEDOM", rng.integers(1, 100))
            metadata.set("PSF_MASKED_FLUX_FRACTION", rng.random())
            self.metadatas.append(metadata)
            brightStarStamp = BrightStarStamp.factory(stamp, psf, wcs, metadata)
            brightStarStamps.append(brightStarStamp)
        self.primary_metadata = PropertyList()
        self.primary_metadata.set("TEST_KEY", "TEST VALUE")
        self.brightStarStamps = BrightStarStamps(brightStarStamps, self.primary_metadata)

    def tearDown(self):
        del self.brightStarStamps
        del self.metadatas
        del self.primary_metadata

    def testBrightStarStamps(self):
        """Test that BrightStarStamps can be serialized and deserialized."""

        with NamedTemporaryFile() as file:
            self.brightStarStamps.writeFits(file.name)
            # Test StampsBase correctly bounces to BrightStarStamps readFits
            brightStarStamps = StampsBase.readFits(file.name)

        primary_metadata = brightStarStamps.metadata
        self.assertEqual(self.primary_metadata["TEST_KEY"], primary_metadata["TEST_KEY"])
        self.assertEqual(len(self.metadatas), primary_metadata["N_STAMPS"])
        for input_metadata, input_stamp, output_stamp in zip(
            self.metadatas, self.brightStarStamps, brightStarStamps
        ):
            output_metadata = output_stamp.metadata
            for key in output_metadata.names():
                input_value = input_metadata[key]
                output_value = output_metadata[key]
                if isinstance(input_value, float):
                    self.assertAlmostEqual(input_value, output_value, places=10)
                else:
                    self.assertEqual(input_value, output_value)
            self.assertMaskedImagesAlmostEqual(input_stamp.stamp_im, output_stamp.stamp_im)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
