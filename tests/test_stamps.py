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
import os

from lsst.meas.algorithms import stamps
from lsst.afw import image as afwImage
from lsst.afw.geom.testUtils import TransformTestBaseClass
from lsst.daf.base import PropertyList
import lsst.geom as geom
import lsst.afw.geom.transformFactory as tF
import lsst.utils.tests

np.random.seed(90210)


def make_stamps(n_stamps=3, use_archive=False):
    stampSize = 25
    # create dummy stamp stamps
    stampImages = [afwImage.MaskedImageF(stampSize, stampSize)
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
    archive_elements = [tF.makeTransform(geom.AffineTransform(np.random.rand(2))) if use_archive else None
                        for _ in range(n_stamps)]
    stamp_list = [stamps.Stamp(stamp_im=stampIm,
                               position=geom.SpherePoint(geom.Angle(ra, geom.degrees),
                                                         geom.Angle(dec, geom.degrees)),
                               archive_element=ae)
                  for stampIm, ra, dec, ae in zip(stampImages, ras, decs, archive_elements)]
    metadata = PropertyList()
    metadata['RA_DEG'] = ras
    metadata['DEC_DEG'] = decs

    return stamps.Stamps(stamp_list, metadata=metadata, use_archive=True)


class TransformTestClass(TransformTestBaseClass):
    """Class for unit tests for Transform<X>To<Y> within meas_algorithms.

    Inherits from `lsst.afw.geom.testUtils.TransformTestBaseClass`.
    """
    def getTestDir(self):
        """Returns the test directory of the `meas_algorithms` package.

        If a similar test is needed in another package, inherit from
        `TransformTestBaseClass` and overwrite this method; see the docstrings
        in the parent class.
        """
        return os.path.join(lsst.utils.getPackageDir("meas_algorithms"), "tests")


class StampsBaseTestCase(lsst.utils.tests.TestCase):
    """Test StampsBase.
    """
    def testReadFitsWithOptionsNotImplementedErrorRaised(self):
        """
        Test that subclasses have their own version
        of this implemented or an error is raised.
        """
        class FakeStampsBase(stamps.StampsBase):
            def __init__(self):
                return

        with self.assertRaises(NotImplementedError):
            FakeStampsBase.readFitsWithOptions('noFile', {})

    def testReadFitsWithOptionsMetadataError(self):
        """Test that error is raised when STAMPCLS returns None
        """
        with tempfile.NamedTemporaryFile() as f:
            ss = make_stamps()
            emptyMetadata = PropertyList()
            stamps.writeFits(
                f.name, [ss[0]], emptyMetadata, None, True, True
            )
            with self.assertRaises(RuntimeError):
                stamps.StampsBase.readFits(f.name)

    def testReadFitsReturnsNewClass(self):
        """Test that readFits will return subclass
        """
        class FakeStampsBase(stamps.StampsBase):
            def __init__(self):
                self._metadata = {}
                return

            @classmethod
            def readFitsWithOptions(cls, filename, options):
                return cls()

            def _refresh_metadata(self):
                self._metadata = {}

        fakeStamps = FakeStampsBase.readFitsWithOptions('noFile', {})
        self.assertEqual(type(fakeStamps), FakeStampsBase)


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
        """Test the class' write and readFits methods for the special case of
           one stamp.
        """
        self.roundtrip(make_stamps(1))

    def testIOsub(self):
        """Test the class' write and readFits when passing on a bounding box.
        """
        bbox = geom.Box2I(geom.Point2I(3, 9), geom.Extent2I(11, 7))
        ss = make_stamps()
        with tempfile.NamedTemporaryFile() as f:
            ss.writeFits(f.name)
            options = {'bbox': bbox}
            subStamps = stamps.Stamps.readFitsWithOptions(f.name, options)
            for s1, s2 in zip(ss, subStamps):
                self.assertEqual(bbox.getDimensions(), s2.stamp_im.getDimensions())
                self.assertMaskedImagesAlmostEqual(s1.stamp_im[bbox], s2.stamp_im)

    def testIOarchive(self):
        """Test the class' write and readFits when Stamps contain Persistables.
        """
        self.roundtripWithArchive(make_stamps(use_archive=True))

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

    def roundtripWithArchive(self, ss):
        """Round trip a Stamps object, including Archive elements, and check values
        """
        transformTest = TransformTestClass()
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
                transformTest.assertTransformsEqual(s1.archive_element, s2.archive_element)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
