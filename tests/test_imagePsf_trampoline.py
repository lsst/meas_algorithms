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

import pickle
import unittest
from copy import deepcopy

import numpy as np

import lsst.utils.tests
from lsst.afw.image import Image, ExposureF
from lsst.afw.typehandling import StorableHelperFactory
from lsst.geom import Box2I, Point2I, Extent2I
from lsst.meas.algorithms import ImagePsf


class DummyImagePsf(ImagePsf):
    _factory = StorableHelperFactory(__name__, "DummyImagePsf")

    def __init__(self, image):
        ImagePsf.__init__(self)
        self.image = image

    # "public" virtual overrides
    def __deepcopy__(self, meta=None):
        return DummyImagePsf(self.image)

    def resized(self, width, height):
        raise NotImplementedError("resized not implemented for DummyImagePsf")

    def isPersistable(self):
        return True

    # "private" virtual overrides are underscored
    def _doComputeKernelImage(self, position=None, color=None):
        return self.image

    def _doComputeBBox(self, position=None, color=None):
        return self.image.getBBox()

    def _getPersistenceName(self):
        return "DummyImagePsf"

    def _getPythonModule(self):
        return __name__

    def _write(self):
        return pickle.dumps(self.image)

    @staticmethod
    def _read(pkl):
        return DummyImagePsf(pickle.loads(pkl))

    def __eq__(self, rhs):
        if isinstance(rhs, DummyImagePsf):
            return np.array_equal(self.image.array, rhs.image.array)
        return False


class ImagePsfTrampolineTestSuite(lsst.utils.tests.TestCase):
    def setUp(self):
        dimensions = Extent2I(7, 7)
        self.bbox = Box2I(Point2I(-dimensions/2), dimensions)
        self.img = Image(self.bbox, dtype=np.float64)
        x, y = np.ogrid[-3:4, -3:4]
        rsqr = x**2 + y**2
        # Some arbitrary circular double Gaussian
        self.img.array[:] = np.exp(-0.5*rsqr**2) + np.exp(-0.5*rsqr**2/4)
        self.img.array /= np.sum(self.img.array)
        self.psf = DummyImagePsf(self.img)

    def testImage(self):
        self.assertImagesEqual(
            self.img,
            self.psf.computeImage()
        )
        self.assertImagesEqual(
            self.img,
            self.psf.computeKernelImage()
        )

    def testBBox(self):
        self.assertEqual(
            self.bbox,
            self.psf.computeBBox()
        )

    def testResized(self):
        with self.assertRaises(NotImplementedError):
            self.psf.resized(9, 9)

    def testClone(self):
        clone1 = deepcopy(self.psf)
        clone2 = self.psf.clone()
        for clone in [clone1, clone2]:
            self.assertIsNot(clone, self.psf)
            self.assertImagesEqual(
                clone.computeImage(),
                self.psf.computeImage()
            )
            self.assertEqual(
                clone.computeApertureFlux(0.5),
                self.psf.computeApertureFlux(0.5)
            )
            self.assertEqual(
                clone.computeShape(),
                self.psf.computeShape()
            )

    def testPersistence(self):
        im = ExposureF(10, 10)
        im.setPsf(self.psf)
        self.assertEqual(im.getPsf(), self.psf)
        with lsst.utils.tests.getTempFilePath(".fits") as tmpFile:
            im.writeFits(tmpFile)
            newIm = ExposureF(tmpFile)
            self.assertEqual(newIm.getPsf(), im.getPsf())


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
