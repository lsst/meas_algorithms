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

from typing import Iterable
from itertools import product
import numpy as np

import lsst.utils.tests
import lsst.geom
import lsst.afw.image as afwImage
from lsst.meas.algorithms.cloughTocher2DInterpolator import (
    CloughTocher2DInterpolateTask,
)
from lsst.meas.algorithms import CloughTocher2DInterpolatorUtils as ctUtils


class CloughTocher2DInterpolateTestCase(lsst.utils.tests.TestCase):
    """Test the CloughTocher2DInterpolateTask."""

    def setUp(self):
        super().setUp()

        self.maskedimage = afwImage.MaskedImageF(100, 121)
        for x in range(100):
            for y in range(121):
                self.maskedimage[x, y] = (3 * y + x * 5, 0, 1.0)

        # Clone the maskedimage so we can compare it after running the task.
        self.reference = self.maskedimage.clone()

        # Set some central pixels as SAT
        sliceX, sliceY = slice(30, 35), slice(40, 45)
        self.maskedimage.mask[sliceX, sliceY] = afwImage.Mask.getPlaneBitMask("SAT")
        self.maskedimage.image[sliceX, sliceY] = np.nan
        # Put nans here to make sure interp is done ok

        # Set an entire column as BAD
        self.maskedimage.mask[54:55, :] = afwImage.Mask.getPlaneBitMask("BAD")
        self.maskedimage.image[54:55, :] = np.nan

        # Set an entire row as BAD
        self.maskedimage.mask[:, 110:111] = afwImage.Mask.getPlaneBitMask("BAD")
        self.maskedimage.image[:, 110:111] = np.nan

        # Set an asymmetric region as BAD
        self.maskedimage.mask[41:42, 63:66] = afwImage.Mask.getPlaneBitMask("SAT")
        self.maskedimage.image[41:42, 63:66] = np.nan
        self.maskedimage.mask[42:43, 63:65] = afwImage.Mask.getPlaneBitMask("SAT")
        self.maskedimage.image[42:43, 63:65] = np.nan
        self.maskedimage.mask[44, 63] = afwImage.Mask.getPlaneBitMask("SAT")
        self.maskedimage.image[44, 63] = np.nan

        # Set a diagonal set of pixels as CR
        for i in range(74, 78):
            self.maskedimage.mask[i, i] = afwImage.Mask.getPlaneBitMask("CR")
            self.maskedimage.image[i, i] = np.nan

        # Set one of the edges as EDGE
        self.maskedimage.mask[0:1, :] = afwImage.Mask.getPlaneBitMask("EDGE")
        self.maskedimage.image[0:1, :] = np.nan

        # Set a smaller streak at the edge
        self.maskedimage.mask[25:28, 0:1] = afwImage.Mask.getPlaneBitMask("EDGE")
        self.maskedimage.image[25:28, 0:1] = np.nan

        # Update the reference image's mask alone, so we can compare them after
        # running the task.
        self.reference.mask.array[:, :] = self.maskedimage.mask.array

        # Create a noise image
        self.noise = self.maskedimage.clone()
        rng = np.random.Generator(np.random.MT19937(5))
        self.noise.image.array[:, :] = rng.normal(size=self.noise.image.array.shape)

    @lsst.utils.tests.methodParametersProduct(n_runs=(1, 2), flipXY=(False, True))
    def test_interpolation(self, n_runs: int, flipXY: bool):
        """Test that the interpolation is done correctly.

        Parameters
        ----------
        n_runs : `int`
            Number of times to run the task. Running the task more than once
            should have no effect.
        flipXY : `bool`
            Whether to set the flipXY config parameter to True.
        """
        config = CloughTocher2DInterpolateTask.ConfigClass()
        config.badMaskPlanes = (
            "BAD",
            "SAT",
            "CR",
            "EDGE",
        )
        config.fillValue = 0.5
        config.flipXY = flipXY
        task = CloughTocher2DInterpolateTask(config)
        for n in range(n_runs):
            task.run(self.maskedimage)

        # Assert that the mask and the variance planes remain unchanged.
        self.assertImagesEqual(self.maskedimage.variance, self.reference.variance)
        self.assertMasksEqual(self.maskedimage.mask, self.reference.mask)

        # Check that the long streak of bad pixels have been replaced with the
        # fillValue, but not the short streak.
        np.testing.assert_array_equal(self.maskedimage.image[0:1, :].array, config.fillValue)
        with self.assertRaises(AssertionError):
            np.testing.assert_array_equal(self.maskedimage.image[25:28, 0:1].array, config.fillValue)

        # Check that interpolated pixels are close to the reference (original),
        # and that none of them is still NaN.
        self.assertTrue(np.isfinite(self.maskedimage.image.array).all())
        self.assertImagesAlmostEqual(
            self.maskedimage.image[1:, :],
            self.reference.image[1:, :],
            rtol=1e-05,
            atol=1e-08,
        )

    @lsst.utils.tests.methodParametersProduct(
        pass_badpix=(True, False),
        pass_goodpix=(True, False),
        flipXY=(False, True),
    )
    def test_interpolation_with_noise(
        self, pass_badpix: bool = True, pass_goodpix: bool = True, flipXY: bool = False
    ):
        """Test that we can reuse the badpix and goodpix.

        Parameters
        ----------
        pass_badpix : `bool`
            Whether to pass the badpix to the task?
        pass_goodpix : `bool`
            Whether to pass the goodpix to the task?
        flipXY : `bool`
            Whether to set the flipXY config parameter to True.
        """

        config = CloughTocher2DInterpolateTask.ConfigClass()
        config.flipXY = flipXY
        config.badMaskPlanes = (
            "BAD",
            "SAT",
            "CR",
            "EDGE",
        )
        task = CloughTocher2DInterpolateTask(config)

        badpix, goodpix = task.run(self.noise)
        task.run(
            self.maskedimage,
            badpix=(badpix if pass_badpix else None),
            goodpix=(goodpix if pass_goodpix else None),
        )

        # Check that the long streak of bad pixels by the edge have been
        # replaced with fillValue, but not the short streak.
        np.testing.assert_array_equal(self.maskedimage.image[0:1, :].array, config.fillValue)
        with self.assertRaises(AssertionError):
            np.testing.assert_array_equal(self.maskedimage.image[25:28, 0:1].array, config.fillValue)

        # Check that interpolated pixels are close to the reference (original),
        # and that none of them is still NaN.
        self.assertTrue(np.isfinite(self.maskedimage.image.array).all())
        self.assertImagesAlmostEqual(
            self.maskedimage.image[1:, :],
            self.reference.image[1:, :],
            rtol=1e-05,
            atol=1e-08,
        )

    def test_interpolation_with_flipXY(self):
        """Test that the interpolation with both values for flipXY."""
        config = CloughTocher2DInterpolateTask.ConfigClass()
        config.badMaskPlanes = (
            "BAD",
            "SAT",
            "CR",
            "EDGE",
        )
        config.flipXY = True
        task = CloughTocher2DInterpolateTask(config)
        badpix_true, goodpix_true = task.run(self.maskedimage)

        config.flipXY = False
        task = CloughTocher2DInterpolateTask(config)
        badpix_false, goodpix_false = task.run(self.maskedimage)

        # Check that the locations of the bad and the good pixels, and the good
        # pixel values themselves are identical.
        np.testing.assert_array_equal(goodpix_false, goodpix_true)
        np.testing.assert_array_equal(badpix_false[:, :2], badpix_true[:, :2])

        # Check that the interpolated values at at least approximately equal.
        np.testing.assert_array_equal(badpix_false[:, 2], badpix_true[:, 2])


class CloughTocher2DInterpolatorUtilsTestCase(CloughTocher2DInterpolateTestCase):
    """Test the CloughTocher2DInterpolatorUtils."""

    @classmethod
    def find_good_pixels_around_bad_pixels(
        cls,
        image: afwImage.MaskedImage,
        maskPlanes: Iterable[str],
        *,
        max_window_extent: lsst.geom.Extent2I,
        badpix: set | None = None,
        goodpix: dict | None = None,
    ):
        """Find the location of bad pixels, and neighboring good pixels.

        Parameters
        ----------
        image : `~lsst.afw.image.MaskedImage`
            Image from which to find the bad and the good pixels.
        maskPlanes : `list` [`str`]
            List of mask planes to consider as bad pixels.
        max_window_extent : `lsst.geom.Extent2I`
            Maximum extent of the window around a bad pixel to consider when
            looking for good pixels.
        badpix : `list` [`tuple` [`int`, `int`]], optional
            A known list of bad pixels. If provided, the function does not look for
            any additional bad pixels, but it verifies that the provided
            coordinates correspond to bad pixels. If an input``badpix`` is not
            found to be bad as specified by ``maskPlanes``, an exception is raised.
        goodpix : `dict` [`tuple` [`int`, `int`], `float`], optional
            A known mapping of the coordinates of good pixels to their values, to
            which any newly found good pixels locations will be added, and the
            values (even for existing items) will be updated.

        Returns
        -------
        badpix : `list` [`tuple` [`int`, `int`]]
            The coordinates of the bad pixels. If ``badpix`` was provided as an
            input argument, the returned quantity is the same as the input.
        goodpix : `dict` [`tuple` [`int`, `int`], `float`]
            Updated mapping of the coordinates of good pixels to their values.

        Raises
        ------
        RuntimeError
            If a pixel passed in as ``goodpix`` is found to be bad as specified by
            ``maskPlanes``.
        ValueError
            If an input ``badpix`` is not found to be bad as specified by
            ``maskPlanes``.
        """

        bbox = image.getBBox()
        if badpix is None:
            iterator = product(range(bbox.minX, bbox.maxX + 1), range(bbox.minY, bbox.maxY + 1))
            badpix = set()
        else:
            iterator = badpix

        if goodpix is None:
            goodpix = {}

        for x, y in iterator:
            if image.mask[x, y] & afwImage.Mask.getPlaneBitMask(maskPlanes):
                if (x, y) in goodpix:
                    raise RuntimeError(
                        f"Pixel ({x}, {y}) is bad as specified by maskPlanes {maskPlanes} but "
                        "passed in as goodpix"
                    )
                badpix.add((x, y))
                window = lsst.geom.Box2I.makeCenteredBox(
                    center=lsst.geom.Point2D(x, y),  # center has to be a Point2D instance.
                    size=max_window_extent,
                )
                # Restrict to the bounding box of the image.
                window.clip(bbox)

                for xx, yy in product(
                    range(window.minX, window.maxX + 1),
                    range(window.minY, window.maxY + 1),
                ):
                    if not (image.mask[xx, yy] & afwImage.Mask.getPlaneBitMask(maskPlanes)):
                        goodpix[(xx, yy)] = image.image[xx, yy]
            elif (x, y) in badpix:
                # If (x, y) is in badpix, but did not get flagged as bad,
                # raise an exception.
                raise ValueError(f"Pixel ({x}, {y}) is not bad as specified by maskPlanes {maskPlanes}")

        return badpix, goodpix

    def test_parity(self, buffer=4):
        """Test that the C++ implementation gives the same results as the
        pure-Python implementation.

        Parameters
        ----------
        buffer : `int`, optional
            Same as the buffer parameter in `findGoodPixelsAroundBadPixels`.
        """
        bpix, gpix = ctUtils.findGoodPixelsAroundBadPixels(
            self.maskedimage, ["BAD", "SAT", "CR", "EDGE"], buffer=buffer
        )
        badpix, goodpix = self.find_good_pixels_around_bad_pixels(
            self.maskedimage,
            ["BAD", "SAT", "CR", "EDGE"],
            max_window_extent=lsst.geom.Extent2I(2 * buffer + 1, 2 * buffer + 1),
        )

        self.assertEqual(len(goodpix), gpix.shape[0])
        for row in gpix:
            x, y, val = int(row[0]), int(row[1]), row[2]
            self.assertEqual(goodpix[(x, y)], val)

        self.assertEqual(set(zip(bpix[:, 0], bpix[:, 1])), badpix)

    def test_findGoodPixelsAroundBadPixels(self):
        """Test the findGoodPixelsAroundBadPixels utility functino."""
        badpix, goodpix = ctUtils.findGoodPixelsAroundBadPixels(
            self.maskedimage,
            ["BAD", "SAT", "CR", "EDGE"],
            buffer=4,
        )

        # Check that badpix and goodpix have no overlaps
        badSet = set(zip(badpix[:, 0], badpix[:, 1]))
        goodSet = set(zip(goodpix[:, 0], goodpix[:, 1]))
        self.assertEqual(len(badSet & goodSet), 0)

        # buffer = 0 should give no goodpix, but same badpix
        badpix0, goodpix0 = ctUtils.findGoodPixelsAroundBadPixels(
            self.maskedimage,
            ["BAD", "SAT", "CR", "EDGE"],
            buffer=0,
        )

        self.assertEqual(len(goodpix0), 0)
        np.testing.assert_array_equal(badpix0, badpix)

        # For large enough buffer, badpix and goodpix should be mutually
        # exclusive and complete. This also checks that edges are handled.
        badpix, goodpix = ctUtils.findGoodPixelsAroundBadPixels(
            self.maskedimage,
            ["BAD", "SAT", "CR", "EDGE"],
            buffer=251,
        )

        self.assertEqual(
            len(badpix) + len(goodpix),
            self.maskedimage.getWidth() * self.maskedimage.getHeight(),
        )

    def test_update_functions(self):
        """Test updateArrayFromImage and updateImageFromArray behave as
        expected.
        """
        badpix, _ = ctUtils.findGoodPixelsAroundBadPixels(
            self.maskedimage,
            ["BAD", "SAT", "CR", "EDGE"],
            buffer=3,
        )

        # Ensure that maskedimage and reference are not the same initially.
        with self.assertRaises(AssertionError):
            self.assertImagesEqual(self.maskedimage.image, self.reference.image)

        # Update badpix values from the reference image
        ctUtils.updateArrayFromImage(badpix, self.reference.image)

        # Update maskedimage from badpix values
        ctUtils.updateImageFromArray(self.maskedimage.image, badpix)

        # maskedimage and reference image should now to be identifical
        self.assertImagesEqual(self.maskedimage.image, self.reference.image)

    @lsst.utils.tests.methodParametersProduct(x0=(0, 23, -53), y0=(0, 47, -31))
    def test_origin(self, x0=23, y0=47):
        """Test that we get consistent results with arbitrary image origins.

        Parameters
        ----------
        x0 : `int`
            The origin of the image along the horizontal axis.
        y0 : `int`
            The origin of the image along the vertical axis.
        """
        # Calling setUp explicitly becomes necessary, as we change in the image
        # in-place and need to reset to the original state when running with
        # a different set of parameters.
        self.setUp()
        badpix0, goodpix0 = ctUtils.findGoodPixelsAroundBadPixels(
            self.maskedimage,
            ["BAD", "SAT", "CR", "EDGE"],
            buffer=4,
        )

        # Check that badpix and goodpix have no overlaps
        badSet = set(zip(badpix0[:, 0], badpix0[:, 1]))
        goodSet = set(zip(goodpix0[:, 0], goodpix0[:, 1]))
        self.assertEqual(len(badSet & goodSet), 0)

        # Set a non-trivial xy0 for the maskedimage
        self.maskedimage.setXY0(lsst.geom.Point2I(x0, y0))
        badpix, goodpix = ctUtils.findGoodPixelsAroundBadPixels(
            self.maskedimage,
            ["BAD", "SAT", "CR", "EDGE"],
            buffer=4,
        )

        # Adjust the x and y columns with origin, so we can compare them.
        badpix0[:, 0] += x0
        goodpix0[:, 0] += x0
        badpix0[:, 1] += y0
        goodpix0[:, 1] += y0

        # The third column (pixel values) must match exactly if the
        # corresponding pixel values are read, regardless of the coordinate.
        np.testing.assert_array_equal(goodpix, goodpix0)
        np.testing.assert_array_equal(badpix, badpix0)

        # Update one of the goodpix arrays from image and check that it is
        # invariant. It would be invariant if it handles the pixel coordinates
        # consistently.
        ctUtils.updateArrayFromImage(goodpix0, self.maskedimage.image)
        np.testing.assert_array_equal(goodpix, goodpix0)

        # There should be some nan values right now.
        self.assertFalse(np.isfinite(self.maskedimage.image.array).all())

        # There should not be any nan values if the image is updated correctly.
        badpix[:, 2] = -99
        ctUtils.updateImageFromArray(self.maskedimage.image, badpix)
        self.assertTrue(np.isfinite(self.maskedimage.image.array).all())


def setup_module(module):
    lsst.utils.tests.init()


class MemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
