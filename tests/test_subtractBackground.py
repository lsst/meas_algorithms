# This file is part of meas_algorithms.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This software is dual licensed under the GNU General Public License and also
# under a 3-clause BSD license. Recipients may choose which of these licenses
# to use; please see the files gpl-3.0.txt and/or bsd_license.txt,
# respectively.  If you choose the GPL option then the following text applies
# (but note that there is still no warranty even if you opt for BSD instead):
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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import unittest
import numpy as np

import lsst.geom
from lsst.meas.algorithms.testUtils import plantSources
import lsst.utils.tests

from lsst.meas.algorithms.subtractBackground import (
    SubtractBackgroundConfig,
    SubtractBackgroundTask,
    backgroundFlatContext,
)


class SubtractBackgroundTaskTestCase(lsst.utils.tests.TestCase):

    def _create_exposure(self, keepSky=True):
        """Return a simulated exposure with some background."""

        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(1024, 1024))
        kwid = 11
        self.sky = 2000

        coordList = [[10, 10, 0, 2.0]]

        exposure = plantSources(
            bbox=bbox,
            kwid=kwid,
            sky=self.sky,
            coordList=coordList,
            addPoissonNoise=True,
        )
        if keepSky:
            # Add sky level back in.
            exposure.image.array += self.sky

        return exposure

    def test_subtractBackground(self):
        exp = self._create_exposure()

        config = SubtractBackgroundConfig()
        task = SubtractBackgroundTask(config=config)

        results = task.run(exp)

        # Check the background was subtracted.
        self.assertFloatsAlmostEqual(np.mean(exp.image.array), 0.0, atol=0.1)

        # Check the background level.
        backgroundImage = results.background.getImage()
        self.assertFloatsAlmostEqual(np.mean(backgroundImage.array), self.sky, atol=0.1)

        # Do a second pass and confirm everything works.
        results2 = task.run(exp, background=results.background)
        self.assertEqual(len(results2.background), 2)

        self.assertFloatsAlmostEqual(np.mean(exp.image.array), 0.0, atol=0.1)
        backgroundImage2 = results2.background.getImage()
        self.assertFloatsAlmostEqual(np.mean(backgroundImage2.array), self.sky, atol=0.1)

    def test_backgroundFlatContext(self):
        exp = self._create_exposure(keepSky=False)

        # Check that doApply=False does nothing inside and out.
        maskedImage = exp.maskedImage.clone()
        with backgroundFlatContext(maskedImage, False):
            self.assertMaskedImagesAlmostEqual(maskedImage, exp.maskedImage)
        self.assertMaskedImagesAlmostEqual(maskedImage, exp.maskedImage)

        # Check that doApply=True does the conversions correctly.
        ratioImage = exp.image.clone()
        ratioImage.array[:, :] = 0.5
        maskedImage = exp.maskedImage.clone()
        with backgroundFlatContext(maskedImage, True, backgroundToPhotometricRatio=ratioImage):
            comparisonImage = exp.maskedImage.clone()
            comparisonImage /= ratioImage
            self.assertMaskedImagesAlmostEqual(maskedImage, comparisonImage)
        self.assertMaskedImagesAlmostEqual(maskedImage, exp.maskedImage)

        # Check that doApply=True with an incorrect ratio image raises.
        with self.assertRaises(ValueError):
            with backgroundFlatContext(maskedImage, True, backgroundToPhotometricRatio=exp.maskedImage):
                pass

        # Check that doApply=True with no ratio image raises.
        with self.assertRaises(RuntimeError):
            with backgroundFlatContext(maskedImage, True):
                pass


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
