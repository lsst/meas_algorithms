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

from lsst.meas.algorithms import brightStarStamps
from lsst.afw import image as afwImage
from lsst.daf.base import PropertySet
import lsst.utils.tests


class BrightStarStampsTestCase(lsst.utils.tests.TestCase):
    """Test BrightStarStamps.
    """
    def setUp(self):
        np.random.seed(12)
        stampSize = (25, 25)
        # create dummy star stamps
        starImages = [afwImage.maskedImage.MaskedImageF(*stampSize)
                      for _ in range(3)]
        for starIm in starImages:
            starImArray = starIm.image.array
            starImArray += np.random.rand(*stampSize)
        ids = ["111", "aaa", "bbb"]
        mags = np.random.rand(3)
        # faint object to test magnitude cuts
        self.faintObjIdx = 1
        mags[self.faintObjIdx] = 18.
        ids[self.faintObjIdx] = "faint"
        fluxes = np.random.rand(3)
        self.starStamps = [brightStarStamps.BrightStarStamp(starStamp=starIm,
                                                            gaiaGMag=mag,
                                                            gaiaId=gaiaId,
                                                            annularFlux=flux)
                           for starIm, mag, gaiaId, flux in zip(starImages, mags, ids, fluxes)]
        self.innerRadius = 40
        self.outerRadius = 50
        self.bss = brightStarStamps.BrightStarStamps(self.starStamps, self.innerRadius, self.outerRadius)

    def tearDown(self):
        del self.bss
        del self.starStamps
        del self.innerRadius
        del self.outerRadius
        del self.faintObjIdx

    def testIO(self):
        """Test the class' write and readFits methods.

        The ``options`` argument to the read method is only used to read
        sub-BBoxes, which is handled by the Butler. Tests of this are done in
        afw.
        """
        with tempfile.NamedTemporaryFile() as f:
            self.bss.writeFits(f.name)
            options = PropertySet()
            bss2 = brightStarStamps.BrightStarStamps.readFitsWithOptions(f.name, options)
            self.assertEqual(len(self.bss), len(bss2))
            for mi1, mi2 in zip(self.bss.getMaskedImages(), bss2.getMaskedImages()):
                self.assertMaskedImagesAlmostEqual(mi1, mi2)
                np.testing.assert_almost_equal(self.bss.getMagnitudes(), bss2.getMagnitudes())
                np.testing.assert_almost_equal(self.bss.getAnnularFluxes(), bss2.getAnnularFluxes())
            for id1, id2 in zip(self.bss.getGaiaIds(), bss2.getGaiaIds()):
                self.assertEqual(id1, id2)

    def testMagnitudeSelection(self):
        brightOnly = self.bss.selectByMag(magMax=7)
        self.assertEqual(len(brightOnly), 2)
        self.assertFalse("faint" in brightOnly.getGaiaIds())

        faintOnly = self.bss.selectByMag(magMin=7)
        self.assertEqual(len(faintOnly), 1)
        self.assertEqual(faintOnly.getGaiaIds()[0], "faint")
        brightObj = self.bss[self.faintObjIdx]
        self.assertMaskedImagesAlmostEqual(brightObj.starStamp, faintOnly.getMaskedImages()[0])

    def testTypeMismatchHandling(self):
        fullStar = self.bss[0]
        # try passing on a dictionary and a maskedImage instead of a
        # BrightStarStamp
        falseStar = {"starStamp": fullStar.starStamp, "gaiaGMag": fullStar.gaiaGMag,
                     "gaiaId": fullStar.gaiaId, "annularFlux": fullStar.annularFlux}
        starIm = fullStar.starStamp
        for wrongType in [falseStar, starIm]:
            # test at initialization
            with self.assertRaises(ValueError):
                _ = brightStarStamps.BrightStarStamps([fullStar, wrongType], innerRadius=self.innerRadius,
                                                      outerRadius=self.outerRadius)
            # test at appending time
            with self.assertRaises(ValueError):
                self.bss.append(wrongType, innerRadius=self.innerRadius, outerRadius=self.outerRadius)

    def testAnnulusMismatch(self):
        """Test an exception is raised if mismatching annulus radii are
        given (as the annularFlux values would then be meaningless)."""
        metadata = self.bss.metadata
        # metadata contains annulus definition; brightStarStamps can be
        # instanciated without specifying annulus radii
        _ = brightStarStamps.BrightStarStamps(self.starStamps, metadata=metadata)
        # or by explicitely passing them on, in which case metadata is optional
        _ = brightStarStamps.BrightStarStamps(self.starStamps, innerRadius=self.innerRadius,
                                              outerRadius=self.outerRadius)
        # Both can be provided as long as values match
        _ = brightStarStamps.BrightStarStamps(self.starStamps, innerRadius=self.innerRadius,
                                              outerRadius=self.outerRadius,
                                              metadata=metadata)
        # An exception should be raised if they do not
        with self.assertRaises(AttributeError):
            _ = brightStarStamps.BrightStarStamps(self.starStamps, innerRadius=self.innerRadius/2,
                                                  metadata=metadata)
        with self.assertRaises(AttributeError):
            _ = brightStarStamps.BrightStarStamps(self.starStamps, outerRadius=self.outerRadius + 1,
                                                  metadata=metadata)
        # or if one tries to concatenate BrightStarStamps with different
        # annulus radii
        bss2 = brightStarStamps.BrightStarStamps(self.starStamps, innerRadius=self.innerRadius/2,
                                                 outerRadius=self.outerRadius)
        with self.assertRaises(AttributeError):
            self.bss.extend(bss2)
        # or append an extra stamp with different annulus radii
        fullStar = self.bss[0]
        with self.assertRaises(AttributeError):
            self.bss.append(fullStar, innerRadius=self.innerRadius/2, outerRadius=self.outerRadius + 1)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
