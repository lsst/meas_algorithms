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
import warnings
import scipy

import lsst.utils.tests
import lsst.afw.image
import lsst.afw.math
from lsst.meas.algorithms.maskStreaks import MaskStreaksTask, LineProfile, Line


def setDetectionMask(maskedImage, forceSlowBin=False, binning=None, detectedPlane="DETECTED",
                     badMaskPlanes=("NO_DATA", "INTRP", "BAD", "SAT", "EDGE"), detectionThreshold=5):
    """Make detection mask and set the mask plane.

    Create a binary image from a masked image by setting all data with signal-to-
    noise below some threshold to zero, and all data above the threshold to one.
    If the binning parameter has been set, this procedure will be preceded by a
    weighted binning of the data in order to smooth the result, after which the
    result is scaled back to the original dimensions. Set the detection mask
    plane with this binary image.

    Parameters
    ----------
    maskedImage : `lsst.afw.image.maskedImage`
        Image to be (optionally) binned and converted.
    forceSlowBin : `bool`, optional
        Force usage of slower binning method to check that the two methods
        give the same result.
    binning : `int`, optional
        Number of pixels by which to bin image.
    detectedPlane : `str`, optional
        Name of mask with pixels that were detected above threshold in image.
    badMaskPlanes : `tuple`, optional
        Names of masks with pixels that are rejected.
    detectionThreshold : `float`, optional
        Boundary in signal-to-noise between non-detections and detections for
        making a binary image from the original input image.
    """
    data = maskedImage.image.array
    weights = 1 / maskedImage.variance.array
    mask = maskedImage.mask

    detectionMask = ((mask.array & mask.getPlaneBitMask(detectedPlane)))
    badPixelMask = mask.getPlaneBitMask(badMaskPlanes)
    badMask = (mask.array & badPixelMask) > 0
    fitMask = detectionMask.astype(bool) & ~badMask

    fitData = np.copy(data)
    fitData[~fitMask] = 0
    fitWeights = weights
    fitWeights[~fitMask] = 0

    if binning:
        # Do weighted binning:
        ymax, xmax = fitData.shape
        if (ymax % binning == 0) and (xmax % binning == 0) and (not forceSlowBin):
            # Faster binning method
            binNumeratorReshape = (fitData * fitWeights).reshape(ymax // binning, binning,
                                                                 xmax // binning, binning)
            binDenominatorReshape = fitWeights.reshape(binNumeratorReshape.shape)
            binnedNumerator = binNumeratorReshape.sum(axis=3).sum(axis=1)
            binnedDenominator = binDenominatorReshape.sum(axis=3).sum(axis=1)
        else:
            # Slower binning method when (image shape mod binsize) != 0
            warnings.warn('Using slow binning method--consider choosing a binsize that evenly divides '
                          f'into the image size, so that {ymax} mod binning == 0 '
                          f'and {xmax} mod binning == 0', stacklevel=2)
            xarray = np.arange(xmax)
            yarray = np.arange(ymax)
            xmesh, ymesh = np.meshgrid(xarray, yarray)
            xbins = np.arange(0, xmax + binning, binning)
            ybins = np.arange(0, ymax + binning, binning)
            numerator = fitWeights * fitData
            binnedNumerator, *_ = scipy.stats.binned_statistic_2d(ymesh.ravel(), xmesh.ravel(),
                                                                  numerator.ravel(), statistic='sum',
                                                                  bins=(ybins, xbins))
            binnedDenominator, *_ = scipy.stats.binned_statistic_2d(ymesh.ravel(), xmesh.ravel(),
                                                                    fitWeights.ravel(), statistic='sum',
                                                                    bins=(ybins, xbins))
        binnedData = np.zeros(binnedNumerator.shape)
        ind = binnedDenominator != 0
        np.divide(binnedNumerator, binnedDenominator, out=binnedData, where=ind)
        binnedWeight = binnedDenominator
        binMask = (binnedData * binnedWeight**0.5) > detectionThreshold
        tmpOutputMask = binMask.repeat(binning, axis=0)[:ymax]
        outputMask = tmpOutputMask.repeat(binning, axis=1)[:, :xmax]
    else:
        outputMask = (fitData * fitWeights**0.5) > detectionThreshold

    # Clear existing Detected Plane:
    maskedImage.mask.array &= ~maskedImage.mask.getPlaneBitMask(detectedPlane)

    # Set Detected Plane with the binary detection mask:
    maskedImage.mask.array[outputMask] |= maskedImage.mask.getPlaneBitMask(detectedPlane)


class TestMaskStreaks(lsst.utils.tests.TestCase):

    def setUp(self):
        self.config = MaskStreaksTask.ConfigClass()
        self.config.dChi2Tolerance = 1e-6
        self.fst = MaskStreaksTask(config=self.config)

        self.testx = 500
        self.testy = 600
        self.exposure = lsst.afw.image.ExposureF(self.testy, self.testx)
        rand = lsst.afw.math.Random(seed=98765)
        lsst.afw.math.randomGaussianImage(self.exposure.image, rand)
        self.exposure.maskedImage.variance.set(1)
        self.maskName = "STREAK"
        self.detectedPlane = "DETECTED"

    def test_binning(self):
        """Test the two binning methods and the no-binning method"""

        binSize = 4
        self.assertEqual(self.testx % binSize, 0)
        self.assertEqual(self.testy % binSize, 0)

        testExposure1 = self.exposure.clone()
        setDetectionMask(testExposure1, binning=binSize, detectedPlane=self.detectedPlane)
        mask1 = testExposure1.getMask()
        reshapeBinning = mask1.array & mask1.getPlaneBitMask(self.detectedPlane)
        testExposure2 = self.exposure.clone()
        with self.assertWarns(Warning):
            setDetectionMask(testExposure2, binning=binSize, forceSlowBin=True)
        mask2 = testExposure2.getMask()
        scipyBinning = mask2.array & mask2.getPlaneBitMask(self.detectedPlane)
        self.assertAlmostEqual(reshapeBinning.tolist(), scipyBinning.tolist())

    def test_canny(self):
        """Test that Canny filter returns binary of equal shape"""

        zeroExposure = lsst.afw.image.ExposureF(self.testy, self.testx)
        cannyZeroExposure = self.fst._cannyFilter(zeroExposure.image.array)
        self.assertEqual(cannyZeroExposure.tolist(), zeroExposure.image.array.tolist())

        exposure = self.exposure.clone()
        setDetectionMask(exposure, detectedPlane=self.detectedPlane)
        mask = exposure.getMask()
        processedImage = mask.array & mask.getPlaneBitMask(self.detectedPlane)
        cannyNonZero = self.fst._cannyFilter(processedImage)
        self.assertEqual(cannyNonZero.tolist(), cannyNonZero.astype(bool).tolist())

    def test_runkht(self):
        """Test the whole thing"""

        # Empty image:
        zeroArray = np.zeros((self.testx, self.testy))
        zeroLines = self.fst._runKHT(zeroArray)
        self.assertEqual(len(zeroLines), 0)
        testExposure = self.exposure.clone()
        result = self.fst.run(testExposure)
        self.assertEqual(len(result.lines), 0)
        resultMask = testExposure.mask.array & testExposure.mask.getPlaneBitMask(self.maskName)
        self.assertEqual(resultMask.tolist(), zeroArray.tolist())

        # Make image with line and check that input line is recovered:
        testExposure1 = self.exposure.clone()
        inputRho = 150
        inputTheta = 45
        inputSigma = 3
        testLine = Line(inputRho, inputTheta, inputSigma)
        lineProfile = LineProfile(testExposure.image.array, testExposure.variance.array**-1,
                                  line=testLine)
        testData = lineProfile.makeProfile(testLine, fitFlux=False)
        testExposure1.image.array = testData * 100
        detectedInd = abs(testData) > 0.1 * (abs(testData)).max()

        testExposure1.mask.addMaskPlane(self.detectedPlane)
        testExposure1.mask.array[detectedInd] |= testExposure1.mask.getPlaneBitMask(self.detectedPlane)
        testExposure2 = testExposure1.clone()
        setDetectionMask(testExposure2, detectedPlane=self.detectedPlane)

        result_withSetDetMask = self.fst.run(testExposure2)
        self.assertEqual(len(result_withSetDetMask.lines), 1)
        self.assertAlmostEqual(inputRho, result_withSetDetMask.lines[0].rho, places=2)
        self.assertAlmostEqual(inputTheta, result_withSetDetMask.lines[0].theta, places=2)
        self.assertAlmostEqual(inputSigma, result_withSetDetMask.lines[0].sigma, places=2)

        result = self.fst.run(testExposure1)
        self.assertEqual(len(result.lines), 1)
        self.assertAlmostEqual(inputRho, result.lines[0].rho, places=2)
        self.assertAlmostEqual(inputTheta, result.lines[0].theta, places=2)
        self.assertAlmostEqual(inputSigma, result.lines[0].sigma, places=2)


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
