import unittest

import numpy as np
import numpy.testing as testing

import lsst.utils.tests
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
from lsst.meas.algorithms import AccumulatorMeanStack


class AccumulatorMeanStackTestCase(lsst.utils.tests.TestCase):
    def make_test_images_to_coadd(self):
        """Make test images to coadd."""
        rng = np.random.Generator(np.random.MT19937(5))

        # Choose an arbitrary number of images and image size.
        # Below we set bad pixels for either 1 or half these images,
        # at arbitrary locations.
        n_image = 10
        image_size = (200, 200)

        # Create noise field images each with constant variance.
        # Create some fake pixel-sized sources at fixed positions as well.
        # Create some bad pixels at random positions, and some others at
        # fixed positions (but not in all images).
        exposures = []
        variances = rng.uniform(size=n_image, low=5.0, high=10.0)

        for i in range(n_image):
            exposure = afwImage.ExposureF(width=image_size[0], height=image_size[1])
            exposure.variance.array[:, :] = variances[i]
            exposure.image.array[:, :] = rng.normal(loc=0.0,
                                                    scale=np.sqrt(variances[i]),
                                                    size=image_size)
            exposure.image.array[50, 50] = rng.normal(loc=100.0,
                                                      scale=np.sqrt(100.0))
            exposure.image.array[100, 100] = rng.normal(loc=200.0,
                                                        scale=np.sqrt(200.0))
            if (i == 2):
                # Create one outlier pixel
                exposure.image.array[150, 125] = 1000.0

            exposure.image.array[150, 150] = np.nan
            exposure.mask.array[150, 150] = afwImage.Mask.getPlaneBitMask(['NO_DATA'])

            # Create two saturated pixels, one above and one below the threshold.
            if (i < 5):
                exposure.mask.array[50, 100] = afwImage.Mask.getPlaneBitMask(['SAT'])
            if i == 1:
                exposure.mask.array[51, 100] = afwImage.Mask.getPlaneBitMask(['SAT'])

            exposure.mask.array[0: 2, :] |= afwImage.Mask.getPlaneBitMask(['EDGE'])

            exposures.append(exposure)

        weights = rng.uniform(size=n_image, low=0.9, high=1.1)

        return exposures, weights

    def make_coadd_exposure(self, exposure_ref):
        """Make a coadd exposure with mask planes."""
        coadd_exposure = afwImage.ExposureF(width=exposure_ref.getWidth(),
                                            height=exposure_ref.getHeight())
        coadd_exposure.mask.addMaskPlane("REJECTED")
        coadd_exposure.mask.addMaskPlane("CLIPPED")
        coadd_exposure.mask.addMaskPlane("SENSOR_EDGE")

        return coadd_exposure

    def make_stats_ctrl(self):
        """Make the reference stats_ctrl."""
        stats_ctrl = afwMath.StatisticsControl()
        stats_ctrl.setAndMask(afwImage.Mask.getPlaneBitMask(["NO_DATA",
                                                             "BAD",
                                                             "SAT",
                                                             "SUSPECT"]))
        stats_ctrl.setNanSafe(True)
        stats_ctrl.setWeighted(True)
        stats_ctrl.setCalcErrorFromInputVariance(True)
        bit = afwImage.Mask.getMaskPlane("SAT")
        stats_ctrl.setMaskPropagationThreshold(bit, 0.2)

        return stats_ctrl

    def make_mask_map(self, stats_ctrl):
        """Make the mask_map."""
        clipped = afwImage.Mask.getPlaneBitMask("CLIPPED")
        edge = afwImage.Mask.getPlaneBitMask("EDGE")
        no_data = afwImage.Mask.getPlaneBitMask("NO_DATA")
        to_reject = stats_ctrl.getAndMask() & (~no_data) & (~edge) & (~clipped)
        mask_map = [(to_reject, afwImage.Mask.getPlaneBitMask("REJECTED")),
                    (edge, afwImage.Mask.getPlaneBitMask("SENSOR_EDGE")),
                    (clipped, clipped)]

        return mask_map

    def test_online_coadd_input_variance_true(self):
        """Test online coaddition with calcErrorFromInputVariance=True."""
        exposures, weights = self.make_test_images_to_coadd()
        coadd_exposure = self.make_coadd_exposure(exposures[0])
        stats_ctrl = self.make_stats_ctrl()
        stats_ctrl.setCalcErrorFromInputVariance(True)
        mask_map = self.make_mask_map(stats_ctrl)

        stats_flags = afwMath.stringToStatisticsProperty("MEAN")
        clipped = afwImage.Mask.getPlaneBitMask("CLIPPED")

        masked_image_list = [exp.maskedImage for exp in exposures]

        afw_masked_image = afwMath.statisticsStack(masked_image_list,
                                                   stats_flags,
                                                   stats_ctrl,
                                                   weights,
                                                   clipped,
                                                   mask_map)

        mask_threshold_dict = AccumulatorMeanStack.stats_ctrl_to_threshold_dict(stats_ctrl)

        # Make the stack with the online accumulator
        stacker = AccumulatorMeanStack(
            coadd_exposure.image.array.shape,
            stats_ctrl.getAndMask(),
            mask_threshold_dict=mask_threshold_dict,
            mask_map=mask_map,
            no_good_pixels_mask=stats_ctrl.getNoGoodPixelsMask(),
            calc_error_from_input_variance=stats_ctrl.getCalcErrorFromInputVariance(),
            compute_n_image=True)

        for exposure, weight in zip(exposures, weights):
            stacker.add_masked_image(exposure.maskedImage, weight=weight)

        stacker.fill_stacked_masked_image(coadd_exposure.maskedImage)

        online_masked_image = coadd_exposure.maskedImage

        # The coadds match at the <1e-5 level.
        testing.assert_array_almost_equal(online_masked_image.image.array,
                                          afw_masked_image.image.array, decimal=5)
        testing.assert_array_almost_equal(online_masked_image.variance.array,
                                          afw_masked_image.variance.array)
        testing.assert_array_equal(online_masked_image.mask.array,
                                   afw_masked_image.mask.array)

    def test_online_coadd_input_variance_false(self):
        """Test online coaddition with calcErrorFromInputVariance=False."""
        exposures, weights = self.make_test_images_to_coadd()
        coadd_exposure = self.make_coadd_exposure(exposures[0])
        stats_ctrl = self.make_stats_ctrl()
        stats_ctrl.setCalcErrorFromInputVariance(False)
        mask_map = self.make_mask_map(stats_ctrl)

        stats_flags = afwMath.stringToStatisticsProperty("MEAN")
        clipped = afwImage.Mask.getPlaneBitMask("CLIPPED")

        masked_image_list = [exp.maskedImage for exp in exposures]

        afw_masked_image = afwMath.statisticsStack(masked_image_list,
                                                   stats_flags,
                                                   stats_ctrl,
                                                   weights,
                                                   clipped,
                                                   mask_map)

        mask_threshold_dict = AccumulatorMeanStack.stats_ctrl_to_threshold_dict(stats_ctrl)

        # Make the stack with the online accumulator

        # By setting no_good_pixels=None we have the same behavior
        # as the default from stats_ctrl.getNoGoodPixelsMask(), but
        # covers the alternate code path.
        stacker = AccumulatorMeanStack(
            coadd_exposure.image.array.shape,
            stats_ctrl.getAndMask(),
            mask_threshold_dict=mask_threshold_dict,
            mask_map=mask_map,
            no_good_pixels_mask=None,
            calc_error_from_input_variance=stats_ctrl.getCalcErrorFromInputVariance(),
            compute_n_image=True)

        for exposure, weight in zip(exposures, weights):
            stacker.add_masked_image(exposure.maskedImage, weight=weight)

        stacker.fill_stacked_masked_image(coadd_exposure.maskedImage)

        online_masked_image = coadd_exposure.maskedImage

        # The coadds match at the <1e-5 level.
        testing.assert_array_almost_equal(online_masked_image.image.array,
                                          afw_masked_image.image.array, decimal=5)
        # The computed variances match at the <1e-4 level.
        testing.assert_array_almost_equal(online_masked_image.variance.array,
                                          afw_masked_image.variance.array, decimal=4)
        testing.assert_array_equal(online_masked_image.mask.array,
                                   afw_masked_image.mask.array)

    def test_online_coadd_image(self):
        """Test online coaddition with regular non-masked images."""
        exposures, weights = self.make_test_images_to_coadd()
        coadd_exposure = self.make_coadd_exposure(exposures[0])
        stats_ctrl = self.make_stats_ctrl()
        stats_ctrl.setAndMask(0)
        stats_ctrl.setCalcErrorFromInputVariance(True)
        mask_map = self.make_mask_map(stats_ctrl)

        stats_flags = afwMath.stringToStatisticsProperty("MEAN")
        clipped = afwImage.Mask.getPlaneBitMask("CLIPPED")

        masked_image_list = [exp.maskedImage for exp in exposures]

        afw_masked_image = afwMath.statisticsStack(masked_image_list,
                                                   stats_flags,
                                                   stats_ctrl,
                                                   weights,
                                                   clipped,
                                                   mask_map)

        # Make the stack with the online accumulator
        stacker = AccumulatorMeanStack(
            coadd_exposure.image.array.shape,
            stats_ctrl.getAndMask(),
            mask_map=mask_map,
            no_good_pixels_mask=stats_ctrl.getNoGoodPixelsMask(),
            calc_error_from_input_variance=stats_ctrl.getCalcErrorFromInputVariance(),
            compute_n_image=True)

        for exposure, weight in zip(exposures, weights):
            stacker.add_image(exposure.image, weight=weight)

        stacker.fill_stacked_image(coadd_exposure.image)

        online_image = coadd_exposure.image

        # The unmasked coadd good pixels should match at the <1e-5 level
        # The masked pixels will not be correct for straight image stacking.
        good_pixels = np.where(afw_masked_image.mask.array == 0)

        testing.assert_array_almost_equal(online_image.array[good_pixels],
                                          afw_masked_image.image.array[good_pixels],
                                          decimal=5)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
