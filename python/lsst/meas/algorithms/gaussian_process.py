import george
from george import kernels
import numpy as np

from CloughTocher2DInterpolatorUtils import ctUtils

class GaussianProcessHODLRSolver():

    def __init__(self, variance=1., correlation_length=1., white_noise=0., mean=0.):
        
        self.variance = variance
        self.correlation_length = correlation_length
        self.white_noise  = white_noise
        self.mean = mean

    def fit(self, x_good):
        kernel = self.variance * kernels.ExpSquaredKernel(self.correlation_length, ndim=2)
        self.gp = george.GP(kernel, mean=self.mean,
                            fit_kernel=False,
                            solver=george.HODLRSolver, seed=42)
        self.gp.compute(x_good, yerr=self.white_noise)

    def predict(self, x_bad, y_good):
        y_pred = self.gp.predict(y_good, x_bad, return_var=False, return_cov=False)
        return y_pred


class InterpolateOverDefectGaussianProcess():

    def __init__(self, maskedImage, defects=["SAT"], fwhm=5, block_size=100):

        self.maskedImage = maskedImage
        self.defects = defects
        self.correlation_length = fwhm
        self.block_size = block_size

    def interpolate_over_defects(self):

        nx = self.maskedImage.getDimensions()[0]
        ny = self.maskedImage.getDimensions()[1]
        for x in range(0, nx, self.block_size):
            for y in range(0, ny, self.block_size):
                sub_nx = min(self.block_size, nx - x)
                sub_ny = min(self.block_size, ny - y)
                sub_masked_image = self.maskedImage[x:x+sub_nx, y:y+sub_ny]
                sub_masked_image = self.interpolate_sub_masked_image(sub_masked_image)
                self.maskedImage[x:x+sub_nx, y:y+sub_ny] = sub_masked_image


    def interpolate_sub_masked_image(self, sub_masked_image):

        cut = self.correlation_length * 5
        bad_pixel, good_pixel = ctUtils.findGoodPixelsAroundBadPixels(sub_masked_image, self.defects, buffer=cut)
        # Do nothing if bad pixel is None.
        if np.shape(bad_pixel)[0] == 0:
            return sub_masked_image
        # Do GP interpolation if bad pixel found.
        else:
            # gp interpolation
            mean = np.mean(good_pixel[:,2:])
            white_noise = np.sqrt(np.mean(sub_masked_image.getVariance().array))
            kernel_amplitude = np.var(good_pixel[:,2:]) - white_noise*2


            gp_hodlr = GaussianProcessHODLRSolver(variance=kernel_amplitude, correlation_length=self.correlation_length,
                                                  white_noise=white_noise, mean=mean)
            gp_hodlr.fit(good_pixel[:,:2])
            gp_predict = gp_hodlr.predict(bad_pixel[:,:2], np.squeeze(good_pixel[:,2:]))
            bad_pixel[:,2:] = gp_predict.reshape(np.shape(bad_pixel[:,2:]))

            # update_value
            ctUtils.updateImageFromArray(sub_masked_image.image, bad_pixel)

            return sub_masked_image