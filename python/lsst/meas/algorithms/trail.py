from __future__ import print_function, division, absolute_import

import copy
import numpy as np
import scipy.ndimage.filters

from lsst.pipe.base import Task, Struct
from lsst.pex.config import Config, Field, ListField, ConfigurableField, ChoiceField, ConfigDictField

import lsst.afw.math as afwMath
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.display as afwDisplay

from .hough import HoughTask, Trail, ConstantProfile, DoubleGaussianProfile


__all__ = ["TrailConfig", "TrailTask"]


class MomentLimitConfig(Config):
    limit = Field(dtype=float, doc="Value to use as a limit")
    style = ChoiceField(dtype=str, default="center", doc="Style of limit to apply",
                        allowed={"lower": "apply lower limit",
                                 "center": "apply upper and lower limit",
                                 "upper": "apply upper limit"})

    def test(self, value):
        norm = value/self.limit
        if self.style == "lower":
            return norm > 1.0
        if self.style == "center":
            return np.abs(norm) < 1.0
        if self.style == "upper":
            return norm < 1.0


class TrailConfig(Config):
    hough = ConfigurableField(target=HoughTask, doc="Identification of trails via Hough Transform")

    mask = Field(dtype=str, default="TRAIL", doc="Mask plane to set for trails")
    badMask = ListField(dtype=str, default=["BAD", "CR", "SAT", "INTRP", "EDGE", "SUSPECT"],
                        doc="Mask planes indicating pixels that shouldn't be searched")
    bins = Field(dtype=int, default=4, doc="Binning factor to apply to image before detection")
    scaleDetected = Field(dtype=float, default=10.0, optional=True,
                          doc="Scale detected pixels by this amount.")
    detectionThreshold = Field(dtype=float, default=2.0,
                               doc=("Threshold (in stdev) for pixels to be considered as detected "
                                    "in the binned image"))
    doBackground = Field(dtype=bool, default=False, doc="Median-ring-filter background subtract")
    sigmaSmooth = Field(dtype=float, default=1.0, doc="Gaussian smooth sigma (binned pixels)")
    kernelWidth = Field(dtype=int, default=11, doc="Width of x-cor kernel in pixels")
    growKernel = Field(dtype=float, default=1.4,
                       doc="Repeat with a kernel larger by this fraction (no repeat if 1.0)")
    maskFluxFraction = Field(dtype=float, default=0.5, doc="Fraction of RMS to mask out to")
    maskWidthFudge = Field(dtype=float, default=3.0, doc="Fudge factor for trail mask width. "
                           "This helps account for intermittent trails and non-Gaussian profile.")
    thetaTolerance = Field(dtype=float, default=0.15,
                           doc="Max theta difference for thetaAlignment() routine.")
    # Selection of pixels
    selection = ConfigDictField(keytype=str, itemtype=MomentLimitConfig, default={},
                                doc=("Selections to apply to x-cor values. Keys may be any of: "
                                     "center, theta, ellip, center_perp, skew, skew_perp, b"))
    luminosityLimit = Field(dtype=float, default=2, doc="Lowest luminosity in Std.Dev.")

    kernelSigma = Field(dtype=float, default=7.0, doc="Gauss sigma to use for x-cor kernel.")
    widths = ListField(dtype=float, default=(1.0, 8.0), doc="*unbinned* width of trail to search for.")
    maxPixels = Field(dtype=int, default=10000, doc="Maximum number of pixels to select")
    apertureFactor = Field(dtype=float, default=3.0, doc="Factor of width to use for measurement aperture")

    def setDefaults(self):
        Config.setDefaults(self)
        self.scaleDetected = None
        for name, limit, style in (
                                   ("center", 2.4, "center"),
                                   ("center_perp", 1.2, "center"),
#                                   ("skew", 20.0, "center"),
#                                   ("skew_perp", 10.0, "center"),
#                                   ("ellip", 0.08, "center"),
                                   ("b", 1.4, "center"),
                                   ):
            self.selection[name] = MomentLimitConfig()
            self.selection[name].limit = limit
            self.selection[name].style = style


class TrailTask(Task):
    def __init__(self, *args, **kwargs):
        Task.__init__(self, *args, **kwargs)
        self.makeSubtask("hough")

    def run(self, exposure):
        """Detect satellite trails in exposure

        @param exposure      The exposure to detect in

        @return trails       A SatelliteTrailList object containing detected trails.
        """
        if False:
            binned = self.binImage(exposure)
            smoothedResults = self.smoothImage(binned)
        else:
            smoothedResults = self.smoothImage(exposure)
        selectionResults = self.selectPixels(smoothedResults,
                                             exposure.getPsf().computeShape().getDeterminantRadius())
        trails = self.findTrails(exposure, selectionResults, )
        self.maskTrails(exposure, trails, selectionResults.width, smoothedResults.rms)
        return trails

    def binImage(self, exposure):
        if self.config.bins == 1:
            return exposure.clone()
        binned = exposure.Factory(afwMath.binImage(exposure.getMaskedImage(), self.config.bins))
        binned.setMetadata(exposure.getMetadata())
        binned.setPsf(exposure.getPsf())
        return binned

    def smoothImage(self, exposure):
        image = exposure.getMaskedImage().getImage()
        array = image.getArray()
        mask = exposure.getMaskedImage().getMask()
        detected = mask.getPlaneBitMask("DETECTED")
        maskVal = mask.getPlaneBitMask(self.config.badMask)

        afwDisplay.getDisplay(frame=1).mtv(exposure, title="Original image")

        # Remove masked and detected pixels so we can get a good measure of the noise
        clipped = exposure.clone()
        clippedMask = clipped.getMaskedImage().getMask().getArray()
        clippedImage = clipped.getMaskedImage().getImage().getArray()
        ignore = clippedMask & (maskVal | detected) > 0  # Ignore these pixels
        clippedImage[ignore] = np.nan

        # Amplify detected pixels (because we care about them a lot)
        if self.config.scaleDetected is not None:
            pixelsToScale = (mask.getArray() & detected > 0)
            pixelsToScale |= array > self.config.detectionThreshold*clippedImage.std()
            array[pixelsToScale] *= self.config.scaleDetected

        # Background subtraction
        if self.config.doBackground:
            bg = medianRing(clippedImage, self.config.kernelWidth, 2.0*self.config.sigmaSmooth)
            array -= bg
            clippedImage -= bg

        smoothed = smooth(image, self.config.sigmaSmooth)
        clippedSmoothed = smooth(clipped.getMaskedImage().getImage(), self.config.sigmaSmooth)

        binned = afwMath.binImage(smoothed, self.config.bins)
        clippedBinned = afwMath.binImage(clippedSmoothed, self.config.bins)

        good = np.isfinite(clippedBinned.getArray())
        quartiles = np.percentile(clippedBinned.getArray()[good], (25.0, 75.0))
        rms = 0.74*(quartiles[1] - quartiles[0])

        displayArray(binned.getArray(), frame=2, title="Smoothed binned image")
        displayArray(clippedBinned.getArray(), frame=3, title="Clipped smoothed binned image")

        return Struct(image=binned.getArray(), rms=rms)

        """

int[2.pi.r.A.exp(-0.5*r^2/sigma^2).dr](0,inf)
z = -0.5*r^2/sigma^2 ==> dz = -r.dr/sigma^2 ==> dr = -sigma^2/r.dz
==> int[-sigma^2.2.pi.A.exp(z).dz](0,-inf)
= [-sigma^2.2.pi.A.exp(z)](0,-inf) = sigma^2.2.pi.A

v/sigma^2/2.pi




        """


    def _makeCalibrationImage(self, psfSigma, width, kernelFactor):
        """Make a fake satellite trail with the PSF to calibrate moments we measure.

        @param psfSigma       Gaussian sigma for a double-Gaussian PSF model (in pixels)
        @param width          Width of the trail in pixels (0.0 for PSF alone, but wider for aircraft trails)

        @return calImg        An afwImage containing a fake satellite/aircraft trail
        """

        kernelWidth = 2*int((kernelFactor*self.config.kernelWidth)//2) + 1

        # tricky.  We have to make a fake trail so it's just like one in the real image

        # To get the binning right, we start with an image 'bins'-times too big
        size = self.config.bins*kernelWidth
        calImage = afwImage.ImageF(size, size)
        calArray = calImage.getArray()

        # Make a trail with the requested (unbinned) width
        calTrail = Trail((self.config.bins*kernelWidth)//2 - 0.5, 0.0*afwGeom.degrees)

        # for wide trails, just add a constant with the stated width
        if False and width > 8.0*psfSigma:
            profile = ConstantProfile(1.0, width)
            insertWidth = width
        # otherwise, use a double gaussian
        else:
            profile = DoubleGaussianProfile(1.0, width/2.0 + psfSigma)
            insertWidth = self.config.bins*kernelWidth # 4.0*(width/2.0 + psfSigma)
        calTrail.insert(calArray, profile, min(insertWidth, size))

        # Now bin and smooth, just as we did the real image
        smoothed = smooth(calImage, self.config.sigmaSmooth)
        calArray = afwMath.binImage(smoothed, self.config.bins).getArray()

        return calArray

    def selectPixelsForWidth(self, width, kernelFactor, psfSigma, moments, rms):
        calImage = self._makeCalibrationImage(psfSigma, width, kernelFactor)
        calMoments = self.calculateMoments(calImage, kernelFactor, isCalibration=True)

        selected = (moments.sumI - calMoments.sumI)/(self.config.luminosityLimit*rms) > 1.0
        displayArray(selected, frame=3, title="Brightness selection")
#        import pdb;pdb.set_trace()
        for name, limit in self.config.selection.items():
            pixels = limit.test(getattr(moments, name) - getattr(calMoments, name))
            displayArray(pixels, frame=3, title="Selection for %s" % (name,))
#            import pdb;pdb.set_trace()
            selected &= limit.test(getattr(moments, name) - getattr(calMoments, name))

        num = selected.sum()
        if num < self.config.maxPixels:
            return selected
        # Take the brightest pixels because they have the most potential to screw things up
        threshold = np.percentile(moments.image[selected], 100.0*(1.0 - self.config.maxPixels/num))
        selected &= moments.image > threshold
        return selected

    def calculateMoments(self, image, kernelFactor, isCalibration=False):
        kernelWidth = 2*int((kernelFactor*self.config.kernelWidth)//2) + 1
        kernelSigma = kernelFactor*self.config.kernelSigma
        return MomentManager(image, kernelWidth=kernelWidth, kernelSigma=kernelSigma,
                             isCalibration=isCalibration)

    def selectPixels(self, smoothed, psfSigma):
        image = smoothed.image
        rms = smoothed.rms

        selection = np.ones(image.shape, dtype=bool)
        for kernelFactor in set((1.0, self.config.growKernel)):
            self.log.debug("Getting moments with kernel grown by factor of %.1f" % (kernelFactor,))
            moments = self.calculateMoments(image, kernelFactor)

            if False:
                displayArray(moments.sumI, frame=7, title="sumI")
                displayArray(moments.imageMoment.ix, frame=8, title="Ix")
                displayArray(moments.imageMoment.iy, frame=9, title="Iy")
                ellip, theta, b = momentToEllipse(moments.imageMoment.ixx, moments.imageMoment.iyy,
                                                  moments.imageMoment.ixy)
                displayArray(ellip, frame=10, title="Ellip")
                displayArray(b, frame=11, title="b")
                displayArray(np.sqrt(moments.imageMoment.ixxx**2 + moments.imageMoment.iyyy**2), frame=12, title="Skew")
                print("KernelFactor=%f bins=%d" % (kernelFactor, self.config.bins))
                import pdb
                pdb.set_trace()

            maxNumPixels = 0
            bestWidth = self.config.widths[0]
            kernelSelection = np.zeros(image.shape, dtype=bool)
            for width in sorted(self.config.widths):
                pixels = self.selectPixelsForWidth(width, kernelFactor, psfSigma, moments, rms)
                numPixels = pixels.sum()
                self.log.debug("For kernelFactor=%f, width=%f, got %d candidate pixels",
                              kernelFactor, width, numPixels)
                if numPixels > maxNumPixels:
                    maxNumPixels = numPixels
                    bestWidth = width
                kernelSelection |= pixels

            if False:
                displayArray(pixels, frame=4, title="Selected pixels for kernel+width")
                import pdb
                pdb.set_trace()

            selection &= kernelSelection

        displayArray(selection, frame=4, title="Selected pixels for kernel")

        # XXX "moments" is specific to the last kernelFactor, and not generic;
        # Should this be done for each kernelFactor?
        if False:
            selection = self.selectThetaAlignment(selection, moments)

        return Struct(pixels=selection, moments=moments, width=bestWidth)

    def selectThetaAlignment(self, pixels, moments):
        xx, yy = np.meshgrid(np.arange(pixels.shape[1], dtype=int), np.arange(pixels.shape[0], dtype=int))
        maxSeparation = min([x/2 for x in pixels.shape])
        thetaMatch, newTheta = thetaAlignment(moments.theta[pixels], xx[pixels], yy[pixels],
                                              tolerance=self.config.thetaTolerance,
                                              limit=3, maxSeparation=maxSeparation)
        moments.theta[pixels] = newTheta
        pixels[pixels] = thetaMatch

        displayArray(pixels, frame=5, title="Selected pixels after theta alignment")

        return pixels

    def findTrails(self, exposure, selection):
        pixels = selection.pixels
        xx, yy = np.meshgrid(np.arange(pixels.shape[1], dtype=int), np.arange(pixels.shape[0], dtype=int))
        theta = selection.moments.theta
        rMax = np.linalg.norm(pixels.shape)
        return self.hough.run(theta[pixels], xx[pixels], yy[pixels], self.config.bins, rMax)

    def maskTrails(self, exposure, trails, width, rms):
        mask = exposure.getMaskedImage().getMask()
        mask.addMaskPlane(self.config.mask)
        maskVal = mask.getPlaneBitMask(self.config.mask)
        psfSigma = exposure.getPsf().computeShape().getDeterminantRadius()
        # rms is from the binned, smoothed image, so to get the noise on the unbinned unsmoothed image, the
        # rms should be increased by a factor of sqrt(bins^2) to account for the binning, and a factor of
        # 2*sqrt(pi)*sigmaSmooth to account for the smoothing.
        imageRms = rms*self.config.bins*self.config.sigmaSmooth*np.sqrt(2.0*np.pi)
        for trail in trails:
            measureWidth = max(width, psfSigma)
            measurements = trail.measure(exposure, 1, self.config.apertureFactor*measureWidth)

            # Determine width that will be sufficient to mask a Gaussian profile down to target level:
            # target = totalFlux/length/sigma/sqrt(2pi)*exp(-0.5*x^2/sigma^2)
            # ==> x = sigma*sqrt(-2*ln(sqrt(2pi)*target*length*sigma/totalFlux))
            sigma = measurements.width
            avgFlux = measurements.flux/trail.length(exposure.getWidth(), exposure.getHeight())
            target = self.config.maskFluxFraction*imageRms
            value = np.sqrt(2*np.pi)*target*sigma/avgFlux
            if value < 1:
                useWidth = sigma*np.sqrt(-2*np.log(value))
                useWidth *= self.config.maskWidthFudge
            else:
                useWidth = width
            useWidth = max(useWidth, psfSigma)
            self.log.info("Masking trail %s (flux=%f width=%f) with width %s", trail, measurements.flux,
                          measurements.width, useWidth)
            trail.setMask(exposure, 2.0*useWidth, maskVal)

        frame = 6
        afwDisplay.getDisplay(frame).mtv(exposure)
        trails.display(exposure.getDimensions(), frame=frame)


def smooth(image, sigma):
    """Smooth an image with a circular Gaussian kernel

    Parameters
    ----------
    image : `lsst.afw.image.Image`
        The image to smooth.
    sigma : float
        Gaussian 'sigma' of the circular smoothing kernel.

    Returns
    -------
    out : `lsst.afw.image`
        Smoothed image.
    """
    out = image.Factory(image.getBBox())
    k = 2*int(6.0*sigma) + 1
    kk = np.arange(k) - k//2
    gauss = (1.0/np.sqrt(2.0*np.pi))*np.exp(-kk*kk/(2.0*sigma))

    # Smooth with separable kernel
    temp = scipy.ndimage.filters.correlate1d(image.getArray(), gauss, mode='reflect')
    scipy.ndimage.filters.correlate1d(temp, gauss, mode='reflect', axis=0, output=out.getArray())

    return out


def medianRing(image, radius, width):
    """Apply a median annulus filter

    Parameters
    ----------
    image : ndarray
        Image pixels as a numpy array.
    radius : float
        Radius of the annulus.
    width : float
        Width of the annulus.

    Returns
    -------
    filtered_image : ndarray
        Image, filtered by median annulus
    """
    k = 2*int(radius + width) + 1
    ring = np.zeros((k, k), dtype=bool)
    a = np.arange(k) - k//2
    x, y = np.meshgrid(a, a)
    r2 = x*x + y*y
    w = (r2 > radius**2) & (r2 < (radius + width)**2)
    ring[w] = True
    return scipy.ndimage.filters.median_filter(image, footprint=ring)


def momentConvolve2d(data, k, sigma, middleOnly=False):
    """Convolve an image with coefficient kernels to compute local 'moments'

    @param data       The input image
    @param k          A vector of indices (e.g. -3,-2,-1,0,1,2,3 )
    @param sigma      Gaussian sigma for an overall smoothing to avoid blowing up higher-order moments
    @param middleOnly Boolean to return the central pixel only (used for calibration images)

    return ImageMoment  A container with attributes for each moment: .i0 .ix .iy .ixx .iyy .ixy etc.

    Each of these convolutions uses a separable kernel, and many share a common convolution
    in at least one dimension.
    """

    # moments are  e.g.   sum(I*x) / sum(I)

    gauss = np.exp(-k**2/(2.0*sigma**2))

    kk = k*k
    k3 = kk*k

    mode = 'reflect'
    correlate1d = scipy.ndimage.filters.correlate1d

    # start with convolutions with our Gaussian in separately in X and Y
    gaussX = correlate1d(data, gauss, mode=mode)
    gaussY = correlate1d(data, gauss, mode=mode, axis=0)

    # zeroth order moment (i.e. a sum), convolve the X gaussian along Y
    sumI = correlate1d(gaussX, gauss, mode=mode, axis=0)
    sumI[np.where(sumI == 0)] = np.finfo(float).tiny

    # normalize up front
    gaussX /= sumI
    gaussY /= sumI

    # Now use gaussX and gaussY to get the moments
    ix = correlate1d(gaussY, gauss*k, mode=mode)
    iy = correlate1d(gaussX, gauss*k, mode=mode, axis=0)
    ixx = correlate1d(gaussY, gauss*kk, mode=mode)
    iyy = correlate1d(gaussX, gauss*kk, mode=mode, axis=0)

    # cross term requires special attention.  Start from scratch.
    ixy0 = correlate1d(data, gauss*k, mode=mode)
    ixy = correlate1d(ixy0, gauss*k, mode=mode, axis=0) /sumI

    # don't bother with 3rd order cross terms
    ix3 = correlate1d(gaussY, gauss*k3, mode=mode)
    iy3 = correlate1d(gaussX, gauss*k3, mode=mode, axis=0)

    values = dict(i0=sumI, ix=ix, iy=iy, ixx=ixx, iyy=iyy, ixy=ixy, ixxx=ix3, iyyy=iy3)
    if middleOnly:
        ny, nx = data.shape
        values = {k: v[ny//2, nx//2] for k, v in values.items()}
    return Struct(**values)


def momentToEllipse(ixx, iyy, ixy, loClip=0.1):
    """Convert moments to ellipse parameters (numpy-safe)

    @param ixx     2nd moment x
    @param iyy     2nd moment y
    @param ixy     2nd moment xy
    @param loClip  Minium value to accept for either A (semi-major) or B (semi-minor)

    @return ellip, theta, B    Ellipticity 1-B/A, Pos.Angle theta, and semi-minor axis B
    """

    tmp = 0.5*(ixx + iyy)
    diff = ixx - iyy
    tmp2 = np.sqrt(0.25*diff**2 + ixy**2)
    a2 = np.clip(tmp + tmp2, loClip, None)
    b2 = np.clip(tmp - tmp2, loClip, None)
    ellip = 1.0 - np.sqrt(b2/a2)
    theta = 0.5*np.arctan2(2.0*ixy, diff)

    return ellip, theta, np.sqrt(b2)


class MomentManager(object):
    """Handle calculation of moments for all pixels in an image.

    We'll try to do this in an on-demand way, so we only calculate something
    if it's being used.
    """

    keys = "sumI", "center", "theta", "ellip", "center_perp", "skew", "skew_perp", "b"

    def __init__(self, image, kernelWidth, kernelSigma, isCalibration=False):
        """Construct

        @param image          The image with moments we want computed.
        @param kernelWidth    The kernel width to use in pixels
        @param kernelSigma    Gaussian sigma for a weight function applied multiplicatively to the kernel
        @param isCalibration  Is this a calibration image?
                                (If so, don't convolve, just get the calib pixel [the center])

        """
        self.image = image
        self.shape = image.shape
        self.isCal = isCalibration
        self.std = image.std()
        self.kernelWidth = kernelWidth
        self.kernelSigma = kernelSigma

        # properties
        self._imageMoment = None
        self._sumI = None
        self._center = None
        self._center_perp = None
        self._ellip = None
        self._theta = None
        self._skew = None
        self._skew_perp = None
        self._b = None

    def _toEllipse(self):
        if (self._ellip is None):
            ixx, iyy, ixy = self.imageMoment.ixx, self.imageMoment.iyy, self.imageMoment.ixy
            self._ellip, self._theta, self._b = momentToEllipse(ixx, iyy, ixy)

    @property
    def imageMoment(self):
        """Compute the convolutions"""
        if self._imageMoment is None:
            kx = np.arange(self.kernelWidth) - self.kernelWidth//2
            self._imageMoment = momentConvolve2d(self.image, kx, self.kernelSigma, middleOnly=self.isCal)
        return self._imageMoment

    @property
    def sumI(self):
        """Get the sum of pixel values"""
        if self._sumI is None:
            self._sumI = 0.0 if self.isCal else self.image
        return self._sumI

    @property
    def center(self):
        """Get the centroid offset (1st moment)"""
        if self._center is None:
            ix, iy = self.imageMoment.ix, self.imageMoment.iy
            self._center = np.sqrt(ix**2 + iy**2)
        return self._center

    @property
    def theta(self):
        """Get the position angle w.r.t. the x-axis in radians."""
        if self._theta is None:
            self._toEllipse()
        return self._theta

    @property
    def center_perp(self):
        """Get the centroid offset (1st moment) perpendicular to the alignment."""
        if self._center_perp is None:
            ix, iy = self.imageMoment.ix, self.imageMoment.iy
            self._center_perp = np.abs(ix*np.sin(self.theta) - iy*np.cos(self.theta))
        return self._center_perp

    @property
    def ellip(self):
        """Get the ellipticity: e = 1 - B/A """
        if self._ellip is None:
            self._toEllipse()
        return self._ellip

    @property
    def skew(self):
        """Get the skewness (3rd moment)"""
        if self._skew is None:
            ixxx, iyyy = self.imageMoment.ixxx, self.imageMoment.iyyy
            self._skew = np.sqrt(ixxx**2 + iyyy**2)
        return self._skew

    @property
    def skew_perp(self):
        """Get the skewness (3rd moment) perpendicular to the alignment."""
        if self._skew_perp is None:
            ixxx, iyyy = self.imageMoment.ixxx, self.imageMoment.iyyy
            self._skew_perp = np.abs(ixxx*np.sin(self.theta) - iyyy*np.cos(self.theta))
        return self._skew_perp

    @property
    def b(self):
        """Get the 'semi-minor axis', B"""
        if self._b is None:
            self._toEllipse()
        return self._b


class MomentLimit(Struct):
    def __init__(self, name, value, limitType):
        """Construct.

        @param name       Name of the limit (must be in MomentManager.keys)
        @param value      The value to use as a limit.
        @param limitType  Is the limit 'lower', 'center', or 'upper' limit?
        """
        Struct.__init__(self, name=name, value=value, limitType=limitType)


class PixelSelector(list):
    """A simple pixel selector.

    Inherit from a list, and we'll contain a list of our MomentLimit objects.
    We'll go through our list, and any pixels which are numerically within the limits
    for all MomentLimit objects "pass" and are kept.

    In the end, we return a boolean image with True set for accepted pixels.
    """

    def __init__(self, momentManager, calMomentManager):
        """Construct

        @param momentManager    MomentManager for the image we're selecting from
        @param calMomentManager The MomentManager for the calibration image.
        """
        list.__init__(self)

        assert(momentManager.kernelWidth == calMomentManager.kernelWidth)
        self.keys = copy.copy(MomentManager.keys)

        self.momentManager = momentManager
        self.calMomentManager = calMomentManager

    def append(self, limit):
        """Overload our parent list's append so that we can verify the MomemntList being appended

        If all is well, we'll call our parent's append()

        @param limit   The MomentLimit object being added to this selector.
        """
        if limit.name not in self.keys:
            raise ValueError("Limit name must be in:" + str(self.keys))
        limitTypes = ('lower', 'center', 'upper')
        if limit.limitType not in limitTypes:
            raise ValueError("Limit limitType must be in:" + str(limitTypes))
        super(PixelSelector, self).append(limit)


    def _test(self, limit):
        """Helper method to determine pass/fail for a specified MomentLimit

        @param limit   The MomentLimit to test.

        @return test   The result image of the test (passing pixels are True)
        """

        val         = getattr(self.momentManager,    limit.name)
        expectation = getattr(self.calMomentManager, limit.name)
        norm        = (val - expectation)/np.abs(limit.value)
        if limit.limitType == 'lower':
            test = norm > 1.0
        elif limit.limitType == 'center':
            test = np.abs(norm) < 1.0
        elif limit.limitType == 'upper':
            test = norm < 1.0
        return test

    def getPixels(self, maxPixels):
        """Check against all MomentLimit and return an image with pixels which passed all tests.

        @param maxPixels   Limit the number of pixels
                           (not implemented here as there's no obvious way to sort them)
        """
        selected = np.ones(self.momentManager.shape, dtype=bool)
        for limit in self:
            test = self._test(limit)

            if False:
                displayArray(test, frame=13, title="%s selection" % (limit.name,))
                import pdb
                pdb.set_trace()

            selected &= test

        num = selected.sum()
        if num < maxPixels:
            return selected

        if False:
            # Take the brightest pixels because they have the most potential to screw things up
            threshold = np.percentile(self.momentManager.image[selected], 100.0*(1.0 - maxPixels/num))
            selected &= self.momentManager.image > threshold
        else:
            where, = np.where(selected)
            seed = 12345
            rng = np.random.RandomState(seed)
            where = rng.shuffle(where)[:maxPixels]
            selected = selected[where]
        return selected


class ProbabilityPixelSelector(PixelSelector):
    """A P-Value based pixel selector.

    This serves the same purpose as the PixelSelector, but computes a p-value for each pixel.
    The MomentLimits are used as 1-sigma thresholds, and the resulting sum of log(p) values
    is computed.  Pixels meeting a specified threshold are kept.
    """

    def __init__(self, *args, **kwargs):
        """Construct.
        """
        self.thresh = kwargs.pop('thresh')
        PixelSelector.__init__(self, *args, **kwargs)
        super(PValuePixelSelector, self).__init__(*args, **kwargs)

        # cache some value in case the user wants to test the same thing twice
        # e.g. an upper limit and a lower limit.
        self.cache = {}
        self.done = []

    def _test(self, limit):
        """Test the significance of pixels for this MomentLimit.

        @param limit  The MomentLimit to test.
        """

        if limit.name in self.cache:
            delta, delta2, neg = self.cache[limit.name]
            expectation = getattr(self.calMomentManager, limit.name)
        else:
            val = getattr(self.momentManager,    limit.name)
            expectation = getattr(self.calMomentManager, limit.name)
            delta = val - expectation
            neg = delta <= 0.0
            delta2 = delta**2
            self.cache[limit.name] = (delta, delta2, neg)

        # If z is normalized, our Gaussian is P = exp(-z**2/2)
        # Or ... z**2 = -2*log(P)

        divByZeroValue = 0.001

        zz = delta2/(limit.value**2) + divByZeroValue

        # Some of these aren't real probabilities, just functions that go to 1 or 0 as needed.
        # This would be trivial with exp() and log() functions, but they're very expensive,
        # so these approximations use only simple arithmatic.

        # Go to 1 for z > 1.  the function is very close to 1 - exp(-x**2)
        # it would go to 1 at large z for both +ve and -ve, so we have to suppress the negative side.
        if limit.limitType == 'lower':
            neg2logp = 1.0/zz
            neg2logp[neg] = 1.0/divByZeroValue

        elif limit.limitType == 'center':
            neg2logp = zz

        # This is the opposite of 'lower'.
        # it keeps the values below z~1 and suppresses those above z=2
        elif limit.limitType == 'upper':
            neg2logp = zz
            neg2logp[neg] = divByZeroValue

        else:
            raise ValueError("Unknown limit type.")

        return neg2logp

    def getPixels(self, maxPixels):
        """Get the pixels which pass all MomentLimit tests.

        @param maxPixels   Return no more than this many 'pass' pixels.  (Sorting by p-value)
        """

        n = 0
        neg2logp = np.zeros(self.momentManager.shape)
        for limit in self:
            n += 1
            val = self._test(limit)
            neg2logp += val
        logp = -0.5*neg2logp

        # logp = -z**2/2 is contributed by each parameter considered
        # So to keep n-sigma points use:
        # 1-sigma:    0.5   # 67% of real candidates accepted
        # 0.7-sigma:  0.25  # 52% accepted
        # 0.5-sigma:  0.125 # 38% accepted
        thresh1 = self.thresh or -0.125*n #-0.5*n

        thresh1 = -0.5*n

        ret = logp > thresh1

        nth = 1.0*maxPixels/logp.size
        if ret.sum() > maxPixels:
            thresh2 = np.percentile(logp, 100*(1.0-nth))
            ret = logp > thresh2
        return ret


def thetaAlignment(theta, x, y, limit=3, tolerance=0.15, maxSeparation=None):
    """A helper function to cull the candidate points.

    @param theta      ndarray of thetas
    @param x          ndarray of x pixel coordinates
    @param y          ndarray of y pixel coordinates

    The basic idea here is that for any pair of points, each has a local measure of theta
    and also vector connecting the two points, which also defines a third theta.
    All of these should agree, so we can eliminate any candidate points for which
    either local theta is more than 'tolerance' different from that defined by the
    dx,dy pixel coordinates.

    This only gets you so far.  With ~1000 candidate points, each one will have ~10 chance
    alignments for a reasonable tolerance.  Though the local theta values are uncertain at the +/-0.1
    level, the pixel coordinate-based thetas are much more precise.  So from those ~10, we can
    search for pairs which have pixelTheta values which line up 'too well'.

    The final step is choosing an delta-angle to define as the smallest separation we'd expect
    to see between any two points.  If we assume the points are uniformly distributed, the probability
    1 point will we be found within 'delta' of another is e^(-2*phi) where fill-factor phi=n*delta/range.
    This might look familiar as the Poisson prob for zero events with rate u=-2*phi:

       P(x=k) = u^k/k! e^(-u)   -->  P(x=0) = e^(-u)

    What we're saying is that if we sit on one of n points in a region 'range',
    the probability of observing 0 points within 'delta' of our position is:

       P(x=0) = e^(-2*n*delta/range)

    The factor of 2 arises because a point may preceed or follow.
    """

    num = len(theta)

    dx = np.subtract.outer(x, x)
    dy = np.subtract.outer(y, y)
    dydx = dy/(dx + 0.01)
    thetaXY = np.arctan(dydx)

    aligned1 = np.abs(thetaXY - theta) < tolerance
    aligned1 |= np.abs(thetaXY - (theta + 1.0*np.pi)) < tolerance
    aligned1 |= np.abs(thetaXY - (theta - 1.0*np.pi)) < tolerance
    aligned2 = np.abs(thetaXY.transpose() - theta).transpose() < tolerance
    aligned2 |= np.abs(thetaXY.transpose() - (theta + 1.0*np.pi)).transpose() < tolerance
    aligned2 |= np.abs(thetaXY.transpose() - (theta - 1.0*np.pi)).transpose() < tolerance
    if maxSeparation:
        dist = dx**2 + dy**2
        closish = dist < maxSeparation**2
        aligned = aligned1 & aligned2 & closish
    else:
        aligned = aligned1 & aligned2

    nNearNeighbours = np.zeros(num)
    newThetas = copy.copy(theta)

    # Using variable names  pCloseNeighbour = e^(2*nCand*closeNeighbourTolerance/tolerance)
    pZeroCloseNeighbour = 0.67
    # compute the closeNeighbourTolerance for which we expect 1 collision
    phi = -0.5*np.log(pZeroCloseNeighbour)

    nCand = aligned.sum(axis=1)
    nCand[nCand == 0] = 1
    closeNeighbourTolerance = phi*tolerance/nCand

    cut = max(limit, 2)
    w, = np.where(nCand > cut)
    for i in w:
        # this will fail near +-pi
        pixelTheta = np.sort(thetaXY[i,aligned[i,:]])
        diff = np.diff(pixelTheta)
        select = (diff < closeNeighbourTolerance[i]) & (diff > 1.0e-8)
        nNearNeighbours[i] = select.sum()
        if nNearNeighbours[i] >= limit:
            pixTheta = 0.5*(pixelTheta[:-1] + pixelTheta[1:])[select]
            idxMedian = len(pixTheta)//2
            newThetas[i] = pixTheta[idxMedian]

    isCandidate = nNearNeighbours >= limit

    return isCandidate, newThetas


def displayArray(array, frame, title=None):
    if False:
        return
    typemap = {np.dtype(t1): t2 for t1, t2 in (("int32", afwImage.ImageI),
                                               ("bool", afwImage.ImageU),
                                               ("float32", afwImage.ImageF),
                                               ("float64", afwImage.ImageD),)}

    dummy = typemap[array.dtype](*reversed(array.shape))
    dummy.getArray()[:] = array
    afwDisplay.getDisplay(frame).mtv(dummy, title=title)
