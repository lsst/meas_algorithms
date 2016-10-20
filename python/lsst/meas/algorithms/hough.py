from __future__ import print_function, division

import copy
import numpy as np
import scipy.ndimage

from lsst.pipe.base import Task, Struct
from lsst.pex.config import Config, Field

import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom


class HoughConfig(Config):
    bins = Field(dtype=int, default=200, doc="Number of bins to use in r,theta space")
    threshold = Field(dtype=int, default=40, doc="Minimum number of 'hits' in a Hough bin for a detection.")
    maxPoints = Field(dtype=int, default=10000,
                      doc="Maximum number of points to allow (solution gets slow for >> 1000)")
    nIter = Field(dtype=int, default=1,
                  doc="Number of times to iterate the solution (more than ~1-3 rarely makes a difference)")
    maxResid = Field(dtype=float, default=5.5,
                     doc="The max coordinate inter-quartile range (pixels) to accept for a solution")
    overlapRange = Field(dtype=float, default=0.2,
                         doc=("Range either side of theta = 0 and 2*pi to be wrapped around."
                              "This allows good solutions to be found very near theta = 0 and theta=2pi"))
    quadraticLimit = Field(dtype=float, default=2.5e-3, doc="Maximum quadratic coefficient for trails")

class HoughTask(Task):
    ConfigClass = HoughConfig

    def twoPiOverlap(self, thetaIn, arrays=None, overlapRange=0.2):
        """
        Take any thetas near zero and *copy* them to above 2*pi.
        Take any thetas near 2*pi and *copy* them to near 0

        Sometimes the theta we're looking for is near 0 or 2*pi and we'd like
        a continue sample of points.  Otherwise, we can end-up with the
        same theta yielding two solutions.
        """
        w_0, = np.where(thetaIn <= overlapRange)
        w_2pi, = np.where(2.0*np.pi - thetaIn <= overlapRange)

        theta = np.append(thetaIn, thetaIn[w_0] + 2.0*np.pi)
        theta = np.append(theta, thetaIn[w_2pi] - 2.0*np.pi)

        # if we're given other arrays, append in the same way
        outArray = []
        if arrays:
            for arr in arrays:
                tmp = np.append(arr, arr[w_0])
                tmp = np.append(tmp, arr[w_2pi])
                outArray.append(tmp)

        return theta, outArray

    def improveCluster(self, theta, x, y):
        """
        @theta  angles
        @x      x pixel coordinates
        @y      y pixel coordinates

        return list of (r_i,theta_i)

        Due to noise in the original image, the theta values are never quite aligned
        with the trail.  The conversion to the normal form r = x*cos(t) + y*sin(t)
        shows that for a given x,y pixel; errors in theta are sinusoids in r(theta).
        So a cluster in r,theta often has streaks of points passing through the true r,theta
        of the trail.  Since we know x,y and the form of r(theta) very well, we
        can compute dr/dtheta = -x*sin(t) + y*cos(t) for each point in r,theta space.  This is a linear
        approximation to r(theta) near the point.  The idea behind this improvement strategy
        is that points coming from the same (real) trail will share a common point
        in r,theta (the premise of the Hough Transform).  For one point, we only know it lies on
        the r(theta) curve, but for two ... the intersection of the curves is a better estimate of r,theta
        for each contributing point.
        """
        rdist = 100
        tdist = 0.15

        t = theta.copy()
        r = x*np.cos(t) + y*np.sin(t)

        dx = np.subtract.outer(x, x)
        dy = np.subtract.outer(y, y)
        dd = np.sqrt(dx**2 + dy**2)

        w0 = (dx == 0)
        dx[w0] = 1.0

        # this is the theta we get if we just draw a line between points in pixel space
        intercept = y - (dy/dx)*x
        sign = np.sign(intercept)
        pixel_theta = np.arctan(dy/dx) + sign*np.pi/2.0
        pixel_theta[(pixel_theta < 0.0)] += np.pi

        # convert each point to a line in r,theta space
        # drdt is slope and 'b' is the intercept
        drdt = -x*np.sin(t) + y*np.cos(t)
        b = r - t*drdt

        # get the distance between points in r,theta space "good" pairs are close together
        isGood = (np.abs(np.subtract.outer(t, t)) < tdist) & (np.abs(np.subtract.outer(r, r)) < rdist)
        isBad = ~isGood

        # solve for the intersections between all lines in r,theta space
        dm = np.subtract.outer(drdt, drdt).transpose()
        parallel = (np.abs(dm) < 0.01)
        dm[parallel] = 1.0
        tt = np.subtract.outer(b, b)/dm
        rr = b + tt*drdt

        tt[isBad] = 0.0
        rr[isBad] = 0.0

        # Create a weight vector to use in a weighted-average
        w = np.zeros(tt.shape) + 1.0e-7

        # weight by the pixel distance (farther is better, less degenerate)
        w[isGood] += dd[isGood]

        # de-weight points that have moved us far from where they started
        dtr = np.abs((tt - t)*(rr - r)) + 1.0
        # de-weight points where the solved theta disagrees with delta-x,delta-y in pixel space
        theta_discrepancy = np.abs(tt - pixel_theta) + 0.01

        w[isGood] /= (dtr*theta_discrepancy)[isGood]

        # for each point, the 'best' theta is the weighted-mean of all places it's line intersects others.
        t_new = np.average(tt, axis=0, weights=w)
        r_new = x*np.cos(t_new) + y*np.sin(t_new)

        # use original values for things that didn't converge
        t0 = (t_new < 1.0e-6)
        t_new[t0] = t[t0]
        r0 = (r_new < 1.0e-6)
        r_new[r0] = r[r0]

        return r_new, t_new, r, x, y

    def hesseBin(self, r0, theta0, rMax=4096):
        """Bin r,theta values to find clusters above a threshold

        @param r0         List of r values
        @param theta0     List of theta values
        @param rMax       Specify the range in r

        @return bin2d, rEdge, tEdge, rsGood, tsGood, idxGood

        bin2d    the 2D binned data
        rEdge    the edges of the bins in r
        tEdge    the edges of the bins in theta
        rsGood   List of best r values for any loci
        tsGood   List of best theta values for any loci
        idxGood  List of indices (from the input arrays) for any points contributing to any loci

        In principal, this function is simple ... bin-up the values of r,theta and see
        if any of the bins have more than 'thresh' points.  Take the input r,thetas which landed
        in any such bin and use them to get the 'best' value of r,theta ... mean, median, whatever.

        """
        r = r0.copy()
        theta = theta0.copy()

        # eliminate any underdesirable r,theta values
        # namely theta~0.0 and r > rMax
        accept = (np.abs(theta) > 1.0e-2) & (np.abs(r) > 1.0*rMax/self.config.bins)

        # This builds a 2d histogram and labels any bins with count level above a threshold
        overlapRange = self.config.overlapRange
        bin2d, rEdge, tEdge = np.histogram2d(r[accept], theta[accept],
                                             bins=(self.config.bins, self.config.bins),
                                             range=((0.0, rMax), (-overlapRange, overlapRange + 2.0*np.pi)))
        locus, numLocus = scipy.ndimage.label(bin2d > self.config.threshold, structure=np.ones((3, 3)))

        # Now check each locus and get the points which contributed to it, or its immediate neighbours
        rs, ts, idx, drs, dts = [], [], [], [], []
        for i in range(numLocus):
            label = i + 1
            loc_r, loc_t = np.where(locus == label)

            iThetaPeak, iRPeak = 0.0, 0.0
            maxVal = 0.0
            for rr, tt in zip(loc_r, loc_t):
                val = bin2d[rr, tt]
                if val > maxVal:
                    maxVal = val
                    iThetaPeak = tt
                    iRPeak = rr

                # iThetaPeak,iRPeak  is the peak count for this label in bin2d

                # get the indices for a 3x3 box with iThetaPeak,iRPeak at the center
                iThetaMin = max(iThetaPeak - 1, 0)
                iThetaMax = min(iThetaPeak + 1, self.config.bins - 1)
                iRMin = max(iRPeak - 1, 0)
                iRMax = min(iRPeak + 1, self.config.bins - 1)

            tlo, thi = tEdge[iThetaMin], tEdge[iThetaMax + 1]
            rlo, rhi = rEdge[iRMin], rEdge[iRMax + 1]

            # for this locus, use the median r,theta for points within the 3x3 box around the peak
            centeredOnPeak = (theta >= tlo) & (theta < thi) & (r >= rlo) & (r < rhi)
            tTmp = np.median(theta[centeredOnPeak])
            dtTmp = theta[centeredOnPeak].std()
            rTmp = np.median(r[centeredOnPeak])
            drTmp = r[centeredOnPeak].std()

            # don't accept theta < 0 or > 2pi
            if tTmp < 0.0 or tTmp > 2.0*np.pi:
                continue

            rs.append(rTmp)
            drs.append(drTmp)
            ts.append(tTmp)
            dts.append(dtTmp)

            # keep a boolean array ID'ing the points which contributed
            w = (theta0 >= tlo) & (theta0 < thi) & (r0 >= rlo) & (r0 < rhi)
            idx.append(w)

        # check for wrapped-theta doubles,
        # - pick the one with the lowest stdev
        # - this is rare, but a bright near-vertical trail can be detected near theta=0 *and* theta=2pi
        # --> the real trail is rarely exactly vertical, so one solution will not converge nicely.
        #     ... the stdev of thetas will be wider by a factor of "a lot", say ~10x
        n = len(rs)
        kill_list = []
        for i in range(n):
            for j in range(i,n):
                dr = abs(rs[i] - rs[j])
                dt = abs(ts[i] - ts[j])
                # if this pair is close in r, but differs by ~2pi in theta, it's the same thing
                # detected twice.
                if dr < 20 and dt > 1.9*np.pi:
                    # the one with great theta-scatter is 'bad'
                    bad = i if dts[i] > dts[j] else j
                    kill_list.append(bad)

        rsGood, tsGood, idxGood = [],[],[]
        for i in range(n):
            if i in kill_list:
                continue
            rsGood.append(rs[i])
            tsGood.append(ts[i])
            idxGood.append(idx[i])

        return bin2d, rEdge, tEdge, rsGood, tsGood, idxGood

    def run(self, thetaIn, xIn, yIn, binning=1.0, rMax=None, seed=12345):
        """Compute the Hough Transform

        @param  thetaIn        The local theta values at pixel locations (should be within 0.2 of true value)
        @param  xIn            The x pixel coordinates corresponding to thetaIn values
        @param  yIn            the y pixel coordinates corresponding to thetaIn values

        @return solutions      A HoughSolutionList with an entry for each locus found.
        """
        rIn, thetaIn = hesseForm(thetaIn, xIn, yIn)

        # wrap the points so we don't have a discontinuity at 0 or 2pi
        theta0, (r0, x0, y0) = self.twoPiOverlap(thetaIn, (rIn, xIn, yIn))

        nPoints = len(r0)
        if nPoints == 0:
            return TrailList()

        rng = np.random.RandomState(seed)
        r, theta, x, y = r0, theta0, x0, y0
        if nPoints > self.config.maxPoints:
            idx = np.arange(nPoints, dtype=int)
            rng.shuffle(idx)
            idx = idx[:self.config.maxPoints]
            r, theta, x, y = r0[idx], theta0[idx], x0[idx], y0[idx]  # noqa: r unused, but keeping consistent

        # improve the r,theta locations
        rNew, thetaNew, _r, _x, _y = self.improveCluster(theta, x, y)
        for i in range(self.config.nIter):
            rNew, thetaNew, _r, _x, _y = self.improveCluster(thetaNew, x, y)

        if rMax is None:
            rMax = np.sqrt(x.max()**2 + y.max()**2)

        # bin the data in r,theta space; get r,theta that pass our threshold as a trail
        bin2d, rEdge, thetaEdge, radii, thetas, indices = self.hesseBin(rNew, thetaNew, rMax=rMax)
        self.log.info("Detected %d peaks in Hough image", len(thetas))

        trails = TrailList()
        for radius, angle, index in zip(radii, thetas, indices):
            residual = x[index]*np.cos(angle) + y[index]*np.sin(angle) - radius

            q1  = np.percentile(residual, 25.0)
            q3  = np.percentile(residual, 75.0)
            iqr = q3 - q1

            # see if there's a significant 2nd-order term in a polynomial fit.
            #
            # Unfortunately, this isn't useful. I'll leave it in with a 'sanity' threshold
            # but to do it right, we need to consider the length of the trail, and the binning
            # used in the image ... maybe more.  Right now, it's rather scale dependent.

            # Ignore potential outliers
            cut = 10.0  # to clip points below 10% and above 90%
            resLo = np.percentile(residual, cut)
            resHi = np.percentile(residual, 100.0 - cut)
            wIqr = (residual > resLo) & (residual < resHi)
            poly = np.polyfit(np.arange(wIqr.sum()), residual[wIqr], 2)
            if not (iqr < self.config.maxResid) & (np.abs(poly[0]) < self.config.quadraticLimit):
                self.log.warn("Rejecting quadratic solution: r=%.1f,theta=%.3f  "
                              "(IQR=%.2f [limit=%.2f]  2nd-order coeff = %.2g [limit=%.2g])" %
                              (radius, angle, iqr, self.config.maxResid, poly[0],
                               self.config.quadraticLimit))
                continue
            trails.append(radius*binning, angle*afwGeom.radians, residIqr=iqr)

        return trails


def hesseForm(thetaIn, x, y):
    """Convert theta, x, y   to Hesse normal form

    @param thetaIn     Local position angle in radians (-pi/2 < thetaIn < pi/2)
    @param x           Local x coordinate in pixels
    @param y           Local y coordinate in pixels

    @return r, theta   Hesse normal form:  r = x*cos(theta) + y*sin(theta)

    The thetaIn we are given is a local position angle wrt to the x-axis.  For an elongated
    shape at position x,y; thetaIn is aligned along the *long* axis.  However, the theta
    used in the Hesse normal form is the angle of the vector normal to the local long axis.
    Theta is therefore different from thetaIn by +/- pi/2.  The sign is determined by the sign
    of the y-intercept.

    The basic geometry is shown at e.g. https://en.wikipedia.org/wiki/Hesse_normal_form
    """

    # if the intercept is y > 0, then add pi/2; otherwise subtract pi/2
    intercept = y - x*np.tan(thetaIn)
    pihalves = np.sign(intercept)*0.5*np.pi
    theta = thetaIn + pihalves

    # now ... -pi < theta < pi ...  convert to 0..2pi range
    theta[(theta < 0.0)] += 2.0*np.pi
    r = x*np.cos(theta) + y*np.sin(theta)

    return r, theta


def angleCompare(theta1, theta2, tolerance):
    if np.abs(theta2 - theta1) < tolerance:
        # no wrap problem
        return True

    # maybe we wrapped, use the dot product to compare
    c1, c2 = np.cos(theta1), np.cos(theta2)
    s1, s2 = np.sin(theta1), np.sin(theta2)
    dot = np.clip(c1*c2 + s1*s2, -1.0, 1.0)
    acos = np.arccos(dot)
    return np.abs(acos) < tolerance


class Trail(object):
    """Hold parameters related to a trail.

    Parameters are stored in Hesse normal form (r, theta).  The trail also
    knows its width, and its flux, and provides method to insert itself into an image
    or to set mask bits in an exposure.
    """

    def __init__(self, r, theta, residIqr=None):
        """Construct a Trail with specified parameters.

        @param r        r from Hesse normal form of the trail
        @param theta    theta from Hesse normal form of the trail
        @param residIqr    The coordinate inter_quart_range for residuals from the solution
        """

        self.r = r
        self.theta = theta
        self.residIqr = residIqr

    @property
    def theta(self):
        return self._theta

    @theta.setter
    def theta(self, value):
        self._theta = value
        radians = value.asRadians()
        self.vx = np.cos(radians)
        self.vy = np.sin(radians)

    def setMask(self, exposure, width, maskVal):
        """Set the mask plane near this trail in an exposure.

        @param exposure      The exposure with mask plane to be set.  Change is in-situ.
        @param maskVal       The bit to use for this mask.

        @return nPixels    The number of pixels set.
        """

        mask = exposure.getMaskedImage().getMask()

        # create a fresh mask and add to that.
        target = type(mask)(mask.getWidth(), mask.getHeight())
        profile = ConstantProfile(maskVal, width)
        self.insert(target, profile, width)

        mask |= target
        return (target.getArray() > 0).sum()

    def endPoints(self, nx, ny):
        """Compute the length of the trail in an nx*ny image.

        @param nx      Image x dimension
        @param ny      Image y dimension
        """

        points = []
        epsilon = 1.0e-8
        if np.abs(self.vy) > epsilon:
            for ix in 0, nx:
                y = (self.r - ix*self.vx)/self.vy
                if (y >= 0) and (y < ny):
                    points.append((ix, y))
        if np.abs(self.vx) > epsilon:
            for iy in 0, ny:
                x = (self.r - iy*self.vy)/self.vx
                if (x >= 0) and (x < nx):
                    points.append((x, iy))
        points = sorted(points, key=lambda x: x[0])
        if len(points) != 2:
            points = ((0, 0), (0, 0))
        return points

    def length(self, nx, ny):
        """Compute the length of the trail in an nx*ny image.

        @param nx      Image x dimension
        @param ny      Image y dimension
        """

        p1, p2 = self.endPoints(nx, ny)
        length = ((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)**0.5
        return length

    def trace(self, nx, ny, offset=0, bins=1):
        """Get x,y values near this trail.

        @param nx      Image width in pixels
        @param ny      Image height in pixels
        @param offset  Distance from trail centerline to return values
        @param bins    Correct for images binned by this amount.
        """

        p1, p2 = self.endPoints(nx, ny)
        xs = p1[0], p2[0]
        ys = p1[1], p2[1]

        dx = max(xs) - min(xs)
        dy = max(ys) - min(ys)
        if dx > dy:
            lo, hi = min(xs), max(xs)
            x = np.arange(lo, hi)
            y = (self.r/bins + offset - x*self.vx)/self.vy
        else:
            lo, hi = min(ys), max(ys)
            y = np.arange(lo, hi)
            x = (self.r/bins + offset - y*self.vy)/self.vx

        return np.rint(x).astype(int), np.rint(y).astype(int)


    def shiftOrigin(self, dx, dy):
        dr       = dx*self.vx + dy*self.vy
        rNew     = self.r - dr
        thetaNew = self.theta.asRadians()
        if rNew < 0.0:
            rNew *= -1.0
            thetaNew += np.pi if thetaNew < np.pi else -np.pi
        trail = copy.copy(self)
        trail.r = rNew
        trail.theta = thetaNew
        return trail

    def residual(self, x, y, bins=1):
        """Get residuals of this fit compared to given x,y coords.

        @param x   array of x pixel coord
        @param y   array of y pixel coord
        """

        dr = x*self.vx + y*self.vy - self.r/bins
        return dr


    def insert(self, exposure, profile, width):
        """Plant this trail in a given exposure.

        @param exposure       The exposure to plant in (accepts ExposureF, ImageF, MaskU or ndarray)
        @param profile        A profile function object.
        @param width          Set pixels to this value.  (Don't plant a Double-Gaussian trail profile).

        This method serves a few purposes.

        (1) To search for a trail with profile similar to a PSF, we plant a PSF-shaped trail
            and measure its parameters for use in CALIBRATING detection limits.
        (2) When we find a trail, our setMask() method calls this method with a maskBit to SET THE MASK.
        (3) For TESTING, we can insert fake trails and try to find them.
        """

        # Handle Exposure, Image, ndarray
        if isinstance(exposure, afwImage.ExposureF):
            image = exposure.getMaskedImage().getImage().getArray()
        elif isinstance(exposure, afwImage.ImageF):
            image = exposure.getArray()
        elif isinstance(exposure, afwImage.MaskU):
            image = exposure.getArray()
        elif isinstance(exposure, np.ndarray):
            image = exposure

        ny, nx = image.shape

        #############################
        # plant the trail
        #############################
        xpart, ypart = np.arange(nx)*self.vx, np.arange(ny)*self.vy
        dot = np.add.outer(ypart, xpart)

        # plant the trail using the distance from our line
        # as the parameter in a 1D DoubleGaussian
        offset = np.abs(dot - self.r)

        # only bother updating the pixels within width
        w = (offset < width/2.0)
        image[w] += profile(offset[w]).astype(image.dtype)
        return image


    def measure(self, exposure, bins, aperture):
        """Measure an aperture flux, a centroid, and a width for this trail in a given exposure.

        @param exposure       The exposure to measure in (accepts ExposureF, ImageF, ndarray)
        @param bins           The binning used in the given exposure.  Needed as r is in pixels.
        @param aperture          The aperture within which to measure.
        """
        # Handle Exposure, Image, ndarray
        if isinstance(exposure, afwImage.ExposureF):
            image = exposure.getMaskedImage().getImage().getArray()
            nx, ny = exposure.getWidth(), exposure.getHeight()
        elif isinstance(exposure, afwImage.ImageF):
            image = exposure.getArray()
            nx, ny = exposure.getWidth(), exposure.getHeight()
        elif isinstance(exposure, np.ndarray):
            image = exposure
            ny, nx = image.shape

        #############################
        # plant the trail
        #############################
        xx, yy = np.meshgrid(np.arange(nx), np.arange(ny))

        # plant the trail using the distance from our line
        # as the parameter in a 1D DoubleGaussian
        dot = xx*self.vx + yy*self.vy
        offset = dot - self.r/bins
        halfWidth = aperture/2.0

        w = (np.abs(offset) < halfWidth) & (np.isfinite(image))
        offset2 = offset[w]**2
        flux = image[w].sum()
        center = bins*(image[w]*offset[w]).sum()/flux
        sigma = np.sqrt((image[w]*offset2).sum()/flux)

        # to suppress noise, iterate with a weight based on first measurement
        weight = np.exp(-0.5*offset2/sigma**2)
        sigma = bins*np.sqrt(((weight*image[w]*offset2).sum()/(weight*image[w]).sum()))

        return Struct(flux=flux, center=center, width=sigma)

    def __str__(self):
        return "%s(r=%.1f,theta=%.3f deg,residIqr=%.2f)" % (self.__class__.__name__, self.r,
                                                        self.theta.asDegrees(), self.residIqr)

    def __repr__(self):
        return "%s(r=%r,theta=%r,residIqr=%r)" % (self.__class__.__name__, self.r, self.theta, self.residIqr)

    def __eq__(self, trail):
        return (self.r == trail.r and self.theta == trail.theta and
                self.width == trail.width and self.flux == trail.flux)

    def isNear(self, trail, drMax, dThetaMax):
        """Fuzzy-compare two trails.

        It's quite possible that the same trail will be detected in a searches for satellites and
        aircraft.  The parameters won't be identical, but they'll be close.

        Use dot product to handle the wrap at 2*pi
        """
        r1, r2 = self.r,     trail.r
        t1, t2 = self.theta, trail.theta
        rTest = np.abs(r2 - r1) < drMax
        tTest = angleCompare(t1, t2, dThetaMax)
        return rTest and tTest

    @staticmethod
    def chooseBest(trail1, trail2):
        """A single place to choose the best trail, if two solutions exist.

        @param trail1   Trail object #1
        @param trail2   Trail object #2
        """
        err1 = trail1.residIqr
        err2 = trail2.residIqr
        return trail1 if (err1 < err2) else trail2

    def display(self, dimensions, frame=1, ctype=None, size=0.5):
        import lsst.afw.display
        lsst.afw.display.getDisplay(frame).line(self.endPoints(*dimensions), origin=afwImage.LOCAL,
                                                ctype=ctype, size=size)


class TrailList(list):
    """A container for Trail objects"""
    def cleanDuplicates(self, drMax, dThetaMax):
        """Go through this list and remove duplicates.  Return a new TrailList."""

        # clean duplicates from each list
        newTrailList = TrailList()
        skip = set()
        for ii, trail1 in enumerate(self):
            if ii in skip:
                continue
            best = trail1
            for jj, trail2 in enumerate(self[ii:]):
                if trail1.isNear(trail2, drMax, dThetaMax):
                    best = trail1.chooseBest(best, trail2)
                    skip.add(jj)
            list.append(newTrailList, best)
        return newTrailList

    def merge(self, trailList, drMax=90.0, dThetaMax=0.15):
        """Merge trails from trailList to this TrailList.  Returns a new TrailList.

        @param trailList     The trailList to merge in
        @param drMax         The max separation in r for identifying duplicates
        @param dthetaMax     The max separation in theta for identifying duplicates
        """

        newList = TrailList()
        tList1 = trailList.cleanDuplicates(drMax, dThetaMax)
        tList2 = self.cleanDuplicates(drMax, dThetaMax)

        # get everything from list 1, and check for duplicates
        for trail1 in tList1:
            best  = trail1
            for trail2 in tList2:
                if trail1.isNear(trail2, drMax, dThetaMax):
                    best = trail1.chooseBest(trail1, trail2)
            list.append(newList, best)

        # get everything from list 2, and throw out duplicates (we already chose the best one)
        for trail2 in tList2:
            haveIt = [trail2.isNear(tNew, drMax, dThetaMax) for tNew in newList]
            if not any(haveIt):
                list.append(newList, trail2)

        return newList

    def append(self, r, theta, residIqr=None):
        """Append a Trail

        @param r        r from Hesse normal form of the trail
        @param theta    theta from Hesse normal form of the trail
        @param residIqr    The coordinate inter_quart_range for residuals from the solution
        """
        list.append(self, Trail(r, theta, residIqr=residIqr))

    def display(self, *args, **kwargs):
        for trail in self:
            trail.display(*args, **kwargs)


class ConstantProfile(object):
    """A constant trail profile"""

    def __init__(self, value, width):
        """Construct

        @param value  The constant value to set
        @param width  The width of the trail
        """
        self.value = value
        self.width = width

    def __call__(self, offset):
        """Return profile value at 'offset'

        @param offset  A distance from the center.
        """
        w  = (offset <= self.width/2.0)
        return self.value*w


class DoubleGaussianProfile(object):
    """A Double Gaussian trail profile"""

    def __init__(self, flux, sigma, fWing=0.1):
        """Construct

        @param flux       The flux of the trail
        @param sigma      Sigma of the inner Gaussian
        @param fWing      Fraction of flux in the wing (outer) Gaussian.
        """
        self.flux  = flux
        self.sigma = sigma
        self.fWing = fWing
        self.fCore = 1.0 - fWing
        self.A1 = 1.0/(np.sqrt(2.0*np.pi)*self.sigma)
        self.A2 = 1.0/(np.sqrt(2.0*np.pi)*(2.0*self.sigma))
        self.coef1 = self.flux*self.fCore*self.A1
        self.coef2 = self.flux*self.fWing*self.A2
        self.twoSigma2Core = 2.0*self.sigma**2
        self.twoSigma2Wing = 2.0*(2.0*self.sigma)**2

    def __call__(self, offset):
        """Return profile value at offset.

        @param offset   Distance from the profile center.
        """
        g1  = np.exp(-offset**2/self.twoSigma2Core)
        g2  = np.exp(-offset**2/self.twoSigma2Wing)
        out = self.coef1*g1 + self.coef2*g2
        return out
