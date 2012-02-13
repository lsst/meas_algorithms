# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import collections
import math

import numpy

import lsst.pex.config as pexConfig
import lsst.afw.detection as afwDetection
import lsst.afw.display.ds9 as ds9
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses as geomEllip
import lsst.afw.cameraGeom as cameraGeom
from . import algorithmsLib
from . import measurement

class SecondMomentStarSelectorConfig(pexConfig.Config):
    fluxLim = pexConfig.Field(
        doc = "specify the minimum psfFlux for good Psf Candidates",
        dtype = float,
        default = 12500.0,
#        minValue = 0.0,
        check = lambda x: x >= 0.0,
    )
    fluxMax = pexConfig.Field(
        doc = "specify the maximum psfFlux for good Psf Candidates (ignored if == 0)",
        dtype = float,
        default = 0.0,
#        minValue = 0.0,
        check = lambda x: x >= 0.0,
    )
    clumpNSigma = pexConfig.Field(
        doc = "candidate PSF's shapes must lie within this many sigma of the average shape",
        dtype = float,
        default = 1.0,
#        minValue = 0.0,
        check = lambda x: x >= 0.0,
    )
    kernelSize = pexConfig.Field(
        doc = "size of the kernel to create",
        dtype = int,
        default = 21,
    )
    borderWidth = pexConfig.Field(
        doc = "number of pixels to ignore around the edge of PSF candidate postage stamps",
        dtype = int,
        default = 0,
    )


Clump = collections.namedtuple('Clump', ['peak', 'x', 'y', 'ixx', 'ixy', 'iyy', 'a', 'b', 'c'])

class SecondMomentStarSelector(object):
    ConfigClass = SecondMomentStarSelectorConfig
    
    _badSourceMask = algorithmsLib.Flags.EDGE | \
        algorithmsLib.Flags.INTERP_CENTER | \
        algorithmsLib.Flags.SATUR_CENTER | \
        algorithmsLib.Flags.PEAKCENTER

    def __init__(self, config):
        """Construct a star selector that uses second moments
        
        This is a naive algorithm and should be used with caution.
        
        @param[in] config: an instance of SecondMomentStarSelectorConfig
        """
        self._kernelSize  = config.kernelSize
        self._borderWidth = config.borderWidth
        self._clumpNSigma = config.clumpNSigma
        self._fluxLim  = config.fluxLim
        self._fluxMax  = config.fluxMax
    

    def selectStars(self, exposure, sourceList):
        """Return a list of PSF candidates that represent likely stars
        
        A list of PSF candidates may be used by a PSF fitter to construct a PSF.
        
        @param[in] exposure: the exposure containing the sources
        @param[in] sourceList: a list of Sources that may be stars
        
        @return psfCandidateList: a list of PSF candidates.
        """
        import lsstDebug
        display = lsstDebug.Info(__name__).display
        displayExposure = lsstDebug.Info(__name__).displayExposure     # display the Exposure + spatialCells

	detector = exposure.getDetector()
	distorter = None
	xy0 = afwGeom.Point2D(0,0)
	if not detector is None:
	    cPix = detector.getCenterPixel()
	    detSize = detector.getSize()
	    xy0.setX(cPix.getX() - int(0.5*detSize.getMm()[0]))
	    xy0.setY(cPix.getY() - int(0.5*detSize.getMm()[1]))
	    distorter = detector.getDistortion()

        mi = exposure.getMaskedImage()
        #
        # Create an Image of Ixx v. Iyy, i.e. a 2-D histogram
        #

	# Use stats on our Ixx/yy values to determine the xMax/yMax range for clump image
	iqqList = []
	for s in sourceList:
	    ixx, iyy = s.getIxx(), s.getIyy()
            # ignore NaN and unrealistically large values
	    if ixx == ixx and ixx < 100.0 and iyy == iyy and iyy < 100.0:
		iqqList.append(s.getIxx())
		iqqList.append(s.getIyy())
        stat = afwMath.makeStatistics(iqqList, afwMath.MEANCLIP | afwMath.STDEVCLIP | afwMath.MAX)
	iqqMean = stat.getValue(afwMath.MEANCLIP)
	iqqStd = stat.getValue(afwMath.STDEVCLIP)
        iqqMax = stat.getValue(afwMath.MAX)
        
        # set the limits to run as high as 5sigma, but rail at iqq=20 (that's fwhm=9pixel ... very high)
	iqqLimit = numpy.max([iqqMean + 5.0*iqqStd, 5.0*iqqMean])
        # if the max value is smaller than our range, use max as the limit, but don't go below 2*mean
        if iqqLimit > iqqMax:
            iqqLimit = numpy.max([2.0*iqqMean, iqqMax])
            
        psfHist = _PsfShapeHistogram(detector=detector, xMax=iqqLimit, yMax=iqqLimit, xy0=xy0)
	
        if display and displayExposure:
            frame = 0
            ds9.mtv(mi, frame=frame, title="PSF candidates")
    
        with ds9.Buffering():
            for source in sourceList:
                if self._isGoodSource(source):
                    if psfHist.insert(source): # n.b. this call has the side effect of inserting
                         ctype = ds9.GREEN # good
                    else:
                         ctype = ds9.MAGENTA # rejected
                else:
                    ctype = ds9.RED         # bad

                if display and displayExposure:
                    ds9.dot("o", source.getXAstrom() - mi.getX0(),
                            source.getYAstrom() - mi.getY0(), frame=frame, ctype=ctype)

        clumps = psfHist.getClumps(display=display)
        
        #
        # Go through and find all the PSF-like objects
        #
        # We'll split the image into a number of cells, each of which contributes only
        # one PSF candidate star
        #
        psfCandidateList = []
    
        # psf candidate shapes must lie within this many RMS of the average shape
        # N.b. if Ixx == Iyy, Ixy = 0 the criterion is
        # dx^2 + dy^2 < self._clumpNSigma*(Ixx + Iyy) == 2*self._clumpNSigma*Ixx
        for source in sourceList:
	    
            Ixx, Ixy, Iyy = source.getIxx(), source.getIxy(), source.getIyy()
	    if distorter:
		xpix, ypix = source.getXAstrom() + xy0.getX(), source.getYAstrom() + xy0.getY()
		p = afwGeom.Point2D(xpix, ypix)
		m = distorter.undistort(p, geomEllip.Quadrupole(Ixx, Iyy, Ixy), detector)
		Ixx, Iyy, Ixy = m.getIXX(), m.getIYY(), m.getIXY()
	    
            x, y = psfHist.momentsToPixel(Ixx, Iyy)
            for clump in clumps:
                dx, dy = (x - clump.x), (y - clump.y)

                if math.sqrt(clump.a*dx*dx + 2*clump.b*dx*dy + clump.c*dy*dy) < 2*self._clumpNSigma:
                    # A test for > would be confused by NaN
                    if not self._isGoodSource(source):
                        continue
                    try:
                        psfCandidate = algorithmsLib.makePsfCandidate(source, exposure)
			
                        # The setXXX methods are class static, but it's convenient to call them on
                        # an instance as we don't know Exposure's pixel type
                        # (and hence psfCandidate's exact type)
                        if psfCandidate.getWidth() == 0:
                            psfCandidate.setBorderWidth(self._borderWidth)
                            psfCandidate.setWidth(self._kernelSize + 2*self._borderWidth)
                            psfCandidate.setHeight(self._kernelSize + 2*self._borderWidth)

                        im = psfCandidate.getImage().getImage()
                        max = afwMath.makeStatistics(im, afwMath.MAX).getValue()
                        if not numpy.isfinite(max):
                            continue

                        source.setFlagForDetection(source.getFlagForDetection() | algorithmsLib.Flags.STAR)
                        psfCandidateList.append(psfCandidate)

                        if display and displayExposure:
                            ds9.dot("o", source.getXAstrom() - mi.getX0(), source.getYAstrom() - mi.getY0(),
                                    size=4, frame=frame, ctype=ds9.CYAN)
                    except:
                        pass
                    break

        return psfCandidateList

    def _isGoodSource(self, source):
        """Should this object be included in the Ixx v. Iyy image?
        """ 

        if source.getFlagForDetection() & self._badSourceMask:
            return False
        if self._fluxLim != None and source.getPsfFlux() < self._fluxLim: # ignore faint objects
            return False
        if self._fluxMax != 0.0 and source.getPsfFlux() > self._fluxMax: # ignore bright objects
            return False

        return True


class _PsfShapeHistogram(object):
    """A class to represent a histogram of (Ixx, Iyy)
    """
    def __init__(self, xSize=32, ySize=32, xMax=30, yMax=30, detector=None, xy0=afwGeom.Point2D(0,0)):
        """Construct a _PsfShapeHistogram

        The maximum seeing FWHM that can be tolerated is [xy]Max/2.35 pixels.
        The 'resolution' of stars vs galaxies/CRs is provided by [xy]Size/[xy]Max.
        A larger (better) resolution may thresh the peaks, but a smaller (worse)
        resolution will allow stars and galaxies/CRs to mix.  The disadvantages of
        a larger (better) resolution can be compensated (some) by using multiple
        histogram peaks.
        
        @input[in] [xy]Size: the size of the psfImage (in pixels)
        @input[in] [xy]Max: the maximum values for I[xy][xy]
        """
        self._xSize, self._ySize = xSize, ySize 
        self._xMax, self._yMax = xMax, yMax
        self._psfImage = afwImage.ImageF(afwGeom.ExtentI(xSize, ySize), 0)
        self._num = 0
	self.detector = detector
	self.xy0 = xy0

    def getImage(self):
        return self._psfImage

    def insert(self, source):
        """Insert source into the histogram."""
	
	ixx, iyy, ixy = source.getIxx(), source.getIyy(), source.getIxy()
	if self.detector:
            distorter = self.detector.getDistortion()
            if distorter:
                p = afwGeom.Point2D(source.getXAstrom()+self.xy0.getX(),
                                    source.getYAstrom() + self.xy0.getY())
                m = distorter.undistort(p, geomEllip.Quadrupole(ixx, iyy, ixy), self.detector)
                ixx, iyy, ixy = m.getIXX(), m.getIYY(), m.getIXY()
	    
        try:
            pixel = self.momentsToPixel(ixx, iyy)
            i = int(pixel[0])
            j = int(pixel[1])
        except:
            return 0

        if i in range(0, self._xSize) and j in range(0, self._ySize):
            if i != 0 or j != 0:
                self._psfImage.set(i, j, self._psfImage.get(i, j) + 1)
                self._num += 1
                return 1                # success

        return 0                        # failure

    def momentsToPixel(self, ixx, iyy):
        #x = math.sqrt(ixx) * self._xSize / self._xMax
        #y = math.sqrt(iyy) * self._ySize / self._yMax
        x = ixx * self._xSize / self._xMax
        y = iyy * self._ySize / self._yMax
        return x, y

    def pixelToMoments(self, x, y):
        """Given a peak position in self._psfImage, return the corresponding (Ixx, Iyy)"""

        #ixx = (x*self._xMax/self._xSize)**2
        #iyy = (y*self._yMax/self._ySize)**2
        ixx = x*self._xMax/self._xSize
        iyy = y*self._yMax/self._ySize
        return ixx, iyy

    def getClumps(self, sigma=1.0, display=False):
        if self._num <= 0:
            raise RuntimeError("No candidate PSF sources")

        psfImage = self.getImage()
        #
        # Embed psfImage into a larger image so we can smooth when measuring it
        #
        width, height = psfImage.getWidth(), psfImage.getHeight()
        largeImg = psfImage.Factory(afwGeom.ExtentI(2*width, 2*height))
        largeImg.set(0)

        bbox = afwGeom.BoxI(afwGeom.PointI(width, height), afwGeom.ExtentI(width, height))
        subLargeImg = psfImage.Factory(largeImg, bbox, afwImage.LOCAL)
        subLargeImg <<= psfImage
        del subLargeImg
        #
        # Now measure that image, looking for the highest peak.  Start by building an Exposure
        #
        msk = afwImage.MaskU(largeImg.getDimensions())
        msk.set(0)
        var = afwImage.ImageF(largeImg.getDimensions())
        var.set(1)
        mpsfImage = afwImage.MaskedImageF(largeImg, msk, var)
        mpsfImage.setXY0(afwGeom.PointI(-width, -height))
        del msk
        del var
        exposure = afwImage.makeExposure(mpsfImage)
        
        #
        # Next run an object detector
        #
        maxVal = afwMath.makeStatistics(psfImage, afwMath.MAX).getValue()
        threshold = maxVal - sigma*math.sqrt(maxVal)
        if threshold <= 0.0:
            threshold = maxVal

        threshold = afwDetection.Threshold(threshold)
            
        ds = afwDetection.FootprintSetF(mpsfImage, threshold, "DETECTED")
        footprints = ds.getFootprints()
        #
        # And measure it.  This policy isn't the one we use to measure
        # Sources, it's only used to characterize this PSF histogram
        #
        psfImageConfig = measurement.MeasureSourcesConfig()
        psfImageConfig.source.astrom = "SDSS"
        psfImageConfig.source.psfFlux = "PSF"
        psfImageConfig.source.apFlux = "NAIVE"
        psfImageConfig.source.modelFlux = None
        psfImageConfig.source.instFlux = None
        psfImageConfig.source.shape = "SDSS"
        psfImageConfig.astrometry.names = ["SDSS"]
        psfImageConfig.photometry.names = ["PSF", "NAIVE"]
        psfImageConfig.photometry["NAIVE"].radius = 3.0
        psfImageConfig.shape.names = ["SDSS"]
        
        gaussianWidth = 1.5                       # Gaussian sigma for detection convolution
        exposure.setPsf(afwDetection.createPsf("DoubleGaussian", 11, 11, gaussianWidth))
        measureSources = psfImageConfig.makeMeasureSources(exposure)
        
        #
        # Show us the Histogram
        #
        if display:
            frame = 1
            dispImage = mpsfImage.Factory(mpsfImage, afwGeom.BoxI(afwGeom.PointI(width, height),
                                                                  afwGeom.ExtentI(width, height)),
                                                                  afwImage.LOCAL)
            ds9.mtv(dispImage,title="PSF Selection Image", frame=frame)


        clumps = list()                 # List of clumps, to return
        e = None                        # thrown exception
        IzzMin = 1.0                    # Minimum value for second moments
        IzzMax = (self._xSize/8.0)**2   # Max value ... clump r < clumpImgSize/8
                                        # diameter should be < 1/4 clumpImgSize
        apFluxes = []
        for i, foot in enumerate(footprints):
            source = afwDetection.Source()
            source.setId(i)
            source.setFootprint(foot)

            try:
                measureSources.measure(source, exposure)
            except Exception, e:
                print "Except:", e
                continue

            x, y = source.getXAstrom(), source.getYAstrom()
            
            # in very distorted tests, some footprints return nan for x,y and can be ignored
            if x != x or y != y:
                continue

            apFluxes.append(source.getApFlux())
            
            val = mpsfImage.getImage().get(int(x) + width, int(y) + height)

            psfClumpIxx = source.getIxx()
            psfClumpIxy = source.getIxy()
            psfClumpIyy = source.getIyy()

            if display:
                if i == 0:
                    ds9.pan(x, y, frame=frame)

                ds9.dot("+", x, y, ctype=ds9.YELLOW, frame=frame)
                ds9.dot("@:%g,%g,%g" % (psfClumpIxx, psfClumpIxy, psfClumpIyy), x, y,
                        ctype=ds9.YELLOW, frame=frame)

            if psfClumpIxx < IzzMin or psfClumpIyy < IzzMin:
                psfClumpIxx = max(psfClumpIxx, IzzMin)
		#psfClumpIxy = 0.0
                psfClumpIyy = max(psfClumpIyy, IzzMin)
                if display:
                    ds9.dot("@:%g,%g,%g" % (psfClumpIxx, psfClumpIxy, psfClumpIyy), x, y,
                            ctype=ds9.RED, frame=frame)

            det = psfClumpIxx*psfClumpIyy - psfClumpIxy*psfClumpIxy
            try:
                a, b, c = psfClumpIyy/det, -psfClumpIxy/det, psfClumpIxx/det
            except ZeroDivisionError:
                a, b, c = 1e4, 0, 1e4

            clumps.append(Clump(peak=val, x=x, y=y, a=a, b=b, c=c,
                                ixx=psfClumpIxx, ixy=psfClumpIxy, iyy=psfClumpIyy))

        if len(clumps) == 0:
            msg = "Failed to determine center of PSF clump"
            if e:
                msg += ": %s" % e
            raise RuntimeError(msg)

        # if it's all we got return it
        if len(clumps) == 1:
            return clumps
        
        # which clump is the best?
        # if we've undistorted the moments, stars should only have 1 clump
        # use the apFlux from the clump measureSources, and take the highest
        # ... this clump has more psf star candidate neighbours than the others.

        # get rid of any that are huge, and thus poorly defined
        goodClumps = []
        for clump in clumps:
            if clump.ixx < IzzMax and clump.iyy < IzzMax:
                goodClumps.append(clump)

        # if culling > IzzMax cost us all clumps, we'll have to take what we have
        if len(goodClumps) == 0:
            goodClumps = clumps
            
        # use the 'brightest' clump
        iBestClump = numpy.argsort(apFluxes)[0]
        clumps = [clumps[iBestClump]]
        return clumps
