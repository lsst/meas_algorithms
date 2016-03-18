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

from lsst.afw.cameraGeom import TAN_PIXELS
from lsst.afw.geom.ellipses import Quadrupole
from lsst.afw.table import SourceCatalog, SourceTable
from lsst.pipe.base import Struct
import lsst.pex.config as pexConfig
import lsst.afw.detection as afwDetection
import lsst.afw.display.ds9 as ds9
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
from . import algorithmsLib
from lsst.meas.base import SingleFrameMeasurementTask, SingleFrameMeasurementConfig
from .starSelector import StarSelectorTask, starSelectorRegistry

class SecondMomentStarSelectorConfig(StarSelectorTask.ConfigClass):
    fluxLim = pexConfig.Field(
        doc = "specify the minimum psfFlux for good Psf Candidates",
        dtype = float,
        default = 12500.0,
        check = lambda x: x >= 0.0,
    )
    fluxMax = pexConfig.Field(
        doc = "specify the maximum psfFlux for good Psf Candidates (ignored if == 0)",
        dtype = float,
        default = 0.0,
        check = lambda x: x >= 0.0,
    )
    clumpNSigma = pexConfig.Field(
        doc = "candidate PSF's shapes must lie within this many sigma of the average shape",
        dtype = float,
        default = 2.0,
        check = lambda x: x >= 0.0,
    )
    histSize = pexConfig.Field(
        doc = "Number of bins in moment histogram",
        dtype = int,
        default = 64,
        check = lambda x: x > 0,
        )
    histMomentMax = pexConfig.Field(
        doc = "Maximum moment to consider",
        dtype = float,
        default = 100.0,
        check = lambda x: x > 0,
        )
    histMomentMaxMultiplier = pexConfig.Field(
        doc = "Multiplier of mean for maximum moments histogram range",
        dtype = float,
        default = 5.0,
        check = lambda x: x > 0,
        )
    histMomentClip = pexConfig.Field(
        doc = "Clipping threshold for moments histogram range",
        dtype = float,
        default = 5.0,
        check = lambda x: x > 0,
        )
    histMomentMinMultiplier = pexConfig.Field(
        doc = "Multiplier of mean for minimum moments histogram range",
        dtype = float,
        default = 2.0,
        check = lambda x: x > 0,
        )

    def setDefaults(self):
        StarSelectorTask.ConfigClass.setDefaults(self)
        self.badFlags = [
            "base_PixelFlags_flag_edge",
            "base_PixelFlags_flag_interpolatedCenter",
            "base_PixelFlags_flag_saturatedCenter",
            "base_PixelFlags_flag_crCenter",
        ]

Clump = collections.namedtuple('Clump', ['peak', 'x', 'y', 'ixx', 'ixy', 'iyy', 'a', 'b', 'c'])

class CheckSource(object):
    """A functor to check whether a source has any flags set that should cause it to be labeled bad."""

    def __init__(self, table, badFlags, fluxLim, fluxMax):
        self.keys = [table.getSchema().find(name).key for name in badFlags]
        self.keys.append(table.getCentroidFlagKey())
        self.fluxLim = fluxLim
        self.fluxMax = fluxMax

    def __call__(self, source):
        for k in self.keys:
            if source.get(k):
                return False
        if self.fluxLim != None and source.getPsfFlux() < self.fluxLim: # ignore faint objects
            return False
        if self.fluxMax != 0.0 and source.getPsfFlux() > self.fluxMax: # ignore bright objects
            return False
        return True

class SecondMomentStarSelectorTask(StarSelectorTask):
    """!A star selector based on second moments

    @warning This is a naive algorithm and should be used with caution
    """
    ConfigClass = SecondMomentStarSelectorConfig
    usesMatches = False # selectStars does not use its matches argument

    def selectStars(self, exposure, sourceCat, matches=None):
        """!Return a list of PSF candidates that represent likely stars
        
        A list of PSF candidates may be used by a PSF fitter to construct a PSF.
        
        @param[in] exposure  the exposure containing the sources
        @param[in] sourceCat  catalog of sources that may be stars (an lsst.afw.table.SourceCatalog)
        @param[in] matches  astrometric matches; ignored by this star selector
        
        @return an lsst.pipe.base.Struct containing:
        - starCat  catalog of selected stars (a subset of sourceCat)
        """
        import lsstDebug
        display = lsstDebug.Info(__name__).display

        displayExposure = lsstDebug.Info(__name__).displayExposure     # display the Exposure + spatialCells
        
        isGoodSource = CheckSource(sourceCat.getTable(), self.config.badFlags, self.config.fluxLim,
                                   self.config.fluxMax)

        detector = exposure.getDetector()

        mi = exposure.getMaskedImage()
        #
        # Create an Image of Ixx v. Iyy, i.e. a 2-D histogram
        #

        # Use stats on our Ixx/yy values to determine the xMax/yMax range for clump image
        iqqList = []
        for s in sourceCat:
            ixx, iyy = s.getIxx(), s.getIyy()
            # ignore NaN and unrealistically large values
            if (ixx == ixx and ixx < self.config.histMomentMax and
                iyy == iyy and iyy < self.config.histMomentMax and
                isGoodSource(s)):
                iqqList.append(s.getIxx())
                iqqList.append(s.getIyy())
        stat = afwMath.makeStatistics(iqqList, afwMath.MEANCLIP | afwMath.STDEVCLIP | afwMath.MAX)
        iqqMean = stat.getValue(afwMath.MEANCLIP)
        iqqStd = stat.getValue(afwMath.STDEVCLIP)
        iqqMax = stat.getValue(afwMath.MAX)

        iqqLimit = max(iqqMean + self.config.histMomentClip*iqqStd,
                       self.config.histMomentMaxMultiplier*iqqMean)
        # if the max value is smaller than our range, use max as the limit, but don't go below N*mean
        if iqqLimit > iqqMax:
            iqqLimit = max(self.config.histMomentMinMultiplier*iqqMean, iqqMax)

        psfHist = _PsfShapeHistogram(detector=detector, xSize=self.config.histSize, ySize=self.config.histSize,
                                     ixxMax=iqqLimit, iyyMax=iqqLimit)

        if display and displayExposure:
            frame = 0
            ds9.mtv(mi, frame=frame, title="PSF candidates")
    
        with ds9.Buffering():
            for source in sourceCat:
                if isGoodSource(source):
                    if psfHist.insert(source): # n.b. this call has the side effect of inserting
                         ctype = ds9.GREEN # good
                    else:
                         ctype = ds9.MAGENTA # rejected
                else:
                    ctype = ds9.RED         # bad

                if display and displayExposure:
                    ds9.dot("o", source.getX() - mi.getX0(),
                            source.getY() - mi.getY0(), frame=frame, ctype=ctype)

        clumps = psfHist.getClumps(display=display)

        #
        # Go through and find all the PSF-like objects
        #
        # We'll split the image into a number of cells, each of which contributes only
        # one PSF candidate star
        #
        starCat = SourceCatalog(sourceCat.schema)

        pixToTanXYTransform = None
        if detector is not None:
            tanSys = detector.makeCameraSys(TAN_PIXELS)
            pixToTanXYTransform = detector.getTransformMap().get(tanSys)
    
        # psf candidate shapes must lie within this many RMS of the average shape
        # N.b. if Ixx == Iyy, Ixy = 0 the criterion is
        # dx^2 + dy^2 < self.config.clumpNSigma*(Ixx + Iyy) == 2*self.config.clumpNSigma*Ixx
        for source in sourceCat:
            if not isGoodSource(source): continue
            Ixx, Ixy, Iyy = source.getIxx(), source.getIxy(), source.getIyy()
            if pixToTanXYTransform:
                p = afwGeom.Point2D(source.getX(), source.getY())
                linTransform = pixToTanXYTransform.linearizeForwardTransform(p).getLinear()
                m = Quadrupole(Ixx, Iyy, Ixy)
                m.transform(linTransform)
                Ixx, Iyy, Ixy = m.getIxx(), m.getIyy(), m.getIxy()
            
            x, y = psfHist.momentsToPixel(Ixx, Iyy)
            for clump in clumps:
                dx, dy = (x - clump.x), (y - clump.y)

                if math.sqrt(clump.a*dx*dx + 2*clump.b*dx*dy + clump.c*dy*dy) < 2*self.config.clumpNSigma:
                    # A test for > would be confused by NaN
                    if not isGoodSource(source):
                        continue
                    try:
                        psfCandidate = algorithmsLib.makePsfCandidate(source, exposure)
                        
                        # The setXXX methods are class static, but it's convenient to call them on
                        # an instance as we don't know Exposure's pixel type
                        # (and hence psfCandidate's exact type)
                        if psfCandidate.getWidth() == 0:
                            psfCandidate.setBorderWidth(self.config.borderWidth)
                            psfCandidate.setWidth(self.config.kernelSize + 2*self.config.borderWidth)
                            psfCandidate.setHeight(self.config.kernelSize + 2*self.config.borderWidth)

                        im = psfCandidate.getMaskedImage().getImage()
                        if not numpy.isfinite(afwMath.makeStatistics(im, afwMath.MAX).getValue()):
                            continue
                        starCat.append(source)

                        if display and displayExposure:
                            ds9.dot("o", source.getX() - mi.getX0(), source.getY() - mi.getY0(),
                                    size=4, frame=frame, ctype=ds9.CYAN)
                    except Exception as err:
                        self.log.error("Failed on source %s: %s" % (source.getId(), err))
                    break

        return Struct(
            starCat = starCat,
        )

class _PsfShapeHistogram(object):
    """A class to represent a histogram of (Ixx, Iyy)
    """
    def __init__(self, xSize=32, ySize=32, ixxMax=30, iyyMax=30, detector=None, xy0=afwGeom.Point2D(0,0)):
        """Construct a _PsfShapeHistogram

        The maximum seeing FWHM that can be tolerated is [xy]Max/2.35 pixels.
        The 'resolution' of stars vs galaxies/CRs is provided by [xy]Size/[xy]Max.
        A larger (better) resolution may thresh the peaks, but a smaller (worse)
        resolution will allow stars and galaxies/CRs to mix.  The disadvantages of
        a larger (better) resolution can be compensated (some) by using multiple
        histogram peaks.
        
        @input[in] [xy]Size: the size of the psfImage (in pixels)
        @input[in] ixxMax, iyyMax: the maximum values for I[xy][xy]
        """
        self._xSize, self._ySize = xSize, ySize 
        self._xMax, self._yMax = ixxMax, iyyMax
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
            tanSys = self.detector.makeCameraSys(TAN_PIXELS)
            if tanSys in self.detector.getTransformMap():
                pixToTanXYTransform = self.detector.getTransformMap()[tanSys]
                p = afwGeom.Point2D(source.getX(), source.getY())
                linTransform = pixToTanXYTransform.linearizeForwardTransform(p).getLinear()
                m = Quadrupole(ixx, iyy, ixy)
                m.transform(linTransform)
                ixx, iyy, ixy = m.getIxx(), m.getIyy(), m.getIxy()
            
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
        largeImg.assign(psfImage, bbox, afwImage.LOCAL) 
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
            
        ds = afwDetection.FootprintSet(mpsfImage, threshold, "DETECTED")
        #
        # And measure it.  This policy isn't the one we use to measure
        # Sources, it's only used to characterize this PSF histogram
        #
        schema = SourceTable.makeMinimalSchema()
        psfImageConfig = SingleFrameMeasurementConfig()
        psfImageConfig.doApplyApCorr = "no"
        psfImageConfig.slots.centroid = "base_SdssCentroid"
        psfImageConfig.slots.psfFlux = None #"base_PsfFlux"
        psfImageConfig.slots.apFlux = "base_CircularApertureFlux_3_0"
        psfImageConfig.slots.modelFlux = None
        psfImageConfig.slots.instFlux = None
        psfImageConfig.slots.calibFlux = None
        psfImageConfig.slots.shape = "base_SdssShape"
        #   Formerly, this code had centroid.sdss, flux.psf, flux.naive,
        #   flags.pixel, and shape.sdss
        psfImageConfig.algorithms.names = ["base_SdssCentroid", "base_CircularApertureFlux", "base_SdssShape"]
        psfImageConfig.algorithms["base_CircularApertureFlux"].radii = [3.0]
        psfImageConfig.validate()
        task = SingleFrameMeasurementTask(schema, config=psfImageConfig)

        sourceCat = SourceCatalog(schema)

        gaussianWidth = 1.5                       # Gaussian sigma for detection convolution
        exposure.setPsf(algorithmsLib.DoubleGaussianPsf(11, 11, gaussianWidth))

        ds.makeSources(sourceCat)
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
        task.run(exposure, sourceCat)   # notes that this is backwards for the new framework
        for i, source in enumerate(sourceCat):
            if source.getCentroidFlag():
                continue
            x, y = source.getX(), source.getY()

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
        # use the apFlux from the clump measurement, and take the highest
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

starSelectorRegistry.register("secondMoment", SecondMomentStarSelectorTask)
