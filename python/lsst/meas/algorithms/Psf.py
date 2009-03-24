# This is not a minimal set of imports
import glob, math, os, sys
from math import *
import eups
import lsst.daf.base as dafBase
import lsst.pex.logging as logging
import lsst.pex.policy as policy
import lsst.afw.detection as afwDetection
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.meas.algorithms as algorithms
import lsst.meas.algorithms.defects as defects
import lsst.meas.algorithms.measureSourceUtils as measureSourceUtils
import lsst.sdqa as sdqa

import lsst.afw.display.ds9 as ds9
import lsst.afw.display.utils as displayUtils

try:
    type(display)
except NameError:
    display = False

class PsfShapeHistogram(object):
    """A class to represent a histogram of (Ixx, Iyy)"""

    def __init__(self, xSize=20, ySize=20, xMax=15, yMax=15):
        self._xSize, self._ySize = xSize, ySize 
        self._xMax, self._yMax = xMax, yMax
        self._psfImage = afwImage.ImageF(self._xSize, self._ySize)
        self._psfImage.set(0)

    def getImage(self):
        return self._psfImage

    def insert(self, source):
        """Insert source into the histogram."""
        i = int(source.getIxx()*self._xSize/self._xMax + 0.5)
        j = int(source.getIyy()*self._ySize/self._yMax + 0.5)
        if i in range(0, self._xSize) and j in range(0, self._ySize):
            if i != 0 or j != 0:
                self._psfImage.set(i, j, self._psfImage.get(i, j) + 1)

    def peakToIxx(self, peakX, peakY):
        """Given a peak position in self._psfImage, return the corresponding (Ixx, Iyy)"""

        xx = peakX*self._xMax/self._xSize
        yy = peakY*self._yMax/self._ySize
        return xx, yy

    def getClump(self):
        psfImage = self.getImage()
        #
        # Embed psfImage into a larger image so we can smooth when measuring it
        #
        width, height = psfImage.getWidth(), psfImage.getHeight()
        largeImg = psfImage.Factory(2*width, 2*height)
        largeImg.set(0)

        bbox = afwImage.BBox(afwImage.PointI(width, height), width, height)
        subLargeImg = psfImage.Factory(largeImg, bbox)
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
        mpsfImage.setXY0(afwImage.PointI(-width, -height))
        del msk
        del var
        exposure = afwImage.makeExposure(mpsfImage)
        
        #
        # Next run an object detector
        #
        max = afwMath.makeStatistics(psfImage, afwMath.MAX).getValue()
        threshold = afwDetection.Threshold(max)
        
        ds = afwDetection.DetectionSetF(mpsfImage, threshold, "DETECTED")
        objects = ds.getFootprints()
        #
        # And measure it.  This policy isn't the one we use to measure
        # Sources, it's only used to characterize this PSF histogram
        #
        psfImagePolicy = policy.Policy()
        psfImagePolicy.add("centroidAlgorithm", "NAIVE")
        psfImagePolicy.add("shapeAlgorithm", "SDSS")
        psfImagePolicy.add("photometryAlgorithm", "NAIVE")
        psfImagePolicy.add("apRadius", 3.0)
        
        sigma = 1
        psf = algorithms.createPSF("DoubleGaussian", 1, 1, sigma)
        measureSources = algorithms.makeMeasureSources(exposure,
                                                       psfImagePolicy,
                                                       psf)
        
        sourceList = afwDetection.SourceSet()

        Imax = None                     # highest peak
        for i in range(len(objects)):
            source = afwDetection.Source()
            sourceList.append(source)
            source.setId(i)
            
            try:
                measureSources.apply(source, objects[i])
            except Exception, e:
                continue
            
            x, y = source.getXAstrom(), source.getYAstrom()
            val = mpsfImage.getImage().get(int(x), int(y))

            if Imax is None or val > Imax:
                Imax = val
                psfClumpX, psfClumpY = x, y
                psfClumpIxx = source.getIxx()
                psfClumpIxy = source.getIxy()
                psfClumpIyy = source.getIyy()
        #
        # Show us the Histogram
        #
        if display:
            frame = 1
            ds9.mtv(mpsfImage.Factory(mpsfImage, afwImage.BBox(afwImage.PointI(width, height), width, height)),
                    frame=frame)
            ds9.dot("+", psfClumpX, psfClumpY, ctype=ds9.YELLOW, frame=frame)
            ds9.dot("@:%g,%g,%g" % (psfClumpIxx, psfClumpIxy, psfClumpIyy), psfClumpX, psfClumpY,
                    ctype=ds9.YELLOW, frame=frame)

            ds9.dot("PSF Image", 0, 0, frame=frame)
        #
        # Check that IxxMin/IyyMin is not too small
        #
        IzzMin = 0.5
        if psfClumpIxx < IzzMin or psfClumpIyy < IzzMin:
            psfClumpIxx, psfClumpIxy, psfClumpIyy = IzzMin, 0, IzzMin
        #
        # Convert psfClump[XY] (in psfImage's coordinates) back to Ixx/Iyy
        #
        psfClumpX, psfClumpY = self.peakToIxx(psfClumpX, psfClumpY)

        return psfClumpX, psfClumpY, psfClumpIxx, psfClumpIxy, psfClumpIyy

def getPsf(exposure, sourceList, moPolicy, sdqaRatings):
    """Return the PSF"""
    #
    # OK, we have all the source.  Let's do something with them
    #
    def goodPsfCandidate(source, fluxLim=moPolicy.get("fluxLim")):
        """Should this object be included in the Ixx v. Iyy image?""" 

        badFlags = algorithms.Flags.EDGE | \
                   algorithms.Flags.INTERP_CENTER | algorithms.Flags.SATUR_CENTER | algorithms.Flags.PEAKCENTER

        if source.getFlagForDetection() & badFlags:
            return False

        if fluxLim != None and source.getPsfFlux() < fluxLim: # ignore faint objects
            return False

        return True            
    #
    mi = exposure.getMaskedImage()
    #
    # Create an Image of Ixx v. Iyy, i.e. a 2-D histogram
    #
    psfHist = PsfShapeHistogram()

    if display:
        frame = 0
        ds9.mtv(mi, frame=frame)

    for source in sourceList:
        if goodPsfCandidate(source):
            psfHist.insert(source)
            
        if display:
            ctype = ds9.GREEN if goodPsfCandidate(source) else ds9.RED
            ds9.dot("o", source.getXAstrom() - mi.getX0(), source.getYAstrom() - mi.getY0(), frame=frame, ctype=ctype)

    psfClumpX, psfClumpY, psfClumpIxx, psfClumpIxy, psfClumpIyy = psfHist.getClump()
    #
    # Go through and find all the PSF-like objects
    #
    # We'll split the image into a number of cells, each of which contributes only
    # one PSF candidate star
    #
    sizePsfCellX = moPolicy.getInt("sizeCellX")
    sizePsfCellY = moPolicy.getInt("sizeCellY")

    psfCellSet = afwMath.SpatialCellSet(afwImage.BBox(afwImage.PointI(mi.getX0(), mi.getY0()),
                                                      mi.getWidth(), mi.getHeight()),
                                        sizePsfCellX, sizePsfCellY)

    psfStars = []

    det = psfClumpIxx*psfClumpIyy - psfClumpIxy*psfClumpIxy
    try:
        a, b, c = psfClumpIyy/det, -psfClumpIxy/det, psfClumpIxx/det
    except ZeroDivisionError:
        a, b, c = 1e4, 0, 1e4

    for source in sourceList:
        Ixx, Ixy, Iyy = source.getIxx(), source.getIxy(), source.getIyy()
        dx, dy = (Ixx - psfClumpX), (Iyy - psfClumpY)

        if math.sqrt(a*dx*dx + 2*b*dx*dy + c*dy*dy) < 2: # A test for > would be confused by NaN's
            if not goodPsfCandidate(source):
                continue

            try:
                psfCellSet.insertCandidate(algorithms.makePsfCandidate(source, mi))

                if display:
                    ds9.dot("x", source.getXAstrom() - mi.getX0(), source.getYAstrom() - mi.getY0(),
                            size=3, frame=frame, ctype=ds9.GREEN)
            except Exception, e:
                continue

            psfStars += [source]
    #
    # setWidth/setHeight are class static, but we'd need to know that the class was <float> to use that info; e.g.
    #     afwMath.SpatialCellImageCandidateF_setWidth(21)
    #
    psfCandidate = algorithms.makePsfCandidate(source, mi)
    psfCandidate.setWidth(21)
    psfCandidate.setHeight(21)
    del psfCandidate
    #
    # Do a PCA decomposition of those PSF candidates
    #
    nEigenComponents = moPolicy.getInt("nEigenComponents")
    spatialOrder  = moPolicy.getInt("spatialOrder")
    nStarPerCell = moPolicy.getInt("nStarPerCell")
    kernelSize = moPolicy.getInt("kernelSize")
    nStarPerCellSpatialFit = moPolicy.getInt("nStarPerCellSpatialFit")
    tolerance = moPolicy.getDouble("tolerance")
    reducedChi2ForPsfCandidates = moPolicy.getDouble("reducedChi2ForPsfCandidates")
    nIterForPsf = moPolicy.getInt("nIterForPsf")

    for iter in range(nIterForPsf):
        #
        # First estimate our PSF
        #
        pair = algorithms.createKernelFromPsfCandidates(psfCellSet, nEigenComponents, spatialOrder,
                                                        kernelSize, nStarPerCell)
        kernel, eigenValues = pair[0], pair[1]; del pair

        pair = algorithms.fitSpatialKernelFromPsfCandidates(kernel, psfCellSet, nStarPerCellSpatialFit, tolerance)
        status, chi2 = pair[0], pair[1]; del pair

        psf = algorithms.createPSF("PCA", kernel)
        #
        # Then clip out bad fits
        #
        psfCandidate = algorithms.makePsfCandidate(source, mi)
        nu = psfCandidate.getWidth()*psfCandidate.getHeight() - 1 # number of degrees of freedom/star for chi^2
        del psfCandidate

        for cell in psfCellSet.getCellList():
            for cand in cell.begin(False): # include bad candidates
                cand = algorithms.cast_PsfCandidateF(cand)

                rchi2 = cand.getChi2()/nu

                if rchi2 > reducedChi2ForPsfCandidates:
                    cand.setStatus(afwMath.SpatialCellCandidate.BAD)
    #
    # Display code for debugging
    #
    if display:
        #
        # Show us the ccandidates
        #
        mos = displayUtils.Mosaic()

        stamps = []; stampInfo = []
        for cell in psfCellSet.getCellList():
            for cand in cell.begin(False): # include bad candidates
                cand = algorithms.cast_PsfCandidateF(cand)

                if not cand.isBad() and display:
                    rchi2 = cand.getChi2()/nu

                    try:
                        im = cand.getImage()
                        stamps.append(im)
                        stampInfo.append("%d %.1f" % (cand.getSource().getId(), rchi2))
                    except Exception, e:
                        pass

        frame = 4
        mos.makeMosaic(stamps, frame=frame)
        mos.drawLabels(stampInfo, frame=frame)
        ds9.dot("PsfCandidates", 0, -3, frame=frame)
        #
        # We have a PSF. Possibly show it to us
        #
        eigenImages = []
        for k in afwMath.cast_LinearCombinationKernel(psf.getKernel()).getKernelList():
            im = afwImage.ImageD(k.getDimensions())
            k.computeImage(im, False)
            eigenImages.append(im)

        frame = 5
        mos.makeMosaic(eigenImages, frame=frame)
        ds9.dot("Eigen Images", 0, 0, frame=frame)

        frame = 6
        psfImages = []
        labels = []
        nx, ny = 3, 3
        for ix in range(nx):
            for iy in range(ny):
                x = (ix + 0.5)*exposure.getWidth()/nx
                y = (iy + 0.5)*exposure.getHeight()/ny

                psfImages.append(psf.getImage(x, y))
                labels.append("PSF(%d,%d)" % (int(x), int(y)))

        mos.makeMosaic(psfImages, frame=frame)
        mos.drawLabels(labels, frame=frame)
    #
    # Generate some stuff for SDQA
    #
    # Count PSF stars
    #
    numGoodStars = 0
    numAvailStars = 0

    for cell in psfCellSet.getCellList():
        numGoodStars += cell.size()

    for cell in psfCellSet.getCellList():
        for cand in cell.begin(False):  # don't ignore BAD stars
            numAvailStars += 1

    sdqaRatings.append(sdqa.SdqaRating("phot.psf.spatialFitChi2", chi2,  -1,
        sdqa.SdqaRating.AMP))
    sdqaRatings.append(sdqa.SdqaRating("phot.psf.numGoodStars", numGoodStars,
        0, sdqa.SdqaRating.AMP))
    sdqaRatings.append(sdqa.SdqaRating("phot.psf.numAvailStars",
        numAvailStars,  0, sdqa.SdqaRating.AMP))
    sdqaRatings.append(sdqa.SdqaRating("phot.psf.spatialLowOrdFlag", 0,  0,
        sdqa.SdqaRating.AMP))

    return (psf, psfCellSet)
