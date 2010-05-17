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
import lsst.meas.algorithms.measureSourceUtils as maUtils
import lsst.sdqa as sdqa

import lsst.afw.display.ds9 as ds9

try:
    type(display)
except NameError:
    display = False
    displayPca = True                   # show the PCA components
    displayIterations = True            # display on each PSF iteration
    
class PsfShapeHistogram(object):
    """A class to represent a histogram of (Ixx, Iyy)"""

    def __init__(self, xSize=40, ySize=40, xMax=30, yMax=30):
        """[xy]Size is the size of the psfImage; [xy]Max are the maximum values for I[xy][xy]"""

        self._xSize, self._ySize = xSize, ySize 
        self._xMax, self._yMax = xMax, yMax
        self._psfImage = afwImage.ImageF(self._xSize, self._ySize)
        self._psfImage.set(0)

    def getImage(self):
        return self._psfImage

    def insert(self, source):
        """Insert source into the histogram."""
        try:
            i = int(source.getIxx()*self._xSize/self._xMax + 0.5)
            j = int(source.getIyy()*self._ySize/self._yMax + 0.5)
        except:
            return

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
        
        ds = afwDetection.FootprintSetF(mpsfImage, threshold, "DETECTED")
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
        e = None                        # thrown exception
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
                    title="PSF Image", frame=frame)
            if Imax is not None:
                ds9.dot("+", psfClumpX, psfClumpY, ctype=ds9.YELLOW, frame=frame)
                ds9.dot("@:%g,%g,%g" % (psfClumpIxx, psfClumpIxy, psfClumpIyy), psfClumpX, psfClumpY,
                        ctype=ds9.YELLOW, frame=frame)
        #
        if Imax is None:
            msg = "Failed to determine center of PSF clump"
            if e:
                msg += ": %s" % e

            raise RuntimeError, msg
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

def getPsf(exposure, sourceList, psfPolicy, sdqaRatings):
    """Return the PSF for the given Exposure and set of Sources, given a Policy

The policy is documented in ip/pipeline/policy/CrRejectDictionary.paf    
    """
    #
    # OK, we have all the source.  Let's do something with them
    #
    def goodPsfCandidate(source, fluxLim=psfPolicy.get("fluxLim")):
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
        ds9.mtv(mi, frame=frame, title="PSF candidates")

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
    sizePsfCellX = psfPolicy.get("sizeCellX")
    sizePsfCellY = psfPolicy.get("sizeCellY")

    psfCellSet = afwMath.SpatialCellSet(afwImage.BBox(afwImage.PointI(mi.getX0(), mi.getY0()),
                                                      mi.getWidth(), mi.getHeight()),
                                        sizePsfCellX, sizePsfCellY)

    psfStars = []

    det = psfClumpIxx*psfClumpIyy - psfClumpIxy*psfClumpIxy
    try:
        a, b, c = psfClumpIyy/det, -psfClumpIxy/det, psfClumpIxx/det
    except ZeroDivisionError:
        a, b, c = 1e4, 0, 1e4

    # psf candidate shapes must lie within this many RMS of the average shape
    # N.b. if Ixx == Iyy, Ixy = 0 the criterion is dx^2 + dy^2 < clumpNSigma*(Ixx + Iyy) == 2*clumpNSigma*Ixx
    clumpNSigma = psfPolicy.get("clumpNSigma")
    for source in sourceList:
        Ixx, Ixy, Iyy = source.getIxx(), source.getIxy(), source.getIyy()
        dx, dy = (Ixx - psfClumpX), (Iyy - psfClumpY)

        if math.sqrt(a*dx*dx + 2*b*dx*dy + c*dy*dy) < 2*clumpNSigma: # A test for > would be confused by NaN
            if not goodPsfCandidate(source):
                continue

            try:
                psfCellSet.insertCandidate(algorithms.makePsfCandidate(source, mi))

                if display:
                    ds9.dot("o", source.getXAstrom() - mi.getX0(), source.getYAstrom() - mi.getY0(),
                            size=4, frame=frame, ctype=ds9.CYAN)
            except Exception, e:
                continue

            psfStars += [source]
    #
    # Do a PCA decomposition of those PSF candidates
    #
    nonLinearSpatialFit = psfPolicy.get("nonLinearSpatialFit")
    nEigenComponents = psfPolicy.get("nEigenComponents")
    spatialOrder  = psfPolicy.get("spatialOrder")
    nStarPerCell = psfPolicy.get("nStarPerCell")
    kernelSize = psfPolicy.get("kernelSize")
    borderWidth = psfPolicy.get("borderWidth")
    nStarPerCellSpatialFit = psfPolicy.get("nStarPerCellSpatialFit")
    constantWeight = psfPolicy.get("constantWeight")
    tolerance = psfPolicy.get("tolerance")
    reducedChi2ForPsfCandidates = psfPolicy.get("reducedChi2ForPsfCandidates")
    nIterForPsf = psfPolicy.get("nIterForPsf")
    #
    # setWidth/setHeight are class static, but we'd need to know that the class was <float> to use that info;
    # e.g.
    #     afwMath.SpatialCellImageCandidateF_setWidth(21)
    #
    psfCandidate = algorithms.makePsfCandidate(source, mi)
    size = kernelSize + 2*borderWidth
    psfCandidate.setBorderWidth(borderWidth)
    psfCandidate.setWidth(size)
    psfCandidate.setHeight(size)

    nu = size*size - 1                  # number of degrees of freedom/star for chi^2    

    del psfCandidate

    for iter in range(nIterForPsf):
        if display and displayPca:      # Build a ImagePca so we can look at its Images (for debugging)
            #
            import lsst.afw.display.utils as displayUtils

            pca = afwImage.ImagePcaF()
            ids = []
            for cell in psfCellSet.getCellList():
                for cand in cell.begin(False): # include bad candidates
                    cand = algorithms.cast_PsfCandidateF(cand)
                    try:
                        im = cand.getImage().getImage()

                        max = afwMath.makeStatistics(im, afwMath.MAX).getValue()
                        if math.isnan(max) or math.isinf(max):
                            continue

                        pca.addImage(im, afwMath.makeStatistics(im, afwMath.SUM).getValue())
                        ids.append(("%d %.1f" % (cand.getSource().getId(), cand.getChi2()/361.0),
                                    ds9.GREEN if cand.getStatus() == afwMath.SpatialCellCandidate.GOOD else
                                    ds9.YELLOW if cand.getStatus() == afwMath.SpatialCellCandidate.UNKNOWN else
                                    ds9.RED))
                    except Exception, e:
                        continue

            mos = displayUtils.Mosaic(); i = 0
            for im in pca.getImageList():
                im = type(im)(im, True)
                try:
                    im /= afwMath.makeStatistics(im, afwMath.MAX).getValue()
                except NotImplementedError:
                    pass
                mos.append(im, ids[i][0], ids[i][1]); i += 1

            mos.makeMosaic(frame=7, title="ImagePca")
            del pca
        #
        # First estimate our PSF
        #
        pair = algorithms.createKernelFromPsfCandidates(psfCellSet, nEigenComponents, spatialOrder,
                                                        kernelSize, nStarPerCell, constantWeight)
        kernel, eigenValues = pair[0], pair[1]; del pair
        #
        # Express eigenValues in units of reduced chi^2 per star
        #
        eigenValues = [l/float(algorithms.countPsfCandidates(psfCellSet, nStarPerCell)*nu)
                       for l in eigenValues]

        #
        # Set the initial amplitudes if we're doing linear fits
        #
        if iter == 0 and not nonLinearSpatialFit:
            for cell in psfCellSet.getCellList():
                for cand in cell.begin(False): # include bad candidates
                    cand = algorithms.cast_PsfCandidateF(cand)
                    try:
                        if cand.getSource().getId() == 115:
                            #import pdb; pdb.set_trace() 
                            pass

                        cand.setAmplitude(afwMath.makeStatistics(cand.getImage().getImage(),
                                                                 afwMath.SUM).getValue())

                        if False and math.isnan(cand.getAmplitude()):
                            ds9.mtv(cand.getImage().getImage(), title="candidate", frame=1)
                            print "amp",  cand.getAmplitude()
                            #import pdb; pdb.set_trace() 
                    except Exception, e:
                        print "RHL", e

        pair = algorithms.fitSpatialKernelFromPsfCandidates(kernel, psfCellSet, nonLinearSpatialFit,
                                                            nStarPerCellSpatialFit, tolerance)
        status, chi2 = pair[0], pair[1]; del pair

        psf = algorithms.createPSF("PCA", kernel)
        #
        # Then clip out bad fits
        #
        for cell in psfCellSet.getCellList():
            for cand in cell.begin(False): # include bad candidates
                cand = algorithms.cast_PsfCandidateF(cand)
                cand.setStatus(afwMath.SpatialCellCandidate.UNKNOWN) # until proven guilty

                rchi2 = cand.getChi2()/nu

                if rchi2 < 0 or rchi2 > reducedChi2ForPsfCandidates*(float(nIterForPsf)/(iter + 1)):
                    cand.setStatus(afwMath.SpatialCellCandidate.BAD)
                    if rchi2 < 0:
                        print "RHL chi^2:", rchi2, cand.getChi2(), nu
                    
        if display and displayIterations:
            if iter > 0:
                ds9.erase(frame=frame)
            maUtils.showPsfSpatialCells(exposure, psfCellSet, nStarPerCell, showChi2=True,
                                        symb="o", ctype=ds9.YELLOW, size=8, frame=frame)
            if nStarPerCellSpatialFit != nStarPerCell:
                maUtils.showPsfSpatialCells(exposure, psfCellSet, nStarPerCellSpatialFit,
                                            symb="o", ctype=ds9.YELLOW, size=10, frame=frame)
            maUtils.showPsfCandidates(exposure, psfCellSet, psf=psf, frame=4)
            maUtils.showPsf(psf, eigenValues, frame=5)
            maUtils.showPsfMosaic(exposure, psf, frame=6)

            if False:
                print "Iter =", iter
                import pdb; pdb.set_trace() 
    #
    # Display code for debugging
    #
    if display:
        maUtils.showPsfSpatialCells(exposure, psfCellSet, nStarPerCell, showChi2=True,
                                    symb="o", ctype=ds9.YELLOW, size=8, frame=frame)
        if nStarPerCellSpatialFit != nStarPerCell:
            maUtils.showPsfSpatialCells(exposure, psfCellSet, nStarPerCellSpatialFit,
                                        symb="o", ctype=ds9.YELLOW, size=10, frame=frame)
        maUtils.showPsfCandidates(exposure, psfCellSet, psf=psf, frame=4)
        maUtils.showPsf(psf, eigenValues, frame=5)
        maUtils.showPsfMosaic(exposure, psf, frame=6)
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
