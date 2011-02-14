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
import lsst.daf.base as dafBase
import lsst.afw.detection as afwDetection
import lsst.afw.display.ds9 as ds9
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.sdqa as sdqa
import algorithmsLib
import utils as maUtils

class PcaPsfDeterminer(object):
    def __init__(self, policy):
        """Construct a PCA PSF Fitter

        @param[in] policy: see policy/PcaPsfDeterminerDictionary.paf
        """
        self._nonLinearSpatialFit    = policy.get("nonLinearSpatialFit")
        self._nEigenComponents       = policy.get("nEigenComponents")
        self._spatialOrder           = policy.get("spatialOrder")
        self._sizeCellX              = policy.get("sizeCellX")
        self._sizeCellY              = policy.get("sizeCellY")
        self._nStarPerCell           = policy.get("nStarPerCell")
        self._kernelSize             = policy.get("kernelSize")
        self._borderWidth            = policy.get("borderWidth")
        self._nStarPerCellSpatialFit = policy.get("nStarPerCellSpatialFit")
        self._constantWeight         = policy.get("constantWeight")
        self._tolerance              = policy.get("tolerance")
        self._reducedChi2ForPsfCandidates = policy.get("reducedChi2ForPsfCandidates")
        self._nIterForPsf            = policy.get("nIterForPsf")

    def determinePsf(self, exposure, psfCandidateList, sdqaRatingSet=None):
        """Determine a PCA PSF model for an exposure given a list of PSF candidates
        
        @param[in] exposure: exposure containing the psf candidates (lsst.afw.image.Exposure)
        @param[in] psfCandidateList: a sequence of PSF candidates (each an lsst.meas.algorithms.PsfCandidate);
            typically obtained by detecting sources and then running them through a star selector
        @param[in,out] sdqaRatingSet: an lsst.sdqa.SdqaRatingSet(); if not None it will gain some SDQA ratings
    
        @return psf: an lsst.meas.algorithms.PcaPsf
        """
        try:
            import lsstDebug
    
            display = lsstDebug.Info(__name__).display
            displayPca = lsstDebug.Info(__name__).displayPca               # show the PCA components
            displayIterations = lsstDebug.Info(__name__).displayIterations # display on each PSF iteration
        except ImportError, e:
            try:
                type(display)
            except NameError:
                display = False
                displayPca = True                   # show the PCA components
                displayIterations = True            # display on each PSF iteration

        mi = exposure.getMaskedImage()
        
        # construct and populate a spatial cell set
        bbox = afwImage.BBox(afwImage.PointI(mi.getX0(), mi.getY0()), mi.getWidth(), mi.getHeight())
        psfCellSet = afwMath.SpatialCellSet(bbox, self._sizeCellX, self._sizeCellY)
        for psfCandidate in psfCandidateList:
            psfCellSet.insertCandidate(psfCandidate)

        if display:
            frame = 0
        
        #
        # Do a PCA decomposition of those PSF candidates
        #
        size = self._kernelSize + 2*self._borderWidth
        nu = size*size - 1                  # number of degrees of freedom/star for chi^2    
    
        for iter in range(self._nIterForPsf):
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
            pair = algorithmsLib.createKernelFromPsfCandidates(psfCellSet, self._nEigenComponents, self._spatialOrder,
                                                            self._kernelSize, self._nStarPerCell, self._constantWeight)
            kernel, eigenValues = pair[0], pair[1]; del pair
            #
            # Express eigenValues in units of reduced chi^2 per star
            #
            eigenValues = [l/float(algorithmsLib.countPsfCandidates(psfCellSet, self._nStarPerCell)*nu)
                           for l in eigenValues]
    
            #
            # Set the initial amplitudes if we're doing linear fits
            #
            if iter == 0 and not self._nonLinearSpatialFit:
                for cell in psfCellSet.getCellList():
                    for cand in cell.begin(False): # include bad candidates
                        cand = algorithmsLib.cast_PsfCandidateF(cand)
                        try:
                            cand.setAmplitude(afwMath.makeStatistics(cand.getImage().getImage(),
                                                                     afwMath.SUM).getValue())
                        except Exception, e:
                            print "RHL", e
    
            pair = algorithmsLib.fitSpatialKernelFromPsfCandidates(kernel, psfCellSet, self._nonLinearSpatialFit,
                                                                self._nStarPerCellSpatialFit, self._tolerance)
            status, chi2 = pair[0], pair[1]; del pair
    
            psf = afwDetection.createPsf("PCA", kernel)
            #
            # Then clip out bad fits
            #
            for cell in psfCellSet.getCellList():
                for cand in cell.begin(False): # include bad candidates
                    cand = algorithmsLib.cast_PsfCandidateF(cand)
                    cand.setStatus(afwMath.SpatialCellCandidate.UNKNOWN) # until proven guilty
    
                    rchi2 = cand.getChi2()/nu
    
                    if rchi2 < 0 or rchi2 > self._reducedChi2ForPsfCandidates*(float(self._nIterForPsf)/(iter + 1)):
                        cand.setStatus(afwMath.SpatialCellCandidate.BAD)
                        if rchi2 < 0:
                            print "RHL chi^2:", rchi2, cand.getChi2(), nu
                        
            if display and displayIterations:
                if iter > 0:
                    ds9.erase(frame=frame)
                maUtils.showPsfSpatialCells(exposure, psfCellSet, self._nStarPerCell, showChi2=True,
                                            symb="o", ctype=ds9.YELLOW, size=8, frame=frame)
                if self._nStarPerCellSpatialFit != self._nStarPerCell:
                    maUtils.showPsfSpatialCells(exposure, psfCellSet, self._nStarPerCellSpatialFit,
                                                symb="o", ctype=ds9.YELLOW, size=10, frame=frame)
                maUtils.showPsfCandidates(exposure, psfCellSet, psf=psf, frame=4, normalize=False)
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
    
        if sdqaRatingSet != None:
            sdqaRatingSet.append(sdqa.SdqaRating("phot.psf.spatialFitChi2", chi2,  -1,
                sdqa.SdqaRating.CCD))
            sdqaRatingSet.append(sdqa.SdqaRating("phot.psf.numGoodStars", numGoodStars,
                0, sdqa.SdqaRating.CCD))
            sdqaRatingSet.append(sdqa.SdqaRating("phot.psf.numAvailStars",
                numAvailStars,  0, sdqa.SdqaRating.CCD))
            sdqaRatingSet.append(sdqa.SdqaRating("phot.psf.spatialLowOrdFlag", 0,  0,
                sdqa.SdqaRating.CCD))
    
        return psf, psfCellSet
