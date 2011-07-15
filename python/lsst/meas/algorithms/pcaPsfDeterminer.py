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
import sys
import lsst.daf.base as dafBase
import lsst.afw.geom as afwGeom
import lsst.afw.detection as afwDetection
import lsst.afw.display.ds9 as ds9
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.sdqa as sdqa
import algorithmsLib
import utils as maUtils
import numpy

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
        self._lambda                 = policy.get("lambda")
        self._reducedChi2ForPsfCandidates = policy.get("reducedChi2ForPsfCandidates")
        self._nIterForPsf            = policy.get("nIterForPsf")

    def _fitPsf(self, exposure, psfCellSet):
        # Determine KL components
        kernel, eigenValues = algorithmsLib.createKernelFromPsfCandidates(
            psfCellSet, exposure.getDimensions(), self._nEigenComponents, self._spatialOrder,
            self._kernelSize, self._nStarPerCell, self._constantWeight)

        # Express eigenValues in units of reduced chi^2 per star
        size = self._kernelSize + 2*self._borderWidth
        nu = size*size - 1                  # number of degrees of freedom/star for chi^2    
        eigenValues = [l/float(algorithmsLib.countPsfCandidates(psfCellSet, self._nStarPerCell)*nu)
                       for l in eigenValues]
        
        # Fit spatial model
        status, chi2 = algorithmsLib.fitSpatialKernelFromPsfCandidates(
            kernel, psfCellSet, self._nonLinearSpatialFit,
            self._nStarPerCellSpatialFit, self._tolerance, self._lambda)
        
        psf = afwDetection.createPsf("PCA", kernel)

        return psf, eigenValues


    def determinePsf(self, exposure, psfCandidateList, sdqaRatingSet=None):
        """Determine a PCA PSF model for an exposure given a list of PSF candidates
        
        @param[in] exposure: exposure containing the psf candidates (lsst.afw.image.Exposure)
        @param[in] psfCandidateList: a sequence of PSF candidates (each an lsst.meas.algorithms.PsfCandidate);
            typically obtained by detecting sources and then running them through a star selector
        @param[in,out] sdqaRatingSet: an lsst.sdqa.SdqaRatingSet(); if not None it will gain some SDQA ratings
    
        @return psf: an lsst.meas.algorithms.PcaPsf
        """
        import lsstDebug
        display = lsstDebug.Info(__name__).display 
        displayExposure = lsstDebug.Info(__name__).displayExposure     # display the Exposure + spatialCells 
        displayPsfCandidates = lsstDebug.Info(__name__).displayPsfCandidates # show the viable candidates 
        displayIterations = lsstDebug.Info(__name__).displayIterations # display on each PSF iteration 
        displayPsfComponents = lsstDebug.Info(__name__).displayPsfComponents # show the PCA components
        displayPsfMosaic = lsstDebug.Info(__name__).displayPsfMosaic   # show mosaic of reconstructed PSF(x,y)
        matchKernelAmplitudes = lsstDebug.Info(__name__).matchKernelAmplitudes # match Kernel amplitudes for spatial plots
        keepMatplotlibPlots = lsstDebug.Info(__name__).keepMatplotlibPlots # Keep matplotlib alive post mortem
        displayPsfSpatialModel = lsstDebug.Info(__name__).displayPsfSpatialModel # Plot spatial model?
        showBadCandidates = lsstDebug.Info(__name__).showBadCandidates # Include bad candidates 
        normalizeResiduals = lsstDebug.Info(__name__).normalizeResiduals # Normalise residuals by object amplitude 
        pause = lsstDebug.Info(__name__).pause                         # Prompt user after each iteration?
         
        if display > 1: 
            pause = True

        mi = exposure.getMaskedImage()
        
        if len(psfCandidateList) == 0:
            raise RuntimeError("No PSF candidates supplied.")

        # construct and populate a spatial cell set
        bbox = mi.getBBox(afwImage.PARENT)
        psfCellSet = afwMath.SpatialCellSet(bbox, self._sizeCellX, self._sizeCellY)
        for psfCandidate in psfCandidateList:
            psfCellSet.insertCandidate(psfCandidate)

        # Set size of image returned around candidate
        psfCandidateList[0].setHeight(self._kernelSize)
        psfCandidateList[0].setWidth(self._kernelSize)

        if display:
            frame = 0
            if displayExposure:
                ds9.mtv(exposure, frame=frame, title="psf determination")
                maUtils.showPsfSpatialCells(exposure, psfCellSet, self._nStarPerCell,
                                            symb="o", ctype=ds9.CYAN, size=4, frame=frame)
        
        #
        # Do a PCA decomposition of those PSF candidates
        #
        size = self._kernelSize + 2*self._borderWidth
        nu = size*size - 1                  # number of degrees of freedom/star for chi^2    
    
        reply = "y"                         # used in interactive mode
        for iter in range(self._nIterForPsf):
            if display and displayPsfCandidates: # Show a mosaic of usable PSF candidates
                #
                import lsst.afw.display.utils as displayUtils
    
                stamps = []
                for cell in psfCellSet.getCellList():
                    for cand in cell.begin(not showBadCandidates): # maybe include bad candidates
                        cand = algorithmsLib.cast_PsfCandidateF(cand)
                            
                        try:
                            im = cand.getImage().getImage()
    
                            chi2 = cand.getChi2()
                            if chi2 > 1e100:
                                chi2Str = ""
                            else:
                                chi2Str = " %.1f" % (chi2)
    
                            stamps.append((cand.getImage(), "%d%s" % (cand.getSource().getId(), chi2Str),
                                           cand.getStatus()))
                        except Exception, e:
                            continue
    
                mos = displayUtils.Mosaic()
                for im, label, status in stamps:
                    im = type(im)(im, True)
                    try:
                        im /= afwMath.makeStatistics(im, afwMath.MAX).getValue()
                    except NotImplementedError:
                        pass
    
                    mos.append(im, label,
                               ds9.GREEN if status == afwMath.SpatialCellCandidate.GOOD else
                               ds9.YELLOW if status == afwMath.SpatialCellCandidate.UNKNOWN else ds9.RED)
                               
    
                mos.makeMosaic(frame=7, title="ImagePca")

            #
            # First, estimate the PSF
            #            
            psf, eigenValues = self._fitPsf(exposure, psfCellSet)

            #
            # In clipping, allow all candidates to be innocent until proven guilty on this iteration
            # 
            for cell in psfCellSet.getCellList():
                for cand in cell.begin(False): # include bad candidates
                    cand = algorithmsLib.cast_PsfCandidateF(cand)
                    cand.setStatus(afwMath.SpatialCellCandidate.UNKNOWN) # until proven guilty

            #
            # Clip out bad fits based on raw chi^2
            #
            minChi2 = self._reducedChi2ForPsfCandidates*1.0*(float(self._nIterForPsf)/(iter + 1))
            if minChi2 < self._reducedChi2ForPsfCandidates:
                minChi2 = self._reducedChi2ForPsfCandidates

            for cell in psfCellSet.getCellList():
                for cand in cell.begin(False): # include bad candidates
                    cand = algorithmsLib.cast_PsfCandidateF(cand)
                    rchi2 = cand.getChi2()  # reduced chi^2 when fitting PSF to candidate
                    if rchi2 < 0 or rchi2 > minChi2:
                        #print "chi^2 clipping %d (%f,%f): %f" % (cand.getSource().getId(), cand.getXCenter(), cand.getYCenter(), rchi2)
                        cand.setStatus(afwMath.SpatialCellCandidate.BAD)
                        if rchi2 < 0:
                            print "RHL chi^2:", rchi2, cand.getChi2(), nu

            #
            # Clip out bad fits based on spatial fitting.
            #
            # This appears to be better at getting rid of sources that have a single dominant kernel component
            # (other than the zeroth; e.g., a nearby contaminant) because the surrounding sources (which help
            # set the spatial model) don't contain that kernel component, and so the spatial modeling
            # downweights the component.
            #

            residuals = list()
            candidates = list()
            kernel = psf.getKernel()
            noSpatialKernel = afwMath.cast_LinearCombinationKernel(psf.getKernel())
            for cell in psfCellSet.getCellList():
                for cand in cell.begin(True):
                    cand = algorithmsLib.cast_PsfCandidateF(cand)
                    candCenter = afwGeom.PointD(cand.getXCenter(), cand.getYCenter())
                    try:
                        im = cand.getImage()
                    except Exception, e:
                        continue

                    # Do fit based on entire postage stamp, by setting everything to 'detected'
                    image = im.Factory(im, True)
                    mask = image.getMask()
                    detected = mask.getPlaneBitMask("DETECTED");
                    mask |= detected

                    fit = algorithmsLib.fitKernelParamsToImage(noSpatialKernel, image, candCenter)
                    params = fit[0]
                    kernels = fit[1]
                    amp = 0.0
                    for p, k in zip(params, kernels):
                        amp += p * afwMath.cast_FixedKernel(k).getSum()

                    #print cand.getSource().getId(), [a / amp for a in params]

                    predict = [kernel.getSpatialFunction(k)(candCenter.getX(), candCenter.getY()) for
                               k in range(kernel.getNKernelParameters())]

                    residuals.append([a / amp - p for a, p in zip(params, predict)])
                    candidates.append(cand)
            
            residuals = numpy.array(residuals)
            for k in range(kernel.getNKernelParameters()):
                rms = residuals[:,k].std()
                #print "RMS for component %d is %f" % (k, rms)
                for i, cand in enumerate(candidates):
                    if numpy.fabs(residuals[i,k]) > 3.0 * rms:
                        #print "Spatial clipping %d (%f,%f) based on %d: %f" % (cand.getSource().getId(), cand.getXCenter(), cand.getYCenter(), k, residuals[i,k])
                        cand.setStatus(afwMath.SpatialCellCandidate.BAD)



            #
            # Display results
            #
            if display and displayIterations:
                if displayExposure:
                    if iter > 0:
                        ds9.erase(frame=frame)
                    maUtils.showPsfSpatialCells(exposure, psfCellSet, self._nStarPerCell, showChi2=True,
                                                symb="o", ctype=ds9.YELLOW, ctypeBad=ds9.RED, size=8, frame=frame)
                    if self._nStarPerCellSpatialFit != self._nStarPerCell:
                        maUtils.showPsfSpatialCells(exposure, psfCellSet, self._nStarPerCellSpatialFit,
                                                    symb="o", ctype=ds9.YELLOW, ctypeBad=ds9.RED,
                                                    size=10, frame=frame)
                while True:
                    try:
                        maUtils.showPsfCandidates(exposure, psfCellSet, psf=psf, frame=4,
                                                  normalize=normalizeResiduals,
                                                  showBadCandidates=showBadCandidates)
                    except:
                        if not showBadCandidates:
                            showBadCandidates = True
                            continue
                    break
    
                if displayPsfComponents:
                    maUtils.showPsf(psf, eigenValues, frame=5)
                if displayPsfMosaic:
                    maUtils.showPsfMosaic(exposure, psf, frame=6)
                if displayPsfSpatialModel:
                    maUtils.plotPsfSpatialModel(exposure, psf, psfCellSet, showBadCandidates=True,
                                                matchKernelAmplitudes=matchKernelAmplitudes,
                                                keepPlots=keepMatplotlibPlots)
    
                if pause:
                    while True:
                        try:
                            reply = raw_input("Next iteration? [ynchpqs] ").strip()
                        except EOFError:
                            reply = "n"
    
                        reply = reply.split()
                        if reply:
                            reply, args = reply[0], reply[1:]
                        else:
                            reply = ""
                            
                        if reply in ("", "c", "h", "n", "p", "q", "s", "y"):
                            if reply == "c":
                                pause = False
                            elif reply == "h":
                                print "c[ontinue without prompting] h[elp] n[o] p[db] q[uit displaying] s[ave fileName] y[es]"
                                continue
                            elif reply == "p":
                                import pdb; pdb.set_trace() 
                            elif reply == "q":
                                display = False
                            elif reply == "s":
                                fileName = args.pop(0)
                                if not fileName:
                                    print "Please provide a filename"
                                    continue
                                
                                print "Saving to %s" % fileName
                                maUtils.saveSpatialCellSet(psfCellSet, fileName=fileName)
                                continue
                            break
                        else:
                            print >> sys.stderr, "Unrecognised response: %s" % reply
    
                    if reply == "n":
                        break

        # One last time, to take advantage of the last iteration
        psf, eigenValues = self._fitPsf(exposure, psfCellSet)

        ##################
        # quick and dirty match to return a sourceSet of objects in the cellSet
        # should be faster than N^2, but not an issue for lists this size
                    
        # put sources in a dict with x,y lookup
        # must disable - Source constructor can't copy all internals and it breaks the pipe
        if False:
            sourceLookup = {}
            for s in sourceList:
                x, y = int(s.getXAstrom()), int(s.getYAstrom())
                key = str(x)+"."+str(y)
                sourceLookup[key] = s
    
            # keep only the good ones
            psfSourceSet = afwDetection.SourceSet()
            for cell in psfCellSet.getCellList():
                for cand in cell.begin(True):  # ignore bad candidates
                    x, y = int(cand.getXCenter()), int(cand.getYCenter())
                    key = str(x)+"."+str(y)
                    psfSourceSet.append(sourceLookup[key])
                    
        #
        # Display code for debugging
        #
        if display and reply != "n":
            if displayExposure:
                maUtils.showPsfSpatialCells(exposure, psfCellSet, self._nStarPerCell, showChi2=True,
                                            symb="o", ctype=ds9.YELLOW, ctypeBad=ds9.RED, size=8, frame=frame)
                if self._nStarPerCellSpatialFit != self._nStarPerCell:
                    maUtils.showPsfSpatialCells(exposure, psfCellSet, self._nStarPerCellSpatialFit,
                                                symb="o", ctype=ds9.YELLOW, ctypeBad=ds9.RED,
                                                size=10, frame=frame)
            maUtils.showPsfCandidates(exposure, psfCellSet, psf=psf, frame=4, normalize=normalizeResiduals,
                                      showBadCandidates=showBadCandidates)
                                      
            maUtils.showPsf(psf, eigenValues, frame=5)
            if displayPsfMosaic:
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
