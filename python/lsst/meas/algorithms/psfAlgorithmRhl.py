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

# This is not a minimal set of imports
import glob, math, os, sys
from math import *
import numpy
import eups
import lsst.daf.base as dafBase
import lsst.pex.logging as logging
import lsst.pex.policy as policy
import lsst.afw.detection as afwDetection
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.meas.algorithms as algorithms
import lsst.meas.algorithms.defects as defects
import lsst.meas.algorithms.utils as maUtils
import lsst.sdqa as sdqa

import lsst.afw.display.ds9 as ds9
    
def getPsf(exposure, sourceList, psfCellSet, psfAlgPolicy, sdqaRatings):
    """Return the PSF for the given Exposure and set of Sources, given a Policy

The policy is documented in ip/pipeline/policy/CrRejectDictionary.paf    
    """
    import lsstDebug

    display = lsstDebug.Info(__name__).display
    displayExposure = lsstDebug.Info(__name__).displayExposure     # display the Exposure + spatialCells
    displayPsfCandidates = lsstDebug.Info(__name__).displayPsfCandidates # show the viable candidates
    displayIterations = lsstDebug.Info(__name__).displayIterations # display on each PSF iteration
    displayPsfMosaic = lsstDebug.Info(__name__).displayPsfMosaic   # show mosaic of reconstructed PSF(x,y)
    showBadCandidates = lsstDebug.Info(__name__).showBadCandidates # Include bad candidates
    normalizeResiduals = lsstDebug.Info(__name__).normalizeResiduals # Normalise residuals by object amplitude
    pause = lsstDebug.Info(__name__).pause                         # Prompt user after each iteration?
    
    if display > 1:
        pause = True
            
    mi = exposure.getMaskedImage()
    
    #
    # Unpack policy
    #
    nonLinearSpatialFit    = psfAlgPolicy.get("nonLinearSpatialFit")
    nEigenComponents       = psfAlgPolicy.get("nEigenComponents")
    spatialOrder           = psfAlgPolicy.get("spatialOrder")
    nStarPerCell           = psfAlgPolicy.get("nStarPerCell")
    kernelSize             = psfAlgPolicy.get("kernelSize")
    borderWidth            = psfAlgPolicy.get("borderWidth")
    nStarPerCellSpatialFit = psfAlgPolicy.get("nStarPerCellSpatialFit")
    constantWeight         = psfAlgPolicy.get("constantWeight")
    tolerance              = psfAlgPolicy.get("tolerance")
    lambda_                = psfAlgPolicy.get("lambda") # lambda's a reserved python keyword
    reducedChi2ForPsfCandidates = psfAlgPolicy.get("reducedChi2ForPsfCandidates")
    nIterForPsf            = psfAlgPolicy.get("nIterForPsf")

    if display:
        frame = 0
        if displayExposure:
            ds9.mtv(exposure, frame=frame, title="psf determination")
            maUtils.showPsfSpatialCells(exposure, psfCellSet, nStarPerCell, showMoments=False,
                                        symb="o", ctype=ds9.CYAN, size=4, frame=frame)
    #
    # Do a PCA decomposition of those PSF candidates
    #
    size = kernelSize + 2*borderWidth
    nu = size*size - 1                  # number of degrees of freedom/star for chi^2    

    reply = "y"                         # used in interactive mode
    for iter in range(nIterForPsf):
        if display and displayPsfCandidates: # Show a mosaic of usable PSF candidates
            #
            import lsst.afw.display.utils as displayUtils

            stamps = []
            for cell in psfCellSet.getCellList():
                for cand in cell.begin(not showBadCandidates): # maybe include bad candidates
                    cand = algorithms.cast_PsfCandidateF(cand)
                        
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
                        cand.setAmplitude(afwMath.makeStatistics(cand.getImage().getImage(),
                                                                 afwMath.SUM).getValue())
                    except Exception, e:
                        print "RHL", e

        pair = algorithms.fitSpatialKernelFromPsfCandidates(kernel, psfCellSet, nonLinearSpatialFit,
                                                            nStarPerCellSpatialFit, tolerance, lambda_)
        status, chi2 = pair[0], pair[1]; del pair

        psf = afwDetection.createPsf("PCA", kernel)
        #
        # Then clip out bad fits
        #
        minChi2 = reducedChi2ForPsfCandidates*1.0*(float(nIterForPsf)/(iter + 1))
        if minChi2 < reducedChi2ForPsfCandidates:
            minChi2 = reducedChi2ForPsfCandidates

        for cell in psfCellSet.getCellList():
            for cand in cell.begin(False): # include bad candidates
                cand = algorithms.cast_PsfCandidateF(cand)
                cand.setStatus(afwMath.SpatialCellCandidate.UNKNOWN) # until proven guilty

                rchi2 = cand.getChi2()  # reduced chi^2 when fitting PSF to candidate

                if rchi2 < 0 or rchi2 > minChi2:
                    cand.setStatus(afwMath.SpatialCellCandidate.BAD)
                    if rchi2 < 0:
                        print "RHL chi^2:", rchi2, cand.getChi2(), nu
                    
        if display and displayIterations:
            if displayExposure:
                if iter > 0:
                    ds9.erase(frame=frame)
                maUtils.showPsfSpatialCells(exposure, psfCellSet, nStarPerCell, showChi2=True,
                                            symb="o", ctype=ds9.YELLOW, size=8, frame=frame)
                if nStarPerCellSpatialFit != nStarPerCell:
                    maUtils.showPsfSpatialCells(exposure, psfCellSet, nStarPerCellSpatialFit,
                                                symb="o", ctype=ds9.YELLOW, size=10, frame=frame)
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

            maUtils.showPsf(psf, eigenValues, frame=5)
            if displayPsfMosaic:
                maUtils.showPsfMosaic(exposure, psf, frame=6)

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
            maUtils.showPsfSpatialCells(exposure, psfCellSet, nStarPerCell, showChi2=True,
                                        symb="o", ctype=ds9.YELLOW, size=8, frame=frame)
            if nStarPerCellSpatialFit != nStarPerCell:
                maUtils.showPsfSpatialCells(exposure, psfCellSet, nStarPerCellSpatialFit,
                                            symb="o", ctype=ds9.YELLOW, size=10, frame=frame)
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

    sdqaRatings.append(sdqa.SdqaRating("phot.psf.spatialFitChi2", chi2,  -1,
        sdqa.SdqaRating.CCD))
    sdqaRatings.append(sdqa.SdqaRating("phot.psf.numGoodStars", numGoodStars,
        0, sdqa.SdqaRating.CCD))
    sdqaRatings.append(sdqa.SdqaRating("phot.psf.numAvailStars",
        numAvailStars,  0, sdqa.SdqaRating.CCD))
    sdqaRatings.append(sdqa.SdqaRating("phot.psf.spatialLowOrdFlag", 0,  0,
        sdqa.SdqaRating.CCD))

    return (psf, psfCellSet, sourceList)
