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

import numpy

import lsst.daf.base as dafBase
import lsst.pex.config as pexConfig
import lsst.pex.exceptions as pexExceptions
import lsst.pex.logging as pexLog
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses as afwEll
import lsst.afw.detection as afwDetection
import lsst.afw.display.ds9 as ds9
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.afw.math as afwMath
from . import algorithmsLib
from . import utils as maUtils

class PcaPsfDeterminerConfig(pexConfig.Config):
    nonLinearSpatialFit = pexConfig.Field(
        doc = "Use non-linear fitter for spatial variation of Kernel",
        dtype = bool,
        default = False,
    )
    nEigenComponents = pexConfig.Field(
        doc = "number of eigen components for PSF kernel creation",
        dtype = int,
        default = 4,
    )
    spatialOrder = pexConfig.Field(
        doc = "specify spatial order for PSF kernel creation",
        dtype = int,
        default = 2,
    )
    sizeCellX = pexConfig.Field(
        doc = "size of cell used to determine PSF (pixels, column direction)",
        dtype = int,
        default = 256,
#        minValue = 10,
        check = lambda x: x >= 10,
    )
    sizeCellY = pexConfig.Field(
        doc = "size of cell used to determine PSF (pixels, row direction)",
        dtype = int,
        default = sizeCellX.default,
#        minValue = 10,
        check = lambda x: x >= 10,
    )
    nStarPerCell = pexConfig.Field(
        doc = "number of stars per psf cell for PSF kernel creation",
        dtype = int,
        default = 3,
    )
    kernelSize = pexConfig.Field(
        doc = "radius of the kernel to create, relative to the square root of the stellar quadrupole moments",
        dtype = int,
        default = 7,
    )
    kernelSizeMin = pexConfig.Field(
        doc = "Minimum radius of the kernel",
        dtype = int,
        default = 25,
    )
    kernelSizeMax = pexConfig.Field(
        doc = "Maximum radius of the kernel",
        dtype = int,
        default = 45,
    )
    borderWidth = pexConfig.Field(
        doc = "Number of pixels to ignore around the edge of PSF candidate postage stamps",
        dtype = int,
        default = 0,
    )
    nStarPerCellSpatialFit = pexConfig.Field(
        doc = "number of stars per psf Cell for spatial fitting",
        dtype = int,
        default = 5,
    )
    constantWeight = pexConfig.Field(
        doc = "Should each PSF candidate be given the same weight, independent of magnitude?",
        dtype = bool,
        default = True,
    )
    nIterForPsf = pexConfig.Field(
        doc = "number of iterations of PSF candidate star list",
        dtype = int,
        default = 3,
    )
    tolerance = pexConfig.Field(
        doc = "tolerance of spatial fitting",
        dtype = float,
        default = 1e-2,
    )
    lam = pexConfig.Field(
        doc = "floor for variance is lam*data",
        dtype = float,
        default = 0.05,
    )
    reducedChi2ForPsfCandidates = pexConfig.Field(
        doc = "for psf candidate evaluation",
        dtype = float,
        default = 2.0,
    )
    spatialReject = pexConfig.Field(
        doc = "Rejection threshold (stdev) for candidates based on spatial fit",
        dtype = float,
        default = 3.0,
    )

class PcaPsfDeterminer(object):
    ConfigClass = PcaPsfDeterminerConfig

    def __init__(self, config):
        """Construct a PCA PSF Fitter

        @param[in] config: instance of PcaPsfDeterminerConfig
        """
        self.config = config
        # N.b. name of component is meas.algorithms.psfDeterminer so you can turn on psf debugging
        # independent of which determiner is active
        self.debugLog = pexLog.Debug("meas.algorithms.psfDeterminer")
        self.warnLog = pexLog.Log(pexLog.getDefaultLog(), "meas.algorithms.psfDeterminer")

    def _fitPsf(self, exposure, psfCellSet, kernelSize, nEigenComponents):
        #
        # Loop trying to use nEigenComponents, but allowing smaller numbers if necessary
        #
        for nEigen in range(nEigenComponents, 0, -1):
            # Determine KL components
            try:
                kernel, eigenValues = algorithmsLib.createKernelFromPsfCandidates(
                    psfCellSet, exposure.getDimensions(), exposure.getXY0(), nEigen,
                    self.config.spatialOrder, kernelSize, self.config.nStarPerCell,
                    bool(self.config.constantWeight))

                break                   # OK, we can get nEigen components
            except pexExceptions.LsstCppException, e:
                if not isinstance(e.message, pexExceptions.LengthErrorException):
                    raise

                if nEigen == 1:         # can't go any lower
                    raise
                    
                msg = e.message.what().strip().split("\n")[-1] # message from exception
                msg = msg.split(":")[-1].strip()               # remove "0: Message: " prefix

                self.warnLog.log(pexLog.Log.WARN, "%s: reducing number of eigen components" % msg)
        #
        # We got our eigen decomposition so let's use it
        #
        # Express eigenValues in units of reduced chi^2 per star
        size = kernelSize + 2*self.config.borderWidth
        nu = size*size - 1                  # number of degrees of freedom/star for chi^2    
        eigenValues = [l/float(algorithmsLib.countPsfCandidates(psfCellSet, self.config.nStarPerCell)*nu)
                       for l in eigenValues]

        # Fit spatial model
        status, chi2 = algorithmsLib.fitSpatialKernelFromPsfCandidates(
            kernel, psfCellSet, bool(self.config.nonLinearSpatialFit),
            self.config.nStarPerCellSpatialFit, self.config.tolerance, self.config.lam)

        psf = algorithmsLib.PcaPsf(kernel)

        return psf, eigenValues, nEigen, chi2


    def determinePsf(self, exposure, psfCandidateList, metadata=None, flagKey=None):
        """Determine a PCA PSF model for an exposure given a list of PSF candidates
        
        @param[in] exposure: exposure containing the psf candidates (lsst.afw.image.Exposure)
        @param[in] psfCandidateList: a sequence of PSF candidates (each an lsst.meas.algorithms.PsfCandidate);
            typically obtained by detecting sources and then running them through a star selector
        @param[in,out] metadata  a home for interesting tidbits of information
        @param[in] flagKey: schema key used to mark sources actually used in PSF determination
    
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
        psfCellSet = afwMath.SpatialCellSet(bbox, self.config.sizeCellX, self.config.sizeCellY)
        sizes = numpy.ndarray(len(psfCandidateList))
        for i, psfCandidate in enumerate(psfCandidateList):
            try:
                psfCellSet.insertCandidate(psfCandidate)
            except Exception, e:
                self.debugLog.debug(2, "Skipping PSF candidate %d of %d: %s" % (i, len(psfCandidateList), e))
                continue
            source = psfCandidate.getSource()

            quad = afwEll.Quadrupole(source.getIxx(), source.getIyy(), source.getIxy())
            axes = afwEll.Axes(quad)
            sizes[i] = axes.getA()
            
        nEigenComponents = self.config.nEigenComponents # initial version

        if self.config.kernelSize >= 15:
            self.debugLog.debug(1, \
                "WARNING: NOT scaling kernelSize by stellar quadrupole moment, but using absolute value")
            actualKernelSize = int(self.config.kernelSize)
        else:
            actualKernelSize = 2 * int(self.config.kernelSize * numpy.sqrt(numpy.median(sizes)) + 0.5) + 1
            if actualKernelSize < self.config.kernelSizeMin:
                actualKernelSize = self.config.kernelSizeMin
            if actualKernelSize > self.config.kernelSizeMax:
                actualKernelSize = self.config.kernelSizeMax
            if display:
                print "Median size=%s" % (numpy.median(sizes),)
        self.debugLog.debug(3, "Kernel size=%s" % (actualKernelSize,))

        # Set size of image returned around candidate
        psfCandidateList[0].setHeight(actualKernelSize)
        psfCandidateList[0].setWidth(actualKernelSize)

        if display:
            frame = 0
            if displayExposure:
                ds9.mtv(exposure, frame=frame, title="psf determination")
                maUtils.showPsfSpatialCells(exposure, psfCellSet, self.config.nStarPerCell,
                                            symb="o", ctype=ds9.CYAN, ctypeUnused=ds9.YELLOW,
                                            size=4, frame=frame)
            
        #
        # Do a PCA decomposition of those PSF candidates
        #
        size = actualKernelSize + 2*self.config.borderWidth
        nu = size*size - 1                  # number of degrees of freedom/star for chi^2    
    
        reply = "y"                         # used in interactive mode
        for iter in range(self.config.nIterForPsf):
            if display and displayPsfCandidates: # Show a mosaic of usable PSF candidates
                #
                import lsst.afw.display.utils as displayUtils

                stamps = []
                for cell in psfCellSet.getCellList():
                    for cand in cell.begin(not showBadCandidates): # maybe include bad candidates
                        cand = algorithmsLib.cast_PsfCandidateF(cand)

                        try:
                            im = cand.getUndistImage().getImage()

                            chi2 = cand.getChi2()
                            if chi2 > 1e100:
                                chi2Str = ""
                            else:
                                chi2Str = " %.1f" % (chi2)

                            stamps.append((cand.getUndistImage().getImage(),
                                           "%d%s" % (cand.getSource().getId(), chi2Str),
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


                mos.makeMosaic(frame=7, title="Psf Candidates")

            # Re-fit until we don't have any candidates with naughty chi^2 values influencing the fit
            cleanChi2 = False # Any naughty (negative/NAN) chi^2 values?
            while not cleanChi2:
                cleanChi2 = True
                #
                # First, estimate the PSF
                #
                psf, eigenValues, nEigenComponents, fitChi2 = \
                    self._fitPsf(exposure, psfCellSet, actualKernelSize, nEigenComponents)
                #
                # In clipping, allow all candidates to be innocent until proven guilty on this iteration.
                # Throw out any prima facie guilty candidates (naughty chi^2 values)
                # 
                for cell in psfCellSet.getCellList():
                    awfulCandidates = []
                    for cand in cell.begin(False): # include bad candidates
                        cand = algorithmsLib.cast_PsfCandidateF(cand)
                        cand.setStatus(afwMath.SpatialCellCandidate.UNKNOWN) # until proven guilty
                        rchi2 = cand.getChi2()
                        if not numpy.isfinite(rchi2) or rchi2 <= 0:
                            # Guilty prima facie
                            awfulCandidates.append(cand)
                            cleanChi2 = False
                            self.debugLog.debug(2, "chi^2=%s; id=%s" %
                                                (cand.getChi2(), cand.getSource().getId()))
                    for cand in awfulCandidates:
                        if display:
                            print "Removing bad candidate: id=%d, chi^2=%f" % \
                                  (cand.getSource().getId(), cand.getChi2())
                        cell.removeCandidate(cand)

            #
            # Clip out bad fits based on raw chi^2
            #
            badCandidates = list()
            for cell in psfCellSet.getCellList():
                for cand in cell.begin(False): # include bad candidates
                    cand = algorithmsLib.cast_PsfCandidateF(cand)
                    rchi2 = cand.getChi2()  # reduced chi^2 when fitting PSF to candidate
                    assert rchi2 > 0
                    if rchi2 > self.config.reducedChi2ForPsfCandidates:
                        badCandidates.append(cand)

            badCandidates.sort(key=lambda x: x.getChi2(), reverse=True)
            numBad = int(len(badCandidates) * (iter + 1) / self.config.nIterForPsf + 0.5)
            for i, c in zip(range(numBad), badCandidates):
                if display:
                    print "Chi^2 clipping %d: %f" % (c.getSource().getId(), c.getChi2())
                c.setStatus(afwMath.SpatialCellCandidate.BAD)

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
                for cand in cell.begin(False):
                    cand = algorithmsLib.cast_PsfCandidateF(cand)
                    candCenter = afwGeom.PointD(cand.getXCenter(), cand.getYCenter())
                    try:
                        im = cand.getMaskedImage(kernel.getWidth(), kernel.getHeight())
                    except Exception, e:
                        continue

                    fit = algorithmsLib.fitKernelParamsToImage(noSpatialKernel, im, candCenter)
                    params = fit[0]
                    kernels = fit[1]
                    amp = 0.0
                    for p, k in zip(params, kernels):
                        amp += p * afwMath.cast_FixedKernel(k).getSum()

                    predict = [kernel.getSpatialFunction(k)(candCenter.getX(), candCenter.getY()) for
                               k in range(kernel.getNKernelParameters())]

                    #print cand.getSource().getId(), [a / amp for a in params], predict

                    residuals.append([a / amp - p for a, p in zip(params, predict)])
                    candidates.append(cand)

            residuals = numpy.array(residuals)

            for k in range(kernel.getNKernelParameters()):
                if False:
                    # Straight standard deviation
                    mean = residuals[:,k].mean()
                    rms = residuals[:,k].std()
                elif False:
                    # Using interquartile range
                    sr = numpy.sort(residuals[:,k])
                    mean = sr[int(0.5*len(sr))] if len(sr) % 2 else \
                           0.5 * (sr[int(0.5*len(sr))] + sr[int(0.5*len(sr))+1])
                    rms = 0.74 * (sr[int(0.75*len(sr))] - sr[int(0.25*len(sr))])
                else:
                    stats = afwMath.makeStatistics(residuals[:,k], afwMath.MEANCLIP | afwMath.STDEVCLIP)
                    mean = stats.getValue(afwMath.MEANCLIP)
                    rms = stats.getValue(afwMath.STDEVCLIP)

                rms = max(1.0e-4, rms)  # Don't trust RMS below this due to numerical issues

                if display:
                    print "Mean for component %d is %f" % (k, mean)
                    print "RMS for component %d is %f" % (k, rms)
                badCandidates = list()
                for i, cand in enumerate(candidates):
                    if numpy.fabs(residuals[i,k] - mean) > self.config.spatialReject * rms:
                        badCandidates.append(i)

                badCandidates.sort(key=lambda x: numpy.fabs(residuals[x,k] - mean), reverse=True)

                numBad = int(len(badCandidates) * (iter + 1) / self.config.nIterForPsf + 0.5)

                for i, c in zip(range(min(len(badCandidates), numBad)), badCandidates):
                    cand = candidates[c]
                    if display:
                        print "Spatial clipping %d (%f,%f) based on %d: %f vs %f" % \
                              (cand.getSource().getId(), cand.getXCenter(), cand.getYCenter(), k,
                               residuals[badCandidates[i],k], self.config.spatialReject * rms)
                    cand.setStatus(afwMath.SpatialCellCandidate.BAD)

            #
            # Display results
            #
            if display and displayIterations:
                if displayExposure:
                    if iter > 0:
                        ds9.erase(frame=frame)
                    maUtils.showPsfSpatialCells(exposure, psfCellSet, self.config.nStarPerCell, showChi2=True,
                                                symb="o", size=8, frame=frame,
                                                ctype=ds9.YELLOW, ctypeBad=ds9.RED, ctypeUnused=ds9.MAGENTA)
                    if self.config.nStarPerCellSpatialFit != self.config.nStarPerCell:
                        maUtils.showPsfSpatialCells(exposure, psfCellSet, self.config.nStarPerCellSpatialFit,
                                                    symb="o", size=10, frame=frame,
                                                    ctype=ds9.YELLOW, ctypeBad=ds9.RED)
                while True:
                    try:
                        maUtils.showPsfCandidates(exposure, psfCellSet, psf=psf, frame=4,
                                                  normalize=normalizeResiduals,
                                                  showBadCandidates=showBadCandidates)
                        maUtils.showPsfCandidates(exposure, psfCellSet, psf=psf, frame=5,
                                                  normalize=normalizeResiduals,
                                                  showBadCandidates=showBadCandidates,
                                                  variance=True)
                    except:
                        if not showBadCandidates:
                            showBadCandidates = True
                            continue
                    break

                if displayPsfComponents:
                    maUtils.showPsf(psf, eigenValues, frame=6)
                if displayPsfMosaic:
                    maUtils.showPsfMosaic(exposure, psf, frame=7)
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
        psf, eigenValues, nEigenComponents, fitChi2 = \
            self._fitPsf(exposure, psfCellSet, actualKernelSize, nEigenComponents)

        #
        # Display code for debugging
        #
        if display and reply != "n":
            if displayExposure:
                maUtils.showPsfSpatialCells(exposure, psfCellSet, self.config.nStarPerCell, showChi2=True,
                                            symb="o", ctype=ds9.YELLOW, ctypeBad=ds9.RED, size=8, frame=frame)
                if self.config.nStarPerCellSpatialFit != self.config.nStarPerCell:
                    maUtils.showPsfSpatialCells(exposure, psfCellSet, self.config.nStarPerCellSpatialFit,
                                                symb="o", ctype=ds9.YELLOW, ctypeBad=ds9.RED,
                                                size=10, frame=frame)
            maUtils.showPsfCandidates(exposure, psfCellSet, psf=psf, frame=4, normalize=normalizeResiduals,
                                      showBadCandidates=showBadCandidates)

            maUtils.showPsf(psf, eigenValues, frame=6)
            if displayPsfMosaic:
                maUtils.showPsfMosaic(exposure, psf, frame=7)
            if displayPsfSpatialModel:
                maUtils.plotPsfSpatialModel(exposure, psf, psfCellSet, showBadCandidates=True,
                                            matchKernelAmplitudes=matchKernelAmplitudes,
                                            keepPlots=keepMatplotlibPlots)
        #
        # Generate some QA information
        #
        # Count PSF stars
        #
        numGoodStars = 0
        numAvailStars = 0

        avgX = 0.0
        avgY = 0.0

        for cell in psfCellSet.getCellList():
            for cand in cell.begin(False):  # don't ignore BAD stars
                numAvailStars += 1

            for cand in cell.begin(True):  # do ignore BAD stars
                cand = algorithmsLib.cast_PsfCandidateF(cand)
                src = cand.getSource()
                if flagKey is not None:
                    src.set(flagKey, True)
                avgX += src.getX()
                avgY += src.getY()
                numGoodStars += 1

        avgX /= numGoodStars
        avgY /= numGoodStars

        if metadata != None:
            metadata.set("spatialFitChi2", fitChi2)
            metadata.set("numGoodStars", numGoodStars)
            metadata.set("numAvailStars", numAvailStars)
            metadata.set("avgX", avgX)
            metadata.set("avgY", avgY)

        psf = algorithmsLib.PcaPsf(psf.getKernel(), afwGeom.Point2D(avgX, avgY))

        return psf, psfCellSet
