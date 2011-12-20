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
import lsst.pex.config as pexConfig
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses as afwEll
import lsst.afw.detection as afwDetection
import lsst.afw.display.ds9 as ds9
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import algorithmsLib
import utils as maUtils
import numpy

class PcaPsfDeteminerConfig(pexConfig.Config):
    nonLinearSpatialFit = pexConfig.Field(
        dtype = bool,
        doc = "Use non-linear fitter for spatial variation of Kernel",
        default = False,
        optional = False,
    )
    nEigenComponents = pexConfig.Field(
        dtype = int,
        doc = "number of eigen components for PSF kernel creation",
        default = 3,
        optional = False,
    )
    spatialOrder = pexConfig.Field(
        dtype = int,
        doc = "specify spatial order for PSF kernel creation",
        default = 2,
        optional = False,
    )
    sizeCellX = pexConfig.Field(
        dtype = int,
        doc = "size of cell used to determine PSF (pixels, column direction)",
        default = 256,
        minValue = 10,
        optional = False,
    )
    sizeCellY = pexConfig.Field(
        dtype = int,
        doc = "size of cell used to determine PSF (pixels, row direction)",
        default = 256,
        minValue = 10,
        optional = False,
    )
    nStarPerCell = pexConfig.Field(
        dtype = int,
        doc = "number of stars per psf cell for PSF kernel creation",
        default = 3,
        optional = False,
    )
    kernelSize = pexConfig.Field(
        dtype = float,
        doc = "radius of the kernel to create, relative to the square root of the stellar quadrupole moments",
        default = 5.0,
        optional = False,
    )
    kernelSizeMin = pexConfig.Field(
        dtype = int,
        doc = "Minimum radius of the kernel",
        default = 13,
        optional = False,
    )
    kernelSizeMax = pexConfig.Field(
        dtype = int,
        doc = "Maximum radius of the kernel",
        default = 45,
        optional = False,
    )
    borderWidth = pexConfig.Field(
        dtype = int,
        doc = "Number of pixels to ignore around the edge of PSF candidate postage stamps",
        default = 0,
        optional = False,
    )
    nStarPerCellSpatialFit = pexConfig.Field(
        dtype = int,
        doc = "number of stars per psf Cell for spatial fitting",
        default = 5,
        optional = False,
    )
    constantWeight = pexConfig.Field(
        dtype = bool,
        doc = "Should each PSF candidate be given the same weight, independent of magnitude?",
        default = True,
        optional = False,
    )
    nIterForPsf = pexConfig.Field(
        dtype = int,
        doc = "number of iterations of PSF candidate star list",
        default = 3,
        optional = False,
    )
    tolerance = pexConfig.Field(
        dtype = float,
        doc = "tolerance of spatial fitting",
        default = 1e-2,
        optional = False,
    )
    lam = pexConfig.Field(
        dtype = float,
        doc = "floor for variance is lam*data",
        default = 0.05,
        optional = False,
    )
    reducedChi2ForPsfCandidates = pexConfig.Field(
        dtype = float,
        doc = "for psf candidate evaluation",
        default = 2.0,
        optional = False,
    )
    spatialReject = pexConfig.Field(
        dtype = float,
        doc = "Rejection threshold (stdev) for candidates based on spatial fit",
        default = 3.0,
        optional = False,
    )


class PcaPsfDeterminer(object):
    ConfigClass = PcaPsfDeteminerConfig

    def __init__(self, config):
        """Construct a PCA PSF Fitter

        @param[in] config: instance of self.ConfigClass = PcaPsfDeteminerConfig
        """
        self.config = config

    def _fitPsf(self, exposure, psfCellSet):
        # Determine KL components
        kernel, eigenValues = algorithmsLib.createKernelFromPsfCandidates(
            psfCellSet, exposure.getDimensions(), self.config.nEigenComponents, self.spatialOrder,
            self.config.kernelSize, self.config.nStarPerCell, self.config.constantWeight)

        # Express eigenValues in units of reduced chi^2 per star
        size = self.config.kernelSize + 2*self.config.borderWidth
        nu = size*size - 1                  # number of degrees of freedom/star for chi^2    
        eigenValues = [l/float(algorithmsLib.countPsfCandidates(psfCellSet, self.config.nStarPerCell)*nu)
                       for l in eigenValues]
        
        # Fit spatial model
        status, chi2 = algorithmsLib.fitSpatialKernelFromPsfCandidates(
            kernel, psfCellSet, self.config.nonLinearSpatialFit,
            self.config.nStarPerCellSpatialFit, self.config.tolerance, self.config.lam)
        
        psf = afwDetection.createPsf("PCA", kernel)

        return psf, eigenValues, chi2


    def determinePsf(self, exposure, psfCandidateList, metadata=None):
        """Determine a PCA PSF model for an exposure given a list of PSF candidates
        
        @param[in] exposure: exposure containing the psf candidates (lsst.afw.image.Exposure)
        @param[in] psfCandidateList: a sequence of PSF candidates (each an lsst.meas.algorithms.PsfCandidate);
            typically obtained by detecting sources and then running them through a star selector
        @param[in,out] metadata  a home for interesting tidbits of information
    
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
                print e
                continue
            source = psfCandidate.getSource()

            quad = afwEll.Quadrupole(source.getIxx(), source.getIyy(), source.getIxy())
            axes = afwEll.Axes(quad)
            sizes[i] = axes.getA()

        if self.config.kernelSize >= 15:
            print "WARNING: NOT scaling kernelSize by stellar quadrupole moment, but using absolute value"
            self.config.kernelSize = int(self.config.kernelSize)
        else:
            self.config.kernelSize = 2 * int(self.config.kernelSize * numpy.sqrt(numpy.median(sizes)) + 0.5) + 1
            if self.config.kernelSize < self.config.kernelSizeMin:
                self.config.kernelSize = self.config.kernelSizeMin
            if self.config.kernelSize > self.config.kernelSizeMax:
                self.config.kernelSize = self.config.kernelSizeMax
            if display:
                print "Median size:", numpy.median(sizes)
                print "Kernel size:", self.config.kernelSize

        # Set size of image returned around candidate
        psfCandidateList[0].setHeight(self.config.kernelSize)
        psfCandidateList[0].setWidth(self.config.kernelSize)

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
        size = self.config.kernelSize + 2*self.config.borderWidth
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
            psf, eigenValues, fitChi2 = self._fitPsf(exposure, psfCellSet)

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
            badCandidates = list()
            for cell in psfCellSet.getCellList():
                for cand in cell.begin(False): # include bad candidates
                    cand = algorithmsLib.cast_PsfCandidateF(cand)
                    rchi2 = cand.getChi2()  # reduced chi^2 when fitting PSF to candidate
                    if rchi2 < 0 or rchi2 > self.config.reducedChi2ForPsfCandidates or numpy.isnan(rchi2):
                        badCandidates.append(cand)
                        if rchi2 < 0:
                            print "RHL chi^2:", cand.getChi2(), nu, cand.getSource().getId()

            badCandidates.sort(key=lam x: x.getChi2(), reverse=True)
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
                        im = cand.getImage()
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
                    if numpy.fabs(residuals[i,k] - mean) > self.spatialReject * rms:
                        badCandidates.append(i)

                badCandidates.sort(key=lam x: numpy.fabs(residuals[x,k] - mean), reverse=True)

                numBad = int(len(badCandidates) * (iter + 1) / self.config.nIterForPsf + 0.5)
                
                for i, c in zip(range(min(len(badCandidates), numBad)), badCandidates):
                    cand = candidates[c]
                    if display:
                        print "Spatial clipping %d (%f,%f) based on %d: %f vs %f" % \
                              (cand.getSource().getId(), cand.getXCenter(), cand.getYCenter(), k,
                               residuals[badCandidates[i],k], self.spatialReject * rms)
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
        psf, eigenValues, fitChi2 = self._fitPsf(exposure, psfCellSet)

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
    
        for cell in psfCellSet.getCellList():
            for cand in cell.begin(False):  # don't ignore BAD stars
                numAvailStars += 1

            for cand in cell.begin(True):  # do ignore BAD stars
                cand = algorithmsLib.cast_PsfCandidateF(cand)
                src = cand.getSource()
                src.setFlagForDetection(src.getFlagForDetection() | algorithmsLib.Flags.PSFSTAR)
                numGoodStars += 1
    
        if metadata != None:
            metadata.set("spatialFitChi2", fitChi2)
            metadata.set("numGoodStars", numGoodStars)
            metadata.set("numAvailStars", numAvailStars)
    
        return psf, psfCellSet
