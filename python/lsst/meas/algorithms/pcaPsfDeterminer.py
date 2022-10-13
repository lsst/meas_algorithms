# This file is part of meas_algorithms.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

__all__ = ["PcaPsfDeterminerConfig", "PcaPsfDeterminerTask"]

import math
import sys

import numpy

import lsst.pex.config as pexConfig
import lsst.pex.exceptions as pexExceptions
import lsst.geom
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses as afwEll
import lsst.afw.display as afwDisplay
import lsst.afw.math as afwMath
from .psfDeterminer import BasePsfDeterminerTask, psfDeterminerRegistry
from .psfCandidate import PsfCandidateF
from .spatialModelPsf import createKernelFromPsfCandidates, countPsfCandidates, \
    fitSpatialKernelFromPsfCandidates, fitKernelParamsToImage
from .pcaPsf import PcaPsf
from . import utils


def numCandidatesToReject(numBadCandidates, numIter, totalIter):
    """Return the number of PSF candidates to be rejected.

    The number of candidates being rejected on each iteration gradually
    increases, so that on the Nth of M iterations we reject N/M of the bad
    candidates.

    Parameters
    ----------
    numBadCandidates : `int`
        Number of bad candidates under consideration.

    numIter : `int`
        The number of the current PSF iteration.

    totalIter : `int`
        The total number of PSF iterations.

    Returns
    -------
    return : `int`
        Number of candidates to reject.
    """
    return int(numBadCandidates*(numIter + 1)//totalIter + 0.5)


class PcaPsfDeterminerConfig(BasePsfDeterminerTask.ConfigClass):
    nonLinearSpatialFit = pexConfig.Field[bool](
        doc="Use non-linear fitter for spatial variation of Kernel",
        default=False,
    )
    nEigenComponents = pexConfig.Field[int](
        doc="number of eigen components for PSF kernel creation",
        default=4,
        check=lambda x: x >= 1
    )
    spatialOrder = pexConfig.Field[int](
        doc="specify spatial order for PSF kernel creation",
        default=2,
    )
    sizeCellX = pexConfig.Field[int](
        doc="size of cell used to determine PSF (pixels, column direction)",
        default=256,
        #        minValue = 10,
        check=lambda x: x >= 10,
    )
    sizeCellY = pexConfig.Field[int](
        doc="size of cell used to determine PSF (pixels, row direction)",
        default=sizeCellX.default,
        #        minValue = 10,
        check=lambda x: x >= 10,
    )
    nStarPerCell = pexConfig.Field[int](
        doc="number of stars per psf cell for PSF kernel creation",
        default=3,
    )
    borderWidth = pexConfig.Field[int](
        doc="Number of pixels to ignore around the edge of PSF candidate postage stamps",
        default=0,
    )
    nStarPerCellSpatialFit = pexConfig.Field[int](
        doc="number of stars per psf Cell for spatial fitting",
        default=5,
    )
    constantWeight = pexConfig.Field[bool](
        doc="Should each PSF candidate be given the same weight, independent of magnitude?",
        default=True,
    )
    nIterForPsf = pexConfig.Field[int](
        doc="number of iterations of PSF candidate star list",
        default=3,
    )
    tolerance = pexConfig.Field[float](
        doc="tolerance of spatial fitting",
        default=1e-2,
    )
    lam = pexConfig.Field[float](
        doc="floor for variance is lam*data",
        default=0.05,
    )
    reducedChi2ForPsfCandidates = pexConfig.Field[float](
        doc="for psf candidate evaluation",
        default=2.0,
    )
    spatialReject = pexConfig.Field[float](
        doc="Rejection threshold (stdev) for candidates based on spatial fit",
        default=3.0,
    )
    pixelThreshold = pexConfig.Field[float](
        doc="Threshold (stdev) for rejecting extraneous pixels around candidate; applied if positive",
        default=0.0,
    )
    doRejectBlends = pexConfig.Field[bool](
        doc="Reject candidates that are blended?",
        default=False,
    )
    doMaskBlends = pexConfig.Field[bool](
        doc="Mask blends in image?",
        default=True,
    )

    def setDefaults(self):
        super().setDefaults()
        self.stampSize = 41


class PcaPsfDeterminerTask(BasePsfDeterminerTask):
    """A measurePsfTask psf estimator.
    """
    ConfigClass = PcaPsfDeterminerConfig

    def _fitPsf(self, exposure, psfCellSet, kernelSize, nEigenComponents):
        PsfCandidateF.setPixelThreshold(self.config.pixelThreshold)
        PsfCandidateF.setMaskBlends(self.config.doMaskBlends)
        #
        # Loop trying to use nEigenComponents, but allowing smaller numbers if necessary
        #
        for nEigen in range(nEigenComponents, 0, -1):
            # Determine KL components
            try:
                kernel, eigenValues = createKernelFromPsfCandidates(
                    psfCellSet, exposure.getDimensions(), exposure.getXY0(), nEigen,
                    self.config.spatialOrder, kernelSize, self.config.nStarPerCell,
                    bool(self.config.constantWeight))

                break                   # OK, we can get nEigen components
            except pexExceptions.LengthError as e:
                if nEigen == 1:         # can't go any lower
                    raise IndexError("No viable PSF candidates survive")

                self.log.warning("%s: reducing number of eigen components", e.what())
        #
        # We got our eigen decomposition so let's use it
        #
        # Express eigenValues in units of reduced chi^2 per star
        size = kernelSize + 2*self.config.borderWidth
        nu = size*size - 1                  # number of degrees of freedom/star for chi^2
        eigenValues = [val/float(countPsfCandidates(psfCellSet, self.config.nStarPerCell)*nu)
                       for val in eigenValues]

        # Fit spatial model
        status, chi2 = fitSpatialKernelFromPsfCandidates(
            kernel, psfCellSet, bool(self.config.nonLinearSpatialFit),
            self.config.nStarPerCellSpatialFit, self.config.tolerance, self.config.lam)

        psf = PcaPsf(kernel)

        return psf, eigenValues, nEigen, chi2

    def determinePsf(self, exposure, psfCandidateList, metadata=None, flagKey=None):
        """Determine a PCA PSF model for an exposure given a list of PSF candidates.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
           Exposure containing the psf candidates.
        psfCandidateList : `list` of `lsst.meas.algorithms.PsfCandidate`
           A sequence of PSF candidates typically obtained by detecting sources
           and then running them through a star selector.
        metadata : `lsst.daf.base import PropertyList` or `None`, optional
           A home for interesting tidbits of information.
        flagKey : `str`, optional
           Schema key used to mark sources actually used in PSF determination.

        Returns
        -------
        psf : `lsst.meas.algorithms.PcaPsf`
           The measured PSF.
        psfCellSet : `lsst.afw.math.SpatialCellSet`
           The PSF candidates.
        """
        import lsstDebug
        display = lsstDebug.Info(__name__).display
        displayExposure = lsstDebug.Info(__name__).displayExposure     # display the Exposure + spatialCells
        displayPsfCandidates = lsstDebug.Info(__name__).displayPsfCandidates  # show the viable candidates
        displayIterations = lsstDebug.Info(__name__).displayIterations  # display on each PSF iteration
        displayPsfComponents = lsstDebug.Info(__name__).displayPsfComponents  # show the PCA components
        displayResiduals = lsstDebug.Info(__name__).displayResiduals         # show residuals
        displayPsfMosaic = lsstDebug.Info(__name__).displayPsfMosaic   # show mosaic of reconstructed PSF(x,y)
        # match Kernel amplitudes for spatial plots
        matchKernelAmplitudes = lsstDebug.Info(__name__).matchKernelAmplitudes
        # Keep matplotlib alive post mortem
        keepMatplotlibPlots = lsstDebug.Info(__name__).keepMatplotlibPlots
        displayPsfSpatialModel = lsstDebug.Info(__name__).displayPsfSpatialModel  # Plot spatial model?
        showBadCandidates = lsstDebug.Info(__name__).showBadCandidates  # Include bad candidates
        # Normalize residuals by object amplitude
        normalizeResiduals = lsstDebug.Info(__name__).normalizeResiduals
        pause = lsstDebug.Info(__name__).pause                         # Prompt user after each iteration?

        if display:
            afwDisplay.setDefaultMaskTransparency(75)
        if display > 1:
            pause = True

        mi = exposure.getMaskedImage()

        if len(psfCandidateList) == 0:
            raise RuntimeError("No PSF candidates supplied.")

        # construct and populate a spatial cell set
        bbox = mi.getBBox()
        psfCellSet = afwMath.SpatialCellSet(bbox, self.config.sizeCellX, self.config.sizeCellY)
        sizes = []
        for i, psfCandidate in enumerate(psfCandidateList):
            if psfCandidate.getSource().getPsfFluxFlag():  # bad measurement
                continue

            try:
                psfCellSet.insertCandidate(psfCandidate)
            except Exception as e:
                self.log.debug("Skipping PSF candidate %d of %d: %s", i, len(psfCandidateList), e)
                continue
            source = psfCandidate.getSource()

            quad = afwGeom.Quadrupole(source.getIxx(), source.getIyy(), source.getIxy())
            axes = afwEll.Axes(quad)
            sizes.append(axes.getA())
        if len(sizes) == 0:
            raise RuntimeError("No usable PSF candidates supplied")
        nEigenComponents = self.config.nEigenComponents  # initial version

        # TODO: DM-36311: Keep only the if block below.
        if self.config.stampSize:
            actualKernelSize = int(self.config.stampSize)
        elif self.config.kernelSize >= 15:
            self.log.warning("NOT scaling kernelSize by stellar quadrupole moment "
                             "because config.kernelSize=%s >= 15; "
                             "using config.kernelSize as the width, instead",
                             self.config.kernelSize)
            actualKernelSize = int(self.config.kernelSize)
        else:
            self.log.warning("scaling kernelSize by stellar quadrupole moment "
                             "because config.kernelSize=%s < 15. This behavior is deprecated.",
                             self.config.kernelSize)
            medSize = numpy.median(sizes)
            actualKernelSize = 2*int(self.config.kernelSize*math.sqrt(medSize) + 0.5) + 1
            if actualKernelSize < self.config.kernelSizeMin:
                actualKernelSize = self.config.kernelSizeMin
            if actualKernelSize > self.config.kernelSizeMax:
                actualKernelSize = self.config.kernelSizeMax

            if display:
                print("Median size=%s" % (medSize,))
        self.log.trace("Kernel size=%s", actualKernelSize)

        if actualKernelSize > psfCandidateList[0].getWidth():
            self.log.warning("Using a region (%d x %d) larger than kernelSize (%d) set while making PSF "
                             "candidates. Consider setting a larger value for kernelSize for "
                             "`makePsfCandidates` to avoid this warning.",
                             actualKernelSize, actualKernelSize, psfCandidateList[0].getWidth())

        if self.config.doRejectBlends:
            # Remove blended candidates completely
            blendedCandidates = []  # Candidates to remove; can't do it while iterating
            for cell, cand in candidatesIter(psfCellSet, False):
                if len(cand.getSource().getFootprint().getPeaks()) > 1:
                    blendedCandidates.append((cell, cand))
                    continue
            if display:
                print("Removing %d blended Psf candidates" % len(blendedCandidates))
            for cell, cand in blendedCandidates:
                cell.removeCandidate(cand)
            if sum(1 for cand in candidatesIter(psfCellSet, False)) == 0:
                raise RuntimeError("All PSF candidates removed as blends")

        if display:
            if displayExposure:
                disp = afwDisplay.Display(frame=0)
                disp.mtv(exposure, title="psf determination")
                utils.showPsfSpatialCells(exposure, psfCellSet, self.config.nStarPerCell, symb="o",
                                          ctype=afwDisplay.CYAN, ctypeUnused=afwDisplay.YELLOW,
                                          size=4, display=disp)

        #
        # Do a PCA decomposition of those PSF candidates
        #
        reply = "y"                         # used in interactive mode
        for iterNum in range(self.config.nIterForPsf):
            if display and displayPsfCandidates:  # Show a mosaic of usable PSF candidates

                stamps = []
                for cell in psfCellSet.getCellList():
                    for cand in cell.begin(not showBadCandidates):  # maybe include bad candidates
                        try:
                            im = cand.getMaskedImage()

                            chi2 = cand.getChi2()
                            if chi2 > 1e100:
                                chi2 = numpy.nan

                            stamps.append((im, "%d%s" %
                                           (utils.splitId(cand.getSource().getId(), True)["objId"], chi2),
                                           cand.getStatus()))
                        except Exception:
                            continue

                if len(stamps) == 0:
                    print("WARNING: No PSF candidates to show; try setting showBadCandidates=True")
                else:
                    mos = afwDisplay.utils.Mosaic()
                    for im, label, status in stamps:
                        im = type(im)(im, True)
                        try:
                            im /= afwMath.makeStatistics(im, afwMath.MAX).getValue()
                        except NotImplementedError:
                            pass

                        mos.append(im, label,
                                   (afwDisplay.GREEN if status == afwMath.SpatialCellCandidate.GOOD else
                                    afwDisplay.YELLOW if status == afwMath.SpatialCellCandidate.UNKNOWN else
                                    afwDisplay.RED))

                    disp8 = afwDisplay.Display(frame=8)
                    mos.makeMosaic(display=disp8, title="Psf Candidates")

            # Re-fit until we don't have any candidates with naughty chi^2 values influencing the fit
            cleanChi2 = False  # Any naughty (negative/NAN) chi^2 values?
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
                    for cand in cell.begin(False):  # include bad candidates
                        cand.setStatus(afwMath.SpatialCellCandidate.UNKNOWN)  # until proven guilty
                        rchi2 = cand.getChi2()
                        if not numpy.isfinite(rchi2) or rchi2 <= 0:
                            # Guilty prima facie
                            awfulCandidates.append(cand)
                            cleanChi2 = False
                            self.log.debug("chi^2=%s; id=%s",
                                           cand.getChi2(), cand.getSource().getId())
                    for cand in awfulCandidates:
                        if display:
                            print("Removing bad candidate: id=%d, chi^2=%f" %
                                  (cand.getSource().getId(), cand.getChi2()))
                        cell.removeCandidate(cand)

            #
            # Clip out bad fits based on reduced chi^2
            #
            badCandidates = list()
            for cell in psfCellSet.getCellList():
                for cand in cell.begin(False):  # include bad candidates
                    rchi2 = cand.getChi2()  # reduced chi^2 when fitting PSF to candidate
                    assert rchi2 > 0
                    if rchi2 > self.config.reducedChi2ForPsfCandidates:
                        badCandidates.append(cand)

            badCandidates.sort(key=lambda x: x.getChi2(), reverse=True)
            numBad = numCandidatesToReject(len(badCandidates), iterNum,
                                           self.config.nIterForPsf)
            for i, c in zip(range(numBad), badCandidates):
                if display:
                    chi2 = c.getChi2()
                    if chi2 > 1e100:
                        chi2 = numpy.nan

                    print("Chi^2 clipping %-4d  %.2g" % (c.getSource().getId(), chi2))
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
            noSpatialKernel = psf.getKernel()
            for cell in psfCellSet.getCellList():
                for cand in cell.begin(False):
                    candCenter = lsst.geom.PointD(cand.getXCenter(), cand.getYCenter())
                    try:
                        im = cand.getMaskedImage(actualKernelSize, actualKernelSize)
                    except Exception:
                        continue

                    fit = fitKernelParamsToImage(noSpatialKernel, im, candCenter)
                    params = fit[0]
                    kernels = fit[1]
                    amp = 0.0
                    for p, k in zip(params, kernels):
                        amp += p*k.getSum()

                    predict = [kernel.getSpatialFunction(k)(candCenter.getX(), candCenter.getY()) for
                               k in range(kernel.getNKernelParameters())]

                    residuals.append([a/amp - p for a, p in zip(params, predict)])
                    candidates.append(cand)

            residuals = numpy.array(residuals)

            for k in range(kernel.getNKernelParameters()):
                if False:
                    # Straight standard deviation
                    mean = residuals[:, k].mean()
                    rms = residuals[:, k].std()
                elif False:
                    # Using interquartile range
                    sr = numpy.sort(residuals[:, k])
                    mean = (sr[int(0.5*len(sr))] if len(sr)%2 else
                            0.5*(sr[int(0.5*len(sr))] + sr[int(0.5*len(sr)) + 1]))
                    rms = 0.74*(sr[int(0.75*len(sr))] - sr[int(0.25*len(sr))])
                else:
                    stats = afwMath.makeStatistics(residuals[:, k], afwMath.MEANCLIP | afwMath.STDEVCLIP)
                    mean = stats.getValue(afwMath.MEANCLIP)
                    rms = stats.getValue(afwMath.STDEVCLIP)

                rms = max(1.0e-4, rms)  # Don't trust RMS below this due to numerical issues

                if display:
                    print("Mean for component %d is %f" % (k, mean))
                    print("RMS for component %d is %f" % (k, rms))
                badCandidates = list()
                for i, cand in enumerate(candidates):
                    if numpy.fabs(residuals[i, k] - mean) > self.config.spatialReject*rms:
                        badCandidates.append(i)

                badCandidates.sort(key=lambda x: numpy.fabs(residuals[x, k] - mean), reverse=True)

                numBad = numCandidatesToReject(len(badCandidates), iterNum,
                                               self.config.nIterForPsf)

                for i, c in zip(range(min(len(badCandidates), numBad)), badCandidates):
                    cand = candidates[c]
                    if display:
                        print("Spatial clipping %d (%f,%f) based on %d: %f vs %f" %
                              (cand.getSource().getId(), cand.getXCenter(), cand.getYCenter(), k,
                               residuals[badCandidates[i], k], self.config.spatialReject*rms))
                    cand.setStatus(afwMath.SpatialCellCandidate.BAD)

            #
            # Display results
            #
            if display and displayIterations:
                if displayExposure:
                    if iterNum > 0:
                        disp.erase()
                    utils.showPsfSpatialCells(exposure, psfCellSet, self.config.nStarPerCell, showChi2=True,
                                              symb="o", size=8, display=disp, ctype=afwDisplay.YELLOW,
                                              ctypeBad=afwDisplay.RED, ctypeUnused=afwDisplay.MAGENTA)
                    if self.config.nStarPerCellSpatialFit != self.config.nStarPerCell:
                        utils.showPsfSpatialCells(exposure, psfCellSet, self.config.nStarPerCellSpatialFit,
                                                  symb="o", size=10, display=disp,
                                                  ctype=afwDisplay.YELLOW, ctypeBad=afwDisplay.RED)
                if displayResiduals:
                    while True:
                        try:
                            disp4 = afwDisplay.Display(frame=4)
                            utils.showPsfCandidates(exposure, psfCellSet, psf=psf, display=disp4,
                                                    normalize=normalizeResiduals,
                                                    showBadCandidates=showBadCandidates)
                            disp5 = afwDisplay.Display(frame=5)
                            utils.showPsfCandidates(exposure, psfCellSet, psf=psf, display=disp5,
                                                    normalize=normalizeResiduals,
                                                    showBadCandidates=showBadCandidates,
                                                    variance=True)
                        except Exception:
                            if not showBadCandidates:
                                showBadCandidates = True
                                continue
                        break

                if displayPsfComponents:
                    disp6 = afwDisplay.Display(frame=6)
                    utils.showPsf(psf, eigenValues, display=disp6)
                if displayPsfMosaic:
                    disp7 = afwDisplay.Display(frame=7)
                    utils.showPsfMosaic(exposure, psf, display=disp7, showFwhm=True)
                    disp7.scale('linear', 0, 1)
                if displayPsfSpatialModel:
                    utils.plotPsfSpatialModel(exposure, psf, psfCellSet, showBadCandidates=True,
                                              matchKernelAmplitudes=matchKernelAmplitudes,
                                              keepPlots=keepMatplotlibPlots)

                if pause:
                    while True:
                        try:
                            reply = input("Next iteration? [ynchpqQs] ").strip()
                        except EOFError:
                            reply = "n"

                        reply = reply.split()
                        if reply:
                            reply, args = reply[0], reply[1:]
                        else:
                            reply = ""

                        if reply in ("", "c", "h", "n", "p", "q", "Q", "s", "y"):
                            if reply == "c":
                                pause = False
                            elif reply == "h":
                                print("c[ontinue without prompting] h[elp] n[o] p[db] q[uit displaying] "
                                      "s[ave fileName] y[es]")
                                continue
                            elif reply == "p":
                                import pdb
                                pdb.set_trace()
                            elif reply == "q":
                                display = False
                            elif reply == "Q":
                                sys.exit(1)
                            elif reply == "s":
                                fileName = args.pop(0)
                                if not fileName:
                                    print("Please provide a filename")
                                    continue

                                print("Saving to %s" % fileName)
                                utils.saveSpatialCellSet(psfCellSet, fileName=fileName)
                                continue
                            break
                        else:
                            print("Unrecognised response: %s" % reply, file=sys.stderr)

                    if reply == "n":
                        break

        # One last time, to take advantage of the last iteration
        psf, eigenValues, nEigenComponents, fitChi2 = \
            self._fitPsf(exposure, psfCellSet, actualKernelSize, nEigenComponents)

        #
        # Display code for debugging
        #
        if display and reply != "n":
            disp = afwDisplay.Display(frame=0)
            if displayExposure:
                utils.showPsfSpatialCells(exposure, psfCellSet, self.config.nStarPerCell, showChi2=True,
                                          symb="o", ctype=afwDisplay.YELLOW, ctypeBad=afwDisplay.RED,
                                          size=8, display=disp)
                if self.config.nStarPerCellSpatialFit != self.config.nStarPerCell:
                    utils.showPsfSpatialCells(exposure, psfCellSet, self.config.nStarPerCellSpatialFit,
                                              symb="o", ctype=afwDisplay.YELLOW, ctypeBad=afwDisplay.RED,
                                              size=10, display=disp)
                if displayResiduals:
                    disp4 = afwDisplay.Display(frame=4)
                    utils.showPsfCandidates(exposure, psfCellSet, psf=psf, display=disp4,
                                            normalize=normalizeResiduals,
                                            showBadCandidates=showBadCandidates)

            if displayPsfComponents:
                disp6 = afwDisplay.Display(frame=6)
                utils.showPsf(psf, eigenValues, display=disp6)

            if displayPsfMosaic:
                disp7 = afwDisplay.Display(frame=7)
                utils.showPsfMosaic(exposure, psf, display=disp7, showFwhm=True)
                disp7.scale("linear", 0, 1)
            if displayPsfSpatialModel:
                utils.plotPsfSpatialModel(exposure, psf, psfCellSet, showBadCandidates=True,
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
                src = cand.getSource()
                if flagKey is not None:
                    src.set(flagKey, True)
                avgX += src.getX()
                avgY += src.getY()
                numGoodStars += 1

        avgX /= numGoodStars
        avgY /= numGoodStars

        if metadata is not None:
            metadata["spatialFitChi2"] = fitChi2
            metadata["numGoodStars"] = numGoodStars
            metadata["numAvailStars"] = numAvailStars
            metadata["avgX"] = avgX
            metadata["avgY"] = avgY

        psf = PcaPsf(psf.getKernel(), lsst.geom.Point2D(avgX, avgY))

        return psf, psfCellSet


def candidatesIter(psfCellSet, ignoreBad=True):
    """Generator for Psf candidates.

    This allows two 'for' loops to be reduced to one.

    Parameters
    ----------
    psfCellSet : `lsst.afw.math.SpatialCellSet`
       SpatialCellSet of PSF candidates.
    ignoreBad : `bool`, optional
       Ignore candidates flagged as BAD?

    Yields
    -------
    cell : `lsst.afw.math.SpatialCell`
       A SpatialCell.
    cand : `lsst.meas.algorithms.PsfCandidate`
      A PsfCandidate.
    """
    for cell in psfCellSet.getCellList():
        for cand in cell.begin(ignoreBad):
            yield (cell, cand)


psfDeterminerRegistry.register("pca", PcaPsfDeterminerTask)
