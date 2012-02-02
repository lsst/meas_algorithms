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

"""Support utilities for Measuring sources"""

import re
import sys

import numpy

import lsst.pex.exceptions as pexExcept
import lsst.daf.base as dafBase
import lsst.afw.detection as afwDet
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.display.ds9 as ds9
import lsst.afw.display.utils as displayUtils
import algorithmsLib

keptPlots = False                       # Have we arranged to keep spatial plots open?


def explainDetectionFlags(flags):
    """Return a string explaining Source's detectionFlags"""

    result = []
    for k, v in getDetectionFlags().items():
        if (flags & v):
            result += [k]

    result.sort()
    return " ".join(result)
    
def getDetectionFlags(key=None):
    """Return a dictionary of Source's detectionFlags"""

    flags = {}
    for k in algorithmsLib.Flags.__dict__.keys():
        if not re.search(r"^[_A-Z0-9]+$", k): # flag names match this re
            continue

        flags[k] = algorithmsLib.Flags.__dict__[k]

    if key:
        return flags.get(key)
    else:
        return flags
    
def showSourceSet(sSet, xy0=(0, 0), frame=0, ctype=ds9.GREEN, symb="+", size=2):
    """Draw the (XAstrom, YAstrom) positions of a set of Sources.  Image has the given XY0"""
    ds9.cmdBuffer.pushSize()

    for s in sSet:
        xc, yc = s.getXAstrom() - xy0[0], s.getYAstrom() - xy0[1]
        
        if symb == "id":
            ds9.dot(str(s.getId()), xc, yc, frame=frame, ctype=ctype, size=size)
        else:
            ds9.dot(symb, xc, yc, frame=frame, ctype=ctype, size=size)

    ds9.cmdBuffer.popSize()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#
# PSF display utilities
#
def showPsfSpatialCells(exposure, psfCellSet, nMaxPerCell=-1, showChi2=False, showMoments=False,
                        symb=None, ctype=None, ctypeUnused=None, ctypeBad=None, size=2, frame=None):
    """Show the SpatialCells.  If symb is something that ds9.dot understands (e.g. "o"), the top nMaxPerCell candidates will be indicated with that symbol, using ctype and size"""

    ds9.cmdBuffer.pushSize()

    origin = [-exposure.getMaskedImage().getX0(), -exposure.getMaskedImage().getY0()]
    for cell in psfCellSet.getCellList():
        displayUtils.drawBBox(cell.getBBox(), origin=origin, frame=frame)

        if nMaxPerCell < 0:
            nMaxPerCell = 0

        i = 0
        goodies = ctypeBad is None
        for cand in cell.begin(goodies):
            if nMaxPerCell > 0:
                i += 1

            cand = algorithmsLib.cast_PsfCandidateF(cand)

            xc, yc = cand.getXCenter() + origin[0], cand.getYCenter() + origin[1]

            if i > nMaxPerCell:
                if not ctypeUnused:
                    continue

            color = ctypeBad if cand.isBad() else ctype

            if symb:
                if i > nMaxPerCell:
                    ct = ctypeUnused
                else:
                    ct = ctype

                ds9.dot(symb, xc, yc, frame=frame, ctype=ct, size=size)

            source = cand.getSource()

            if showChi2:
                rchi2 = cand.getChi2()
                if rchi2 > 1e100:
                    rchi2 = numpy.nan
                ds9.dot("%d %.1f" % (source.getId(), rchi2),
                        xc - size, yc - size - 4, frame=frame, ctype=color, size=size)

            if showMoments:
                ds9.dot("%.2f %.2f %.2f" % (source.getIxx(), source.getIxy(), source.getIyy()),
                        xc-size, yc + size + 4, frame=frame, ctype=color, size=size)

    ds9.cmdBuffer.popSize()

def showPsfCandidates(exposure, psfCellSet, psf=None, frame=None, normalize=True, showBadCandidates=True,
                      variance=None, chi=None):
    """Display the PSF candidates.
If psf is provided include PSF model and residuals;  if normalize is true normalize the PSFs (and residuals)
If chi is True, generate a plot of residuals/sqrt(variance), i.e. chi
"""
    if chi is None:
        if variance is not None:        # old name for chi
            chi = variance
    #
    # Show us the ccandidates
    #
    mos = displayUtils.Mosaic()
    #
    candidateCenters = []
    candidateCentersBad = []
    candidateIndex = 0
    for cell in psfCellSet.getCellList():
        for cand in cell.begin(False): # include bad candidates
            cand = algorithmsLib.cast_PsfCandidateF(cand)

            rchi2 = cand.getChi2()
            if rchi2 > 1e100:
                rchi2 = numpy.nan

            if not showBadCandidates and cand.isBad():
                continue

            if psf:
                im_resid = displayUtils.Mosaic(gutter=0, background=-5, mode="x")

                try:
                    im = cand.getImage()
                    im = type(im)(im, True)
                    im.setXY0(cand.getImage().getXY0())
                except:
                    continue

                if not variance:
                    im_resid.append(type(im)(im, True))

                chi2 = algorithmsLib.subtractPsf(psf, im, cand.getXCenter(), cand.getYCenter())
                
                resid = im
                if variance:
                    resid = resid.getImage()
                    var = im.getVariance()
                    var = type(var)(var, True)
                    numpy.sqrt(var.getArray(), var.getArray()) # inplace sqrt
                    resid /= var
                    
                im_resid.append(resid)

                # Fit the PSF components directly to the data (i.e. ignoring the spatial model)
                im = cand.getImage()

                im = type(im)(im, True)
                im.setXY0(cand.getImage().getXY0())

                noSpatialKernel = afwMath.cast_LinearCombinationKernel(psf.getKernel())
                candCenter = afwGeom.PointD(cand.getXCenter(), cand.getYCenter())
                fit = algorithmsLib.fitKernelParamsToImage(noSpatialKernel, im, candCenter)
                params = fit[0]
                kernels = afwMath.KernelList(fit[1])
                outputKernel = afwMath.LinearCombinationKernel(kernels, params)

                outImage = afwImage.ImageD(outputKernel.getDimensions())
                outputKernel.computeImage(outImage, False)

                im -= outImage.convertF()
                resid = im
                if variance:
                    resid = resid.getImage()
                    resid /= var
                im_resid.append(resid)

                im = im_resid.makeMosaic()
            else:
                im = cand.getImage()

            if normalize:
                im /= afwMath.makeStatistics(im, afwMath.MAX).getValue()

            if psf:
                lab = "%d chi^2 %.1f" % (cand.getSource().getId(), rchi2)
                ctype = ds9.RED if cand.isBad() else ds9.GREEN
            else:
                lab = "%d flux %8.3g" % (cand.getSource().getId(), cand.getSource().getPsfFlux())
                ctype = ds9.GREEN

            mos.append(im, lab, ctype)

            if False and numpy.isnan(rchi2):
                ds9.mtv(cand.getImage().getImage(), title="candidate", frame=1)
                print "amp",  cand.getAmplitude()

            im = cand.getImage()
            center = (candidateIndex, cand.getXCenter() - im.getX0(), cand.getYCenter() - im.getY0())
            candidateIndex += 1
            if cand.isBad():
                candidateCentersBad.append(center)
            else:
                candidateCenters.append(center)

    if variance:
        title = "chi(Psf fit)"
    else:
        title = "Psf Candidates"
    mosaicImage = mos.makeMosaic(frame=frame, title=title)

    ds9.cmdBuffer.pushSize()

    for centers, color in ((candidateCenters, ds9.GREEN), (candidateCentersBad, ds9.RED)):
        for cen in centers:
            bbox = mos.getBBox(cen[0])
            ds9.dot("+", cen[1] + bbox.getMinX(), cen[2] + bbox.getMinY(), frame=frame, ctype=color)

    ds9.cmdBuffer.popSize()

    return mosaicImage

def plotPsfSpatialModel(exposure, psf, psfCellSet, showBadCandidates=True, numSample=128,
                        matchKernelAmplitudes=False, keepPlots=True):
    """Plot the PSF spatial model."""

    try:
        import numpy
        import matplotlib.pyplot as plt
        import matplotlib.colors
    except ImportError, e:
        print "Unable to import numpy and matplotlib: %s" % e
        return
    
    noSpatialKernel = afwMath.cast_LinearCombinationKernel(psf.getKernel())
    candPos = list()
    candFits = list()
    badPos = list()
    badFits = list()
    candAmps = list()
    badAmps = list()
    for cell in psfCellSet.getCellList():
        for cand in cell.begin(False):
            cand = algorithmsLib.cast_PsfCandidateF(cand)
            if not showBadCandidates and cand.isBad():
                continue
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

            targetFits = badFits if cand.isBad() else candFits
            targetPos = badPos if cand.isBad() else candPos
            targetAmps = badAmps if cand.isBad() else candAmps

            targetFits.append([x / amp for x in params])
            targetPos.append(candCenter)
            targetAmps.append(amp)

    numCandidates = len(candFits)
    numBasisFuncs = noSpatialKernel.getNBasisKernels()

    xGood = numpy.array([pos.getX() for pos in candPos])
    yGood = numpy.array([pos.getY() for pos in candPos])
    zGood = numpy.array(candFits)
    ampGood = numpy.array(candAmps)

    xBad = numpy.array([pos.getX() for pos in badPos])
    yBad = numpy.array([pos.getY() for pos in badPos])
    zBad = numpy.array(badFits)
    ampBad = numpy.array(badAmps)
    numBad = len(badPos)

    xRange = numpy.linspace(0, exposure.getWidth(), num=numSample)
    yRange = numpy.linspace(0, exposure.getHeight(), num=numSample)

    kernel = psf.getKernel()
    for k in range(kernel.getNKernelParameters()):
        func = kernel.getSpatialFunction(k)
        dfGood = zGood[:,k] - numpy.array([func(pos.getX(), pos.getY()) for pos in candPos])
        yMin = dfGood.min()
        yMax = dfGood.max()
        if numBad > 0:
            dfBad = zBad[:,k] - numpy.array([func(pos.getX(), pos.getY()) for pos in badPos])
            yMin = min([yMin, dfBad.min()])
            yMax = max([yMax, dfBad.max()])
        yMin -= 0.05 * (yMax - yMin)
        yMax += 0.05 * (yMax - yMin)

        fRange = numpy.ndarray((len(xRange), len(yRange)))
        for j, yVal in enumerate(yRange):
            for i, xVal in enumerate(xRange):
                fRange[j][i] = func(xVal, yVal)

        fig = plt.figure(k)

        fig.clf()
        try:
            fig.canvas._tkcanvas._root().lift() # == Tk's raise, but raise is a python reserved word
        except:                                 # protect against API changes
            pass

        fig.suptitle('Kernel component %d' % k)

        ax = fig.add_axes((0.1, 0.05, 0.35, 0.35))
        ax.set_autoscale_on(False)
        ax.set_xbound(lower=0, upper=exposure.getHeight())
        ax.set_ybound(lower=yMin, upper=yMax)
        ax.plot(yGood, dfGood, 'b+')
        if numBad > 0:
            ax.plot(yBad, dfBad, 'r+')
        ax.axhline(0.0)
        ax.set_title('Residuals as a function of y')
        ax = fig.add_axes((0.1, 0.55, 0.35, 0.35))
        ax.set_autoscale_on(False)
        ax.set_xbound(lower=0, upper=exposure.getWidth())
        ax.set_ybound(lower=yMin, upper=yMax)
        ax.plot(xGood, dfGood, 'b+')
        if numBad > 0:
            ax.plot(xBad, dfBad, 'r+')
        ax.axhline(0.0)
        ax.set_title('Residuals as a function of x')

        ax = fig.add_axes((0.55, 0.05, 0.35, 0.35))

        vmin = None
        vmax = None

        if matchKernelAmplitudes:
            if k == 0:
                vmin = 0.0
                vmax = fRange.max() # + 0.05 * numpy.fabs(fRange.max())
            else:
                pass
        else:
            vmin = fRange.min() # - 0.05 * numpy.fabs(fRange.min())
            vmax = fRange.max() # + 0.05 * numpy.fabs(fRange.max())

        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        im = ax.imshow(fRange, aspect='auto', norm=norm,
                       extent=[0, exposure.getWidth()-1, 0, exposure.getHeight()-1])
        ax.set_title('Spatial polynomial')
        plt.colorbar(im, orientation='horizontal', ticks=[vmin, vmax])

        ax = fig.add_axes((0.55, 0.55, 0.35, 0.35))
        if False:
            ax.scatter(xGood, yGood, c=dfGood, marker='o')
            ax.scatter(xBad, yBad, c=dfBad, marker='x')
            ax.set_xbound(lower=0, upper=exposure.getWidth())
            ax.set_ybound(lower=0, upper=exposure.getHeight())
            ax.set_title('Spatial residuals')
            plt.colorbar(im, orientation='horizontal')
        else:
            ax.plot(-2.5*numpy.log10(candAmps), zGood[:,k], 'b+')
            if numBad > 0:
                ax.plot(-2.5*numpy.log10(badAmps), zBad[:,k], 'r+')
            ax.set_ybound(lower=-1.0, upper=1.0)
            ax.set_title('Flux variation')

        fig.show()

    global keptPlots
    if keepPlots and not keptPlots:
        # Keep plots open when done
        def show():
            print "%s: Please close plots when done." % __name__
            try:
                plt.show()
            except:
                pass
            print "Plots closed, exiting..."
        import atexit
        atexit.register(show)
        keptPlots = True

def showPsf(psf, eigenValues=None, XY=None, normalize=True, frame=None):
    """Display a PSF's eigen images

    If normalize is True, set the largest absolute value of each eigenimage to 1.0 (n.b. sum == 0.0 for i > 0)
    """

    if eigenValues:
        coeffs = eigenValues
    elif XY is not None:
        coeffs = psf.getLocalKernel(afwGeom.PointD(XY[0], XY[1])).getKernelParameters()
    else:
        coeffs = None

    mos = displayUtils.Mosaic()
    i = 0
    for k in afwMath.cast_LinearCombinationKernel(psf.getKernel()).getKernelList():
        im = afwImage.ImageD(k.getDimensions())
        k.computeImage(im, False)
        if normalize:
            im /= numpy.max(numpy.abs(im.getArray()))
            
        if coeffs:
            mos.append(im, "%g" % (coeffs[i]/coeffs[0]))
            i += 1
        else:
            mos.append(im)

    mos.makeMosaic(frame=frame, title="Eigen Images")

    return mos

def showPsfMosaic(exposure, psf=None, nx=7, ny=None, frame=None):
    """Show a mosaic of Psf images.  exposure may be an Exposure (optionally with PSF), or a tuple (width, height)
    """
    mos = displayUtils.Mosaic()

    try:                                # maybe it's a real Exposure
        width, height = exposure.getWidth(), exposure.getHeight()
        if not psf:
            psf = exposure.getPsf()

            centroider = algorithmsLib.makeMeasureAstrometry(exposure)
            centroider.addAlgorithm("GAUSSIAN")
    except AttributeError:
        centroider = None
        try:                            # OK, maybe a list [width, height]
            width, height = exposure[0], exposure[1]
        except TypeError:               # I guess not
            raise RuntimeError, ("Unable to extract width/height from object of type %s" % type(exposure))

    if not ny:
        ny = int(nx*float(height)/width + 0.5)
        if not ny:
            ny = 1
    
    centroider = algorithmsLib.makeMeasureAstrometry(exposure)
    centroider.addAlgorithm(algorithmsLib.GaussianAstrometryControl())

    centers = []
    for iy in range(ny):
        for ix in range(nx):
            x = int(ix*(width-1)/(nx-1))
            y = int(iy*(height-1)/(ny-1))

            im = psf.computeImage(afwGeom.PointD(x, y)).convertF()
            mos.append(im, "PSF(%d,%d)" % (x, y))
    
            exp = afwImage.makeExposure(afwImage.makeMaskedImage(im))
            w, h = im.getWidth(), im.getHeight()
            cen = afwGeom.Point2D(im.getX0() + w//2, im.getY0() + h//2)
            src = afwDet.Source()
            foot = afwDet.Footprint(exp.getBBox())
            src.setFootprint(foot)

            c = centroider.measure(src, exp, cen).find()
                centers.append((c.getX() - im.getX0(), c.getY() - im.getY0()))

    mos.makeMosaic(frame=frame, title="Model Psf", mode=nx)

    if centers and frame is not None:
        ds9.cmdBuffer.pushSize()

        i = 0
        for cen in centers:
            bbox = mos.getBBox(i); i += 1
            ds9.dot("+", cen[0] + bbox.getMinX(), cen[1] + bbox.getMinY(), frame=frame)

        ds9.cmdBuffer.popSize()

    return mos

def showPsfResiduals(exposure, sourceSet, magType="psf", scale=10, frame=None, showAmps=False):
    mimIn = exposure.getMaskedImage()
    mimIn = mimIn.Factory(mimIn, True)  # make a copy to subtract from
    
    psf = exposure.getPsf()
    psfWidth, psfHeight = psf.getLocalKernel().getDimensions()
    #
    # Make the image that we'll paste our residuals into.  N.b. they can overlap the edges
    #
    w, h = int(mimIn.getWidth()/scale), int(mimIn.getHeight()/scale)

    im = mimIn.Factory(w + psfWidth, h + psfHeight)

    cenPos = []
    for s in sourceSet:
        x, y = s.getXAstrom(), s.getYAstrom()
        
        sx, sy = int(x/scale + 0.5), int(y/scale + 0.5)

        smim = im.Factory(im, afwGeom.BoxI(afwGeom.PointI(sx, sy), afwGeom.ExtentI(psfWidth, psfHeight)),
                         afwImage.PARENT)
        sim = smim.getImage()

        try:
            if magType == "ap":
                flux = s.getApFlux()
            elif magType == "model":
                flux = s.getModelFlux()
            elif magType == "psf":
                flux = s.getPsfFlux()
            else:
                raise RuntimeError("Unknown flux type %s" % magType)
            
            algorithmsLib.subtractPsf(psf, mimIn, x, y, flux)
        except Exception, e:
            print e

        try:
            expIm = mimIn.getImage().Factory(mimIn.getImage(),
                                             afwGeom.BoxI(afwGeom.PointI(int(x) - psfWidth//2,
                                                                         int(y) - psfHeight//2),
                                                          afwGeom.ExtentI(psfWidth, psfHeight)),
                                             afwImage.PARENT)
        except pexExcept.LsstCppException:
            continue

        cenPos.append([x - expIm.getX0() + sx, y - expIm.getY0() + sy])

        sim += expIm

    if frame is not None:
        ds9.mtv(im, frame=frame)
        ds9.cmdBuffer.pushSize()
        for x, y in cenPos:
            ds9.dot("+", x, y, frame=frame)
        ds9.cmdBuffer.popSize()

        if showAmps:
            nx, ny = namp
            for i in range(nx):
                for j in range(ny):
                    xc = numpy.array([0, 1, 1, 0, 0])
                    yc = numpy.array([0, 0, 1, 1, 0])

                    corners = []
                    for k in range(len(xc)):
                        corners.append([psfWidth//2 + w/nx*(i + xc[k]), psfHeight//2 + h/ny*(j + yc[k])])

                    ds9.line(corners, frame=frame)

    return im

def writeSourceSetAsCsv(sourceSet, fd=sys.stdout):
    """Write a SourceSet as a CSV file"""

    if not sourceSet:
        raise RuntimeError, "Please provide at least one Source"

    source = sourceSet[0]

    measurementTypes = (("astrometry", source.getAstrometry),
                        ("photometry", source.getPhotometry),
                        ("shape", source.getShape),
                        )

    print >> fd, "#misc::id:int:1"

    for measureType, getWhat in measurementTypes:
        a = getWhat()

        for value in a.getValues():
            for s in value.getSchema():
                if s.getType() == s.LONG:
                    typeName = "long"
                else:
                    typeName = "double"

                if s.isArray():
                    n = s.getDimen()
                else:
                    n = 1
                print >> fd, "#%s:%s:%s:%s:%d" % (measureType, value.getAlgorithm(), s.getName(),
                                                  typeName, n)


    for source in sourceSet:
        out = "%d" % (source.getId())

        for a in (source.getAstrometry(),
                  source.getPhotometry(),
                  source.getShape(),
                  ):
            for value in a.getValues():
                for sch in value.getSchema():
                    if out:
                        out += ", "
                    out += str(value.get(sch.getName()))

        print >> fd, out


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def saveSpatialCellSet(psfCellSet, fileName="foo.fits", frame=None):
    """Write the contents of a SpatialCellSet to a many-MEF fits file"""
    
    mode = "w"
    for cell in psfCellSet.getCellList():
        for cand in cell.begin(False):  # include bad candidates
            cand = algorithmsLib.cast_PsfCandidateF(cand)

            dx = afwImage.positionToIndex(cand.getXCenter(), True)[1]
            dy = afwImage.positionToIndex(cand.getYCenter(), True)[1]
            im = afwMath.offsetImage(cand.getImage(), -dx, -dy, "lanczos5")

            md = dafBase.PropertySet()
            md.set("CELL", cell.getLabel())
            md.set("ID", cand.getId())
            md.set("XCENTER", cand.getXCenter())
            md.set("YCENTER", cand.getYCenter())
            md.set("BAD", cand.isBad())
            md.set("AMPL", cand.getAmplitude())
            md.set("FLUX", cand.getSource().getPsfFlux())
            md.set("CHI2", cand.getSource().getChi2())

            im.writeFits(fileName, md, mode)
            mode = "a"

            if frame is not None:
                ds9.mtv(im, frame=frame)
