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

"""Support utilities for Measuring sources"""

__all__ = ["splitId", "showSourceSet", "showPsfSpatialCells",
           "showPsfCandidates", "makeSubplots", "plotPsfSpatialModel",
           "showPsf", "showPsfMosaic", "showPsfResiduals", "saveSpatialCellSet"]

import math
import numpy
import logging

import lsst.pex.exceptions as pexExcept
import lsst.daf.base as dafBase
import lsst.geom
import lsst.afw.geom as afwGeom
import lsst.afw.detection as afwDet
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.afw.display as afwDisplay
import lsst.afw.display.utils as displayUtils
import lsst.meas.base as measBase
from . import subtractPsf, fitKernelParamsToImage

keptPlots = False                       # Have we arranged to keep spatial plots open?

afwDisplay.setDefaultMaskTransparency(75)

_LOG = logging.getLogger(__name__)


def splitId(oid, asDict=True):

    objId = int((oid & 0xffff) - 1)      # Should be the value set by apps code

    if asDict:
        return dict(objId=objId)
    else:
        return [objId]


def showSourceSet(sSet, xy0=(0, 0), display=None, ctype=afwDisplay.GREEN, symb="+", size=2):
    """Draw the (XAstrom, YAstrom) positions of a set of Sources.  Image has the given XY0"""

    if not display:
        display = afwDisplay.Display()
    with display.Buffering():
        for s in sSet:
            xc, yc = s.getXAstrom() - xy0[0], s.getYAstrom() - xy0[1]

            if symb == "id":
                display.dot(str(splitId(s.getId(), True)["objId"]), xc, yc, ctype=ctype, size=size)
            else:
                display.dot(symb, xc, yc, ctype=ctype, size=size)

#
# PSF display utilities
#


def showPsfSpatialCells(exposure, psfCellSet, nMaxPerCell=-1, showChi2=False, showMoments=False,
                        symb=None, ctype=None, ctypeUnused=None, ctypeBad=None, size=2, display=None):
    """Show the SpatialCells.

    If symb is something that afwDisplay.Display.dot() understands (e.g. "o"),
    the top nMaxPerCell candidates will be indicated with that symbol, using
    ctype and size.
    """

    if not display:
        display = afwDisplay.Display()
    with display.Buffering():
        origin = [-exposure.getMaskedImage().getX0(), -exposure.getMaskedImage().getY0()]
        for cell in psfCellSet.getCellList():
            displayUtils.drawBBox(cell.getBBox(), origin=origin, display=display)

            if nMaxPerCell < 0:
                nMaxPerCell = 0

            i = 0
            goodies = ctypeBad is None
            for cand in cell.begin(goodies):
                if nMaxPerCell > 0:
                    i += 1

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

                    display.dot(symb, xc, yc, ctype=ct, size=size)

                source = cand.getSource()

                if showChi2:
                    rchi2 = cand.getChi2()
                    if rchi2 > 1e100:
                        rchi2 = numpy.nan
                    display.dot("%d %.1f" % (splitId(source.getId(), True)["objId"], rchi2),
                                xc - size, yc - size - 4, ctype=color, size=2)

                if showMoments:
                    display.dot("%.2f %.2f %.2f" % (source.getIxx(), source.getIxy(), source.getIyy()),
                                xc-size, yc + size + 4, ctype=color, size=size)
    return display


def showPsfCandidates(exposure, psfCellSet, psf=None, display=None, normalize=True, showBadCandidates=True,
                      fitBasisComponents=False, variance=None, chi=None):
    """Display the PSF candidates.

    If psf is provided include PSF model and residuals;  if normalize is true normalize the PSFs
    (and residuals)

    If chi is True, generate a plot of residuals/sqrt(variance), i.e. chi

    If fitBasisComponents is true, also find the best linear combination of the PSF's components
    (if they exist)
    """
    if not display:
        display = afwDisplay.Display()

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
        for cand in cell.begin(False):  # include bad candidates
            rchi2 = cand.getChi2()
            if rchi2 > 1e100:
                rchi2 = numpy.nan

            if not showBadCandidates and cand.isBad():
                continue

            if psf:
                im_resid = displayUtils.Mosaic(gutter=0, background=-5, mode="x")

                try:
                    im = cand.getMaskedImage()  # copy of this object's image
                    xc, yc = cand.getXCenter(), cand.getYCenter()

                    margin = 0 if True else 5
                    w, h = im.getDimensions()
                    bbox = lsst.geom.BoxI(lsst.geom.PointI(margin, margin), im.getDimensions())

                    if margin > 0:
                        bim = im.Factory(w + 2*margin, h + 2*margin)

                        stdev = numpy.sqrt(afwMath.makeStatistics(im.getVariance(), afwMath.MEAN).getValue())
                        afwMath.randomGaussianImage(bim.getImage(), afwMath.Random())
                        bim.getVariance().set(stdev**2)

                        bim.assign(im, bbox)
                        im = bim
                        xc += margin
                        yc += margin

                    im = im.Factory(im, True)
                    im.setXY0(cand.getMaskedImage().getXY0())
                except Exception:
                    continue

                if not variance:
                    im_resid.append(im.Factory(im, True))

                if True:                # tweak up centroids
                    mi = im
                    psfIm = mi.getImage()
                    config = measBase.SingleFrameMeasurementTask.ConfigClass()
                    config.slots.centroid = "base_SdssCentroid"

                    schema = afwTable.SourceTable.makeMinimalSchema()
                    measureSources = measBase.SingleFrameMeasurementTask(schema, config=config)
                    catalog = afwTable.SourceCatalog(schema)

                    extra = 10          # enough margin to run the sdss centroider
                    miBig = mi.Factory(im.getWidth() + 2*extra, im.getHeight() + 2*extra)
                    miBig[extra:-extra, extra:-extra, afwImage.LOCAL] = mi
                    miBig.setXY0(mi.getX0() - extra, mi.getY0() - extra)
                    mi = miBig
                    del miBig

                    exp = afwImage.makeExposure(mi)
                    exp.setPsf(psf)

                    footprintSet = afwDet.FootprintSet(mi,
                                                       afwDet.Threshold(0.5*numpy.max(psfIm.getArray())),
                                                       "DETECTED")
                    footprintSet.makeSources(catalog)

                    if len(catalog) == 0:
                        raise RuntimeError("Failed to detect any objects")

                    measureSources.run(catalog, exp)
                    if len(catalog) == 1:
                        source = catalog[0]
                    else:               # more than one source; find the once closest to (xc, yc)
                        dmin = None  # an invalid value to catch logic errors
                        for i, s in enumerate(catalog):
                            d = numpy.hypot(xc - s.getX(), yc - s.getY())
                            if i == 0 or d < dmin:
                                source, dmin = s, d
                    xc, yc = source.getCentroid()

                # residuals using spatial model
                try:
                    subtractPsf(psf, im, xc, yc)
                except Exception:
                    continue

                resid = im
                if variance:
                    resid = resid.getImage()
                    var = im.getVariance()
                    var = var.Factory(var, True)
                    numpy.sqrt(var.getArray(), var.getArray())  # inplace sqrt
                    resid /= var

                im_resid.append(resid)

                # Fit the PSF components directly to the data (i.e. ignoring the spatial model)
                if fitBasisComponents:
                    im = cand.getMaskedImage()

                    im = im.Factory(im, True)
                    im.setXY0(cand.getMaskedImage().getXY0())

                    try:
                        noSpatialKernel = psf.getKernel()
                    except Exception:
                        noSpatialKernel = None

                    if noSpatialKernel:
                        candCenter = lsst.geom.PointD(cand.getXCenter(), cand.getYCenter())
                        fit = fitKernelParamsToImage(noSpatialKernel, im, candCenter)
                        params = fit[0]
                        kernels = afwMath.KernelList(fit[1])
                        outputKernel = afwMath.LinearCombinationKernel(kernels, params)

                        outImage = afwImage.ImageD(outputKernel.getDimensions())
                        outputKernel.computeImage(outImage, False)

                        im -= outImage.convertF()
                        resid = im

                        if margin > 0:
                            bim = im.Factory(w + 2*margin, h + 2*margin)
                            afwMath.randomGaussianImage(bim.getImage(), afwMath.Random())
                            bim *= stdev

                            bim.assign(resid, bbox)
                            resid = bim

                        if variance:
                            resid = resid.getImage()
                            resid /= var

                        im_resid.append(resid)

                im = im_resid.makeMosaic()
            else:
                im = cand.getMaskedImage()

            if normalize:
                im /= afwMath.makeStatistics(im, afwMath.MAX).getValue()

            objId = splitId(cand.getSource().getId(), True)["objId"]
            if psf:
                lab = "%d chi^2 %.1f" % (objId, rchi2)
                ctype = afwDisplay.RED if cand.isBad() else afwDisplay.GREEN
            else:
                lab = "%d flux %8.3g" % (objId, cand.getSource().getPsfInstFlux())
                ctype = afwDisplay.GREEN

            mos.append(im, lab, ctype)

            if False and numpy.isnan(rchi2):
                display.mtv(cand.getMaskedImage().getImage(), title="showPsfCandidates: candidate")
                print("amp", cand.getAmplitude())

            im = cand.getMaskedImage()
            center = (candidateIndex, xc - im.getX0(), yc - im.getY0())
            candidateIndex += 1
            if cand.isBad():
                candidateCentersBad.append(center)
            else:
                candidateCenters.append(center)

    if variance:
        title = "chi(Psf fit)"
    else:
        title = "Stars & residuals"
    mosaicImage = mos.makeMosaic(display=display, title=title)

    with display.Buffering():
        for centers, color in ((candidateCenters, afwDisplay.GREEN), (candidateCentersBad, afwDisplay.RED)):
            for cen in centers:
                bbox = mos.getBBox(cen[0])
                display.dot("+", cen[1] + bbox.getMinX(), cen[2] + bbox.getMinY(), ctype=color)

    return mosaicImage


def makeSubplots(fig, nx=2, ny=2, Nx=1, Ny=1, plottingArea=(0.1, 0.1, 0.85, 0.80),
                 pxgutter=0.05, pygutter=0.05, xgutter=0.04, ygutter=0.04,
                 headroom=0.0, panelBorderWeight=0, panelColor='black'):
    """Return a generator of a set of subplots, a set of Nx*Ny panels of nx*ny plots.  Each panel is fully
    filled by row (starting in the bottom left) before the next panel is started.  If panelBorderWidth is
    greater than zero a border is drawn around each panel, adjusted to enclose the axis labels.

    E.g.
    subplots = makeSubplots(fig, 2, 2, Nx=1, Ny=1, panelColor='k')
    ax = subplots.next(); ax.text(0.3, 0.5, '[0, 0] (0,0)')
    ax = subplots.next(); ax.text(0.3, 0.5, '[0, 0] (1,0)')
    ax = subplots.next(); ax.text(0.3, 0.5, '[0, 0] (0,1)')
    ax = subplots.next(); ax.text(0.3, 0.5, '[0, 0] (1,1)')
    fig.show()

    Parameters
    ----------
    fig : `matplotlib.pyplot.figure`
        The matplotlib figure to draw
    nx : `int`
        The number of plots in each row of each panel
    ny : `int`
        The number of plots in each column of each panel
    Nx : `int`
        The number of panels in each row of the figure
    Ny : `int`
        The number of panels in each column of the figure
    plottingArea : `tuple`
        (x0, y0, x1, y1) for the part of the figure containing all the panels
    pxgutter : `float`
        Spacing between columns of panels in units of (x1 - x0)
    pygutter : `float`
        Spacing between rows of panels in units of (y1 - y0)
    xgutter : `float`
        Spacing between columns of plots within a panel in units of (x1 - x0)
    ygutter : `float`
        Spacing between rows of plots within a panel in units of (y1 - y0)
    headroom : `float`
        Extra spacing above each plot for e.g. a title
    panelBorderWeight : `int`
        Width of border drawn around panels
    panelColor : `str`
        Colour of border around panels
    """

    log = _LOG.getChild("makeSubplots")
    try:
        import matplotlib.pyplot as plt
    except ImportError as e:
        log.warning("Unable to import matplotlib: %s", e)
        return

    # Make show() call canvas.draw() too so that we know how large the axis labels are.  Sigh
    try:
        fig.__show
    except AttributeError:
        fig.__show = fig.show

        def myShow(fig):
            fig.__show()
            fig.canvas.draw()

        import types
        fig.show = types.MethodType(myShow, fig)
    #
    # We can't get the axis sizes until after draw()'s been called, so use a callback  Sigh^2
    #
    axes = {}                           # all axes in all the panels we're drawing: axes[panel][0] etc.
    #

    def on_draw(event):
        """
        Callback to draw the panel borders when the plots are drawn to the canvas
        """
        if panelBorderWeight <= 0:
            return False

        for p in axes.keys():
            bboxes = []
            for ax in axes[p]:
                bboxes.append(ax.bbox.union([label.get_window_extent() for label in
                                             ax.get_xticklabels() + ax.get_yticklabels()]))

            ax = axes[p][0]

            # this is the bbox that bounds all the bboxes, again in relative
            # figure coords

            bbox = ax.bbox.union(bboxes)

            xy0, xy1 = ax.transData.inverted().transform(bbox)
            x0, y0 = xy0
            x1, y1 = xy1
            w, h = x1 - x0, y1 - y0
            # allow a little space around BBox
            x0 -= 0.02*w
            w += 0.04*w
            y0 -= 0.02*h
            h += 0.04*h
            h += h*headroom
            # draw BBox
            ax.patches = []             # remove old ones
            rec = ax.add_patch(plt.Rectangle((x0, y0), w, h, fill=False,
                                             lw=panelBorderWeight, edgecolor=panelColor))
            rec.set_clip_on(False)

        return False

    fig.canvas.mpl_connect('draw_event', on_draw)
    #
    # Choose the plotting areas for each subplot
    #
    x0, y0 = plottingArea[0:2]
    W, H = plottingArea[2:4]
    w = (W - (Nx - 1)*pxgutter - (nx*Nx - 1)*xgutter)/float(nx*Nx)
    h = (H - (Ny - 1)*pygutter - (ny*Ny - 1)*ygutter)/float(ny*Ny)
    #
    # OK!  Time to create the subplots
    #
    for panel in range(Nx*Ny):
        axes[panel] = []
        px = panel%Nx
        py = panel//Nx
        for window in range(nx*ny):
            x = nx*px + window%nx
            y = ny*py + window//nx
            ax = fig.add_axes((x0 + xgutter + pxgutter + x*w + (px - 1)*pxgutter + (x - 1)*xgutter,
                               y0 + ygutter + pygutter + y*h + (py - 1)*pygutter + (y - 1)*ygutter,
                               w, h), frame_on=True, facecolor='w')
            axes[panel].append(ax)
            yield ax


def plotPsfSpatialModel(exposure, psf, psfCellSet, showBadCandidates=True, numSample=128,
                        matchKernelAmplitudes=False, keepPlots=True):
    """Plot the PSF spatial model."""

    log = _LOG.getChild("plotPsfSpatialModel")
    try:
        import matplotlib.pyplot as plt
        import matplotlib as mpl
    except ImportError as e:
        log.warning("Unable to import matplotlib: %s", e)
        return

    noSpatialKernel = psf.getKernel()
    candPos = list()
    candFits = list()
    badPos = list()
    badFits = list()
    candAmps = list()
    badAmps = list()
    for cell in psfCellSet.getCellList():
        for cand in cell.begin(False):
            if not showBadCandidates and cand.isBad():
                continue
            candCenter = lsst.geom.PointD(cand.getXCenter(), cand.getYCenter())
            try:
                im = cand.getMaskedImage()
            except Exception:
                continue

            fit = fitKernelParamsToImage(noSpatialKernel, im, candCenter)
            params = fit[0]
            kernels = fit[1]
            amp = 0.0
            for p, k in zip(params, kernels):
                amp += p * k.getSum()

            targetFits = badFits if cand.isBad() else candFits
            targetPos = badPos if cand.isBad() else candPos
            targetAmps = badAmps if cand.isBad() else candAmps

            targetFits.append([x / amp for x in params])
            targetPos.append(candCenter)
            targetAmps.append(amp)

    xGood = numpy.array([pos.getX() for pos in candPos]) - exposure.getX0()
    yGood = numpy.array([pos.getY() for pos in candPos]) - exposure.getY0()
    zGood = numpy.array(candFits)

    xBad = numpy.array([pos.getX() for pos in badPos]) - exposure.getX0()
    yBad = numpy.array([pos.getY() for pos in badPos]) - exposure.getY0()
    zBad = numpy.array(badFits)
    numBad = len(badPos)

    xRange = numpy.linspace(0, exposure.getWidth(), num=numSample)
    yRange = numpy.linspace(0, exposure.getHeight(), num=numSample)

    kernel = psf.getKernel()
    nKernelComponents = kernel.getNKernelParameters()
    #
    # Figure out how many panels we'll need
    #
    nPanelX = int(math.sqrt(nKernelComponents))
    nPanelY = nKernelComponents//nPanelX
    while nPanelY*nPanelX < nKernelComponents:
        nPanelX += 1

    fig = plt.figure(1)
    fig.clf()
    try:
        fig.canvas._tkcanvas._root().lift()  # == Tk's raise, but raise is a python reserved word
    except Exception:                                  # protect against API changes
        pass
    #
    # Generator for axes arranged in panels
    #
    mpl.rcParams["figure.titlesize"] = "x-small"
    subplots = makeSubplots(fig, 2, 2, Nx=nPanelX, Ny=nPanelY, xgutter=0.06, ygutter=0.06, pygutter=0.04)

    for k in range(nKernelComponents):
        func = kernel.getSpatialFunction(k)
        dfGood = zGood[:, k] - numpy.array([func(pos.getX(), pos.getY()) for pos in candPos])
        yMin = dfGood.min()
        yMax = dfGood.max()
        if numBad > 0:
            dfBad = zBad[:, k] - numpy.array([func(pos.getX(), pos.getY()) for pos in badPos])
            yMin = min([yMin, dfBad.min()])
            yMax = max([yMax, dfBad.max()])
        yMin -= 0.05 * (yMax - yMin)
        yMax += 0.05 * (yMax - yMin)

        yMin = -0.01
        yMax = 0.01

        fRange = numpy.ndarray((len(xRange), len(yRange)))
        for j, yVal in enumerate(yRange):
            for i, xVal in enumerate(xRange):
                fRange[j][i] = func(xVal, yVal)

        ax = next(subplots)

        ax.set_autoscale_on(False)
        ax.set_xbound(lower=0, upper=exposure.getHeight())
        ax.set_ybound(lower=yMin, upper=yMax)
        ax.plot(yGood, dfGood, 'b+')
        if numBad > 0:
            ax.plot(yBad, dfBad, 'r+')
        ax.axhline(0.0)
        ax.set_title('Residuals(y)')

        ax = next(subplots)

        if matchKernelAmplitudes and k == 0:
            vmin = 0.0
            vmax = 1.1
        else:
            vmin = fRange.min()
            vmax = fRange.max()

        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        im = ax.imshow(fRange, aspect='auto', origin="lower", norm=norm,
                       extent=[0, exposure.getWidth()-1, 0, exposure.getHeight()-1])
        ax.set_title('Spatial poly')
        plt.colorbar(im, orientation='horizontal', ticks=[vmin, vmax])

        ax = next(subplots)
        ax.set_autoscale_on(False)
        ax.set_xbound(lower=0, upper=exposure.getWidth())
        ax.set_ybound(lower=yMin, upper=yMax)
        ax.plot(xGood, dfGood, 'b+')
        if numBad > 0:
            ax.plot(xBad, dfBad, 'r+')
        ax.axhline(0.0)
        ax.set_title('K%d Residuals(x)' % k)

        ax = next(subplots)

        photoCalib = exposure.getPhotoCalib()
        # If there is no calibration factor, use 1.0.
        if photoCalib.getCalibrationMean() <= 0:
            photoCalib = afwImage.PhotoCalib(1.0)

        ampMag = [photoCalib.instFluxToMagnitude(candAmp) for candAmp in candAmps]
        ax.plot(ampMag, zGood[:, k], 'b+')
        if numBad > 0:
            badAmpMag = [photoCalib.instFluxToMagnitude(badAmp) for badAmp in badAmps]
            ax.plot(badAmpMag, zBad[:, k], 'r+')

        ax.set_title('Flux variation')

    fig.show()

    global keptPlots
    if keepPlots and not keptPlots:
        # Keep plots open when done
        def show():
            print("%s: Please close plots when done." % __name__)
            try:
                plt.show()
            except Exception:
                pass
            print("Plots closed, exiting...")
        import atexit
        atexit.register(show)
        keptPlots = True


def showPsf(psf, eigenValues=None, XY=None, normalize=True, display=None):
    """Display a PSF's eigen images

    If normalize is True, set the largest absolute value of each eigenimage to 1.0 (n.b. sum == 0.0 for i > 0)
    """

    if eigenValues:
        coeffs = eigenValues
    elif XY is not None:
        coeffs = psf.getLocalKernel(lsst.geom.PointD(XY[0], XY[1])).getKernelParameters()
    else:
        coeffs = None

    mos = displayUtils.Mosaic(gutter=2, background=-0.1)
    for i, k in enumerate(psf.getKernel().getKernelList()):
        im = afwImage.ImageD(k.getDimensions())
        k.computeImage(im, False)
        if normalize:
            im /= numpy.max(numpy.abs(im.getArray()))

        if coeffs:
            mos.append(im, "%g" % (coeffs[i]/coeffs[0]))
        else:
            mos.append(im)

    if not display:
        display = afwDisplay.Display()
    mos.makeMosaic(display=display, title="Kernel Basis Functions")

    return mos


def showPsfMosaic(exposure, psf=None, nx=7, ny=None, showCenter=True, showEllipticity=False,
                  showFwhm=False, stampSize=0, display=None, title=None):
    """Show a mosaic of Psf images.  exposure may be an Exposure (optionally with PSF),
    or a tuple (width, height)

    If stampSize is > 0, the psf images will be trimmed to stampSize*stampSize
    """

    scale = 1.0
    if showFwhm:
        showEllipticity = True
        scale = 2*math.log(2)         # convert sigma^2 to HWHM^2 for a Gaussian

    mos = displayUtils.Mosaic()

    try:                                # maybe it's a real Exposure
        width, height = exposure.getWidth(), exposure.getHeight()
        x0, y0 = exposure.getXY0()
        if not psf:
            psf = exposure.getPsf()
    except AttributeError:
        try:                            # OK, maybe a list [width, height]
            width, height = exposure[0], exposure[1]
            x0, y0 = 0, 0
        except TypeError:               # I guess not
            raise RuntimeError("Unable to extract width/height from object of type %s" % type(exposure))

    if not ny:
        ny = int(nx*float(height)/width + 0.5)
        if not ny:
            ny = 1

    centroidName = "SdssCentroid"
    shapeName = "base_SdssShape"

    schema = afwTable.SourceTable.makeMinimalSchema()
    schema.getAliasMap().set("slot_Centroid", centroidName)
    schema.getAliasMap().set("slot_Centroid_flag", centroidName+"_flag")

    control = measBase.SdssCentroidControl()
    centroider = measBase.SdssCentroidAlgorithm(control, centroidName, schema)

    sdssShape = measBase.SdssShapeControl()
    shaper = measBase.SdssShapeAlgorithm(sdssShape, shapeName, schema)
    table = afwTable.SourceTable.make(schema)

    table.defineCentroid(centroidName)
    table.defineShape(shapeName)

    bbox = None
    if stampSize > 0:
        w, h = psf.computeImage(lsst.geom.PointD(0, 0)).getDimensions()
        if stampSize <= w and stampSize <= h:
            bbox = lsst.geom.BoxI(lsst.geom.PointI((w - stampSize)//2, (h - stampSize)//2),
                                  lsst.geom.ExtentI(stampSize, stampSize))

    centers = []
    shapes = []
    for iy in range(ny):
        for ix in range(nx):
            x = int(ix*(width-1)/(nx-1)) + x0
            y = int(iy*(height-1)/(ny-1)) + y0

            im = psf.computeImage(lsst.geom.PointD(x, y)).convertF()
            imPeak = psf.computePeak(lsst.geom.PointD(x, y))
            im /= imPeak
            if bbox:
                im = im.Factory(im, bbox)
            lab = "PSF(%d,%d)" % (x, y) if False else ""
            mos.append(im, lab)

            exp = afwImage.makeExposure(afwImage.makeMaskedImage(im))
            exp.setPsf(psf)
            w, h = im.getWidth(), im.getHeight()
            centerX = im.getX0() + w//2
            centerY = im.getY0() + h//2
            src = table.makeRecord()
            spans = afwGeom.SpanSet(exp.getBBox())
            foot = afwDet.Footprint(spans)
            foot.addPeak(centerX, centerY, 1)
            src.setFootprint(foot)

            try:
                centroider.measure(src, exp)
                centers.append((src.getX() - im.getX0(), src.getY() - im.getY0()))

                shaper.measure(src, exp)
                shapes.append((src.getIxx(), src.getIxy(), src.getIyy()))
            except Exception:
                pass

    if not display:
        display = afwDisplay.Display()
    mos.makeMosaic(display=display, title=title if title else "Model Psf", mode=nx)

    if centers and display:
        with display.Buffering():
            for i, (cen, shape) in enumerate(zip(centers, shapes)):
                bbox = mos.getBBox(i)
                xc, yc = cen[0] + bbox.getMinX(), cen[1] + bbox.getMinY()
                if showCenter:
                    display.dot("+", xc, yc, ctype=afwDisplay.BLUE)

                if showEllipticity:
                    ixx, ixy, iyy = shape
                    ixx *= scale
                    ixy *= scale
                    iyy *= scale
                    display.dot("@:%g,%g,%g" % (ixx, ixy, iyy), xc, yc, ctype=afwDisplay.RED)

    return mos


def showPsfResiduals(exposure, sourceSet, magType="psf", scale=10, display=None):
    mimIn = exposure.getMaskedImage()
    mimIn = mimIn.Factory(mimIn, True)  # make a copy to subtract from

    psf = exposure.getPsf()
    psfWidth, psfHeight = psf.getLocalKernel(psf.getAveragePosition()).getDimensions()
    #
    # Make the image that we'll paste our residuals into.  N.b. they can overlap the edges
    #
    w, h = int(mimIn.getWidth()/scale), int(mimIn.getHeight()/scale)

    im = mimIn.Factory(w + psfWidth, h + psfHeight)

    cenPos = []
    for s in sourceSet:
        x, y = s.getX(), s.getY()

        sx, sy = int(x/scale + 0.5), int(y/scale + 0.5)

        smim = im.Factory(im, lsst.geom.BoxI(lsst.geom.PointI(sx, sy),
                                             lsst.geom.ExtentI(psfWidth, psfHeight)))
        sim = smim.getImage()

        try:
            if magType == "ap":
                flux = s.getApInstFlux()
            elif magType == "model":
                flux = s.getModelInstFlux()
            elif magType == "psf":
                flux = s.getPsfInstFlux()
            else:
                raise RuntimeError("Unknown flux type %s" % magType)

            subtractPsf(psf, mimIn, x, y, flux)
        except Exception as e:
            print(e)

        try:
            expIm = mimIn.getImage().Factory(mimIn.getImage(),
                                             lsst.geom.BoxI(lsst.geom.PointI(int(x) - psfWidth//2,
                                                                             int(y) - psfHeight//2),
                                                            lsst.geom.ExtentI(psfWidth, psfHeight)),
                                             )
        except pexExcept.Exception:
            continue

        cenPos.append([x - expIm.getX0() + sx, y - expIm.getY0() + sy])

        sim += expIm

    if display:
        display = afwDisplay.Display()
        display.mtv(im, title="showPsfResiduals: image")
        with display.Buffering():
            for x, y in cenPos:
                display.dot("+", x, y)

    return im


def saveSpatialCellSet(psfCellSet, fileName="foo.fits", display=None):
    """Write the contents of a SpatialCellSet to a many-MEF fits file"""

    mode = "w"
    for cell in psfCellSet.getCellList():
        for cand in cell.begin(False):  # include bad candidates
            dx = afwImage.positionToIndex(cand.getXCenter(), True)[1]
            dy = afwImage.positionToIndex(cand.getYCenter(), True)[1]
            im = afwMath.offsetImage(cand.getMaskedImage(), -dx, -dy, "lanczos5")

            md = dafBase.PropertySet()
            md.set("CELL", cell.getLabel())
            md.set("ID", cand.getId())
            md.set("XCENTER", cand.getXCenter())
            md.set("YCENTER", cand.getYCenter())
            md.set("BAD", cand.isBad())
            md.set("AMPL", cand.getAmplitude())
            md.set("FLUX", cand.getSource().getPsfInstFlux())
            md.set("CHI2", cand.getSource().getChi2())

            im.writeFits(fileName, md, mode)
            mode = "a"

            if display:
                display.mtv(im, title="saveSpatialCellSet: image")
