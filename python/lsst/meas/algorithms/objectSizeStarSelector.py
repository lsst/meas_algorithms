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
import collections
import math, sys

import numpy
try:
    import matplotlib.pyplot as pyplot
    fig = None
except ImportError:
    pyplot = None

import lsst.pex.config as pexConfig
import lsst.afw.detection as afwDetection
import lsst.afw.display.ds9 as ds9
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses as geomEllip
import lsst.afw.cameraGeom as cameraGeom
from . import algorithmsLib
from . import measurement
from lsst.meas.algorithms.starSelectorRegistry import starSelectorRegistry

class ObjectSizeStarSelectorConfig(pexConfig.Config):
    fluxMin = pexConfig.Field(
        doc = "specify the minimum psfFlux for good Psf Candidates",
        dtype = float,
        default = 12500.0,
#        minValue = 0.0,
        check = lambda x: x >= 0.0,
    )
    fluxMax = pexConfig.Field(
        doc = "specify the maximum psfFlux for good Psf Candidates (ignored if == 0)",
        dtype = float,
        default = 0.0,
        check = lambda x: x >= 0.0,
    )
    kernelSize = pexConfig.Field(
        doc = "size of the kernel to create",
        dtype = int,
        default = 21,
    )
    borderWidth = pexConfig.Field(
        doc = "number of pixels to ignore around the edge of PSF candidate postage stamps",
        dtype = int,
        default = 0,
    )
    badFlags = pexConfig.ListField(
        doc = "List of flags which cause a source to be rejected as bad",
        dtype = str,
        default = ["initial.flags.pixel.edge",
                   "initial.flags.pixel.interpolated.center",
                   "initial.flags.pixel.saturated.center",
                   "initial.flags.pixel.cr.center",
                   ]
        )
    histSize = pexConfig.Field(
        doc = "Number of bins in size histogram",
        dtype = int,
        default = 64,
        )
    widthMin = pexConfig.Field(
        doc = "minimum width to include in histogram",
        dtype = float,
        default = 0.0,
        check = lambda x: x >= 0.0,
        )
    widthMax = pexConfig.Field(
        doc = "maximum width to include in histogram",
        dtype = float,
        default = 10.0,
        check = lambda x: x >= 0.0,
        )

class EventHandler(object):
    """A class to handle key strokes with matplotlib displays"""
    def __init__(self, axes, xs, ys, x, y, frames=[0]):
        self.axes = axes
        self.xs = xs
        self.ys = ys
        self.x = x
        self.y = y
        self.frames = frames

        self.cid = self.axes.figure.canvas.mpl_connect('key_press_event', self)

    def __call__(self, ev):
        if ev.inaxes != self.axes:
            return
        
        if ev.key and ev.key in ("p"):
            dist = numpy.hypot(self.xs - ev.xdata, self.ys - ev.ydata)
            dist[numpy.where(numpy.isnan(dist))] = 1e30
            dmin = min(dist)

            which = numpy.where(dist == min(dist))

            x = self.x[which][0]
            y = self.y[which][0]
            for frame in self.frames:
                ds9.pan(x, y, frame=frame)
            ds9.cmdBuffer.flush()
        else:
            pass

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def assignClusters(yvec, centers):
    """Return a vector of centerIds based on their distance to the centers"""
    minDist = numpy.nan*numpy.ones_like(yvec)
    clusterId = numpy.empty_like(yvec, dtype=int)

    for i, mean in enumerate(centers):
        dist = abs(yvec - mean)
        if i == 0:
            update = dist == dist       # True for all points
        else:
            update = dist < minDist

        minDist[update] = dist[update]
        clusterId[update] = i

    return clusterId

def kcenters(yvec, nCluster=3,  useMedian=False):
    """A classic k-means algorithm, clustering yvec into nCluster clusters

    Return the set of centres, and the cluster ID for each of the points

    If useMedian is true, use the median of the cluster as its centre, rather than
    the traditional mean
    """
    mean0 = sorted(yvec)[len(yvec)//10] # guess
    centers = mean0*numpy.arange(1, nCluster + 1)
        
    func = numpy.median if useMedian else numpy.mean

    clusterId = numpy.zeros_like(yvec, dtype=int) - 1 # which cluster the points are assigned to
    while True:
        oclusterId = clusterId
        clusterId = assignClusters(yvec, centers)

        if numpy.all(clusterId == oclusterId):
            break

        for i in range(nCluster):
            centers[i] = func(yvec[clusterId == i])

    return centers, clusterId

def improveCluster(yvec, centers, clusterId, nsigma=2.0, niter=10, clusterNum=0):
    """Improve our estimate of one of the clusters (clusterNum) by sigma-clipping around its median"""

    nMember = sum(clusterId == clusterNum)
    if nMember < 5:
        return clusterId
    for iter in range(niter):
        old_nMember = nMember
        
        inCluster0 = clusterId == clusterNum
        yv = yvec[inCluster0]
        
        centers[clusterNum] = numpy.median(yv)
        stdev = numpy.std(yv)

        syv = sorted(yv)
        stdev_iqr = 0.741*(syv[int(0.75*nMember)] - syv[int(0.25*nMember)])

        sd = stdev if stdev < stdev_iqr else stdev_iqr

        if False:
            print "sigma(iqr) = %.3f, sigma = %.3f" % (stdev_iqr, numpy.std(yv))
        newCluster0 = abs(yvec - centers[clusterNum]) < nsigma*sd
        clusterId[numpy.logical_and(inCluster0, newCluster0)] = clusterNum
        clusterId[numpy.logical_and(inCluster0, numpy.logical_not(newCluster0))] = -1
        
        nMember = sum(clusterId == clusterNum)
        if nMember == old_nMember:
            break

    return clusterId

def plot(mag, width, centers, clusterId, marker="o", markersize=2, markeredgewidth=0, ltype='-',
         clear=True):

    global fig
    if not fig:
        fig = pyplot.figure()
        newFig = True
    else:
        newFig = False
        if clear:
            fig.clf()

    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80));

    xmin = sorted(mag)[int(0.05*len(mag))]
    xmax = sorted(mag)[int(0.95*len(mag))]

    axes.set_xlim(-17.5, -13)
    axes.set_xlim(xmin - 0.1*(xmax - xmin), xmax + 0.1*(xmax - xmin))
    axes.set_ylim(0, 10)

    colors = ["r", "g", "b", "c", "m", "k",]
    for k, mean in enumerate(centers):
        if k == 0:
            axes.plot(axes.get_xlim(), (mean, mean,), "k%s" % ltype)

        l = (clusterId == k)
        axes.plot(mag[l], width[l], marker, markersize=markersize, markeredgewidth=markeredgewidth,
                  color=colors[k%len(colors)])

    l = (clusterId == -1)
    axes.plot(mag[l], width[l], marker, markersize=markersize, markeredgewidth=markeredgewidth,
              color='k')

    if newFig:
        axes.set_xlabel("model")
        axes.set_ylabel(r"$\sqrt{I_{xx} + I_{yy}}$")

    return fig
        
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class ObjectSizeStarSelector(object):
    ConfigClass = ObjectSizeStarSelectorConfig

    def __init__(self, config, schema=None, key=None):
        """Construct a star selector that uses second moments
        
        This is a naive algorithm and should be used with caution.
        
        @param[in] config: An instance of ObjectSizeStarSelectorConfig
        @param[in,out] schema: An afw.table.Schema to register the selector's flag field.
                               If None, the sources will not be modified.
        @param[in] key: An existing Flag Key to use instead of registering a new field.
        """
        self._kernelSize  = config.kernelSize
        self._borderWidth = config.borderWidth
        self._widthMin = config.widthMin
        self._widthMax = config.widthMax
        self._fluxMin  = config.fluxMin
        self._fluxMax  = config.fluxMax
        self._badFlags = config.badFlags
        self._histSize = config.histSize
        if key is not None:
            self._key = key
            if schema is not None and key not in schema:
                raise LookupError("The key passed to the star selector is not present in the schema")
        elif schema is not None:
            self._key = schema.addField("classification.objectSize.star", type="Flag",
                                        doc="selected as a star by ObjectSizeStarSelector")
        else:
            self._key = None
            
    def selectStars(self, exposure, catalog, matches=None):
        """Return a list of PSF candidates that represent likely stars
        
        A list of PSF candidates may be used by a PSF fitter to construct a PSF.
        
        @param[in] exposure: the exposure containing the sources
        @param[in] catalog: a SourceCatalog containing sources that may be stars
        @param[in] matches: astrometric matches; ignored by this star selector
        
        @return psfCandidateList: a list of PSF candidates.
        """
        import lsstDebug
        display = lsstDebug.Info(__name__).display
        displayExposure = lsstDebug.Info(__name__).displayExposure     # display the Exposure + spatialCells
        plotMagSize = lsstDebug.Info(__name__).plotMagSize             # display the magnitude-size relation
        dumpData = lsstDebug.Info(__name__).dumpData                   # dump data to pickle file?

	detector = exposure.getDetector()
	distorter = None
	xy0 = afwGeom.Point2D(0,0)
	if not detector is None:
	    cPix = detector.getCenterPixel()
	    detSize = detector.getSize()
	    xy0.setX(cPix.getX() - int(0.5*detSize.getMm()[0]))
	    xy0.setY(cPix.getY() - int(0.5*detSize.getMm()[1]))
	    distorter = detector.getDistortion()
        #
        # Look at the distribution of stars in the magnitude-size plane
        #
        flux = catalog.get("initial.flux.gaussian")
        #mag = -2.5*numpy.log10(flux)

        xx = numpy.empty(len(catalog))
        xy = numpy.empty_like(xx)
        yy = numpy.empty_like(xx)
        for i, source in enumerate(catalog):
            Ixx, Ixy, Iyy = source.getIxx(), source.getIxy(), source.getIyy()
            if distorter:
                xpix, ypix = source.getX() + xy0.getX(), source.getY() + xy0.getY()
                p = afwGeom.Point2D(xpix, ypix)
                m = distorter.undistort(p, geomEllip.Quadrupole(Ixx, Iyy, Ixy), detector)
                Ixx, Ixy, Iyy = m.getIxx(), m.getIxy(), m.getIyy()

            xx[i], xy[i], yy[i] = Ixx, Ixy, Iyy
            
        width = numpy.sqrt(xx + yy)

        bad = reduce(lambda x, y: numpy.logical_or(x, catalog.get(y)), self._badFlags, False)
        bad = numpy.logical_or(bad, flux < self._fluxMin)
        bad = numpy.logical_or(bad, numpy.logical_not(numpy.isfinite(width)))
        bad = numpy.logical_or(bad, numpy.logical_not(numpy.isfinite(flux)))
        bad = numpy.logical_or(bad, width < self._widthMin)
        bad = numpy.logical_or(bad, width > self._widthMax)
        if self._fluxMax > 0:
            bad = numpy.logical_or(bad, flux > self._fluxMax)
        good = numpy.logical_not(bad)

        #mag = mag[good]
        mag = -2.5*numpy.log10(flux[good])
        width = width[good]
        #
        # Look for the maximum in the size histogram, then search upwards for the minimum that separates
        # the initial peak (of, we presume, stars) from the galaxies
        #
        if dumpData:
            import os, cPickle as pickle
            _ii = 0
            while True:
                pickleFile = os.path.expanduser(os.path.join("~", "widths-%d.pkl" % _ii))
                if not os.path.exists(pickleFile):
                    break
                _ii += 1

            with open(pickleFile, "wb") as fd:
                pickle.dump(mag, fd, -1)
                pickle.dump(width, fd, -1)

        centers, clusterId = kcenters(width, nCluster=4, useMedian=True)

        if display and plotMagSize and pyplot:
            fig = plot(mag, width, centers, clusterId,
                       marker="+", markersize=3, markeredgewidth=None, ltype=':', clear=True)
        else:
            fig = None
        
        clusterId = improveCluster(width, centers, clusterId)
        
        if display and plotMagSize and pyplot:
            plot(mag, width, centers, clusterId, marker="x", markersize=3, markeredgewidth=None)
        
        stellar = (clusterId == 0)
        #
        # We know enough to plot, if so requested
        #
        frame = 0

        if fig:
            if display and displayExposure:
                ds9.mtv(exposure.getMaskedImage(), frame=frame, title="PSF candidates")

                global eventHandler
                eventHandler = EventHandler(fig.get_axes()[0], mag, width,
                                            catalog.getX()[good], catalog.getY()[good], frames=[frame])

            fig.show()

            #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

            while True:
                try:
                    reply = raw_input("continue? [c h(elp) q(uit) p(db)] ").strip()
                except EOFError:
                    reply = "y"

                if reply:
                    if reply[0] == "h":
                        print """\
    We cluster the points; red are the stellar candidates and the other colours are other clusters.
    Points labelled + are rejects from the cluster (only for cluster 0).

    At this prompt, you can continue with almost any key; 'p' enters pdb, and 'h' prints this text

    If displayExposure is true, you can put the cursor on a point and hit 'p' to see it in ds9.
    """
                    elif reply[0] == "p":
                        import pdb; pdb.set_trace()
                    elif reply[0] == 'q':
                        sys.exit(1)
                    else:
                        break
        
        if display and displayExposure:
            mi = exposure.getMaskedImage()
    
            with ds9.Buffering():
                for i, source in enumerate(catalog):
                    if good[i]:
                        ctype = ds9.GREEN # star candidate
                    else:
                        ctype = ds9.RED # not star
			
                    ds9.dot("+", source.getX() - mi.getX0(),
                            source.getY() - mi.getY0(), frame=frame, ctype=ctype)
        #
        # Time to use that stellar classification to generate psfCandidateList
        #
        with ds9.Buffering():
            psfCandidateList = []
            for isStellar, source in zip(stellar, [s for g, s in zip(good, catalog) if g]):
                if not isStellar:
                    continue
                
                try:
                    psfCandidate = algorithmsLib.makePsfCandidate(source, exposure)
                    # The setXXX methods are class static, but it's convenient to call them on
                    # an instance as we don't know Exposure's pixel type
                    # (and hence psfCandidate's exact type)
                    if psfCandidate.getWidth() == 0:
                        psfCandidate.setBorderWidth(self._borderWidth)
                        psfCandidate.setWidth(self._kernelSize + 2*self._borderWidth)
                        psfCandidate.setHeight(self._kernelSize + 2*self._borderWidth)

                    im = psfCandidate.getMaskedImage().getImage()
                    vmax = afwMath.makeStatistics(im, afwMath.MAX).getValue()
                    if not numpy.isfinite(vmax):
                        continue
                    if self._key is not None:
                        source.set(self._key, True)
                    psfCandidateList.append(psfCandidate)

                    if display and displayExposure:
                        ds9.dot("o", source.getX() - mi.getX0(), source.getY() - mi.getY0(),
                                size=4, frame=frame, ctype=ds9.CYAN)
                except Exception as err:
                    pass # FIXME: should log this!

        return psfCandidateList

starSelectorRegistry.register("objectSize", ObjectSizeStarSelector)
