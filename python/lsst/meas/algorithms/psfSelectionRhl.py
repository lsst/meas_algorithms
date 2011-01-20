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
    
class PsfShapeHistogram(object):
    """A class to represent a histogram of (Ixx, Iyy)"""

    def __init__(self, xSize=40, ySize=40, xMax=30, yMax=30):
        """[xy]Size is the size of the psfImage; [xy]Max are the maximum values for I[xy][xy]"""

        self._xSize, self._ySize = xSize, ySize 
        self._xMax, self._yMax = xMax, yMax
        self._psfImage = afwImage.ImageF(self._xSize, self._ySize)
        self._psfImage.set(0)
        self._num = 0

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
                self._num += 1

    def peakToIxx(self, peakX, peakY):
        """Given a peak position in self._psfImage, return the corresponding (Ixx, Iyy)"""

        xx = peakX*self._xMax/self._xSize
        yy = peakY*self._yMax/self._ySize
        return xx, yy

    def getClump(self, display=False):
        if self._num <= 0:
            raise RuntimeError("No candidate PSF sources")
        
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
        psfImagePolicy = policy.Policy(policy.PolicyString(
            """#<?cfg paf policy?>
            source: {
                astrom: SDSS
                psfFlux: PSF
                apFlux: NAIVE
                shape: SDSS
            }
            astrometry: {
                SDSS: {
                    enabled: true
                }
            }
            photometry: {
                PSF: {
                    enabled: true
                }
                NAIVE: {
                    radius: 3.0
                }
            }
            shape: {
                SDSS: {
                    enabled: true
                }
            }
            """))
        
        sigma = 1
        exposure.setPsf(afwDetection.createPsf("DoubleGaussian", 11, 11, sigma))
        measureSources = algorithms.makeMeasureSources(exposure, psfImagePolicy)
        
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
                print "Except:", e
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



def selectPsfSources(exposure, sourceList, psfPolicy):
    """Get a list of suitable stars to construct a PSF."""

    try:
        import lsstDebug
        display = lsstDebug.Info(__name__).display
    except ImportError, e:
        try:
            type(display)
        except NameError:
            display = False
    
    #
    # Unpack policy
    #
    kernelSize   = psfPolicy.get("kernelSize")
    borderWidth  = psfPolicy.get("borderWidth")
    sizePsfCellX = psfPolicy.get("sizeCellX")
    sizePsfCellY = psfPolicy.get("sizeCellY")
    clumpNSigma  = psfPolicy.get("clumpNSigma")
    #
    # OK, we have all the source.  Let's do something with them
    #
    def goodPsfCandidate(source, fluxLim=psfPolicy.get("fluxLim")):
        """Should this object be included in the Ixx v. Iyy image?""" 

        badFlags = algorithms.Flags.EDGE | \
                   algorithms.Flags.INTERP_CENTER | \
                   algorithms.Flags.SATUR_CENTER | \
                   algorithms.Flags.PEAKCENTER

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
            ds9.dot("o", source.getXAstrom() - mi.getX0(),
                    source.getYAstrom() - mi.getY0(), frame=frame, ctype=ctype)

    psfClumpX, psfClumpY, psfClumpIxx, psfClumpIxy, psfClumpIyy = psfHist.getClump(display=display)
    #
    # Go through and find all the PSF-like objects
    #
    # We'll split the image into a number of cells, each of which contributes only
    # one PSF candidate star
    #
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
    for source in sourceList:
        Ixx, Ixy, Iyy = source.getIxx(), source.getIxy(), source.getIyy()
        dx, dy = (Ixx - psfClumpX), (Iyy - psfClumpY)

        if math.sqrt(a*dx*dx + 2*b*dx*dy + c*dy*dy) < 2*clumpNSigma: # A test for > would be confused by NaN
            if not goodPsfCandidate(source):
                continue

            try:
                cand = algorithms.makePsfCandidate(source, mi)
                #
                # The setXXX methods are class static, but it's convenient to call them on
                # an instance as we don't know Exposure's pixel type (and hence cand's exact type)
                if cand.getWidth() == 0:
                    cand.setBorderWidth(borderWidth)
                    cand.setWidth(kernelSize + 2*borderWidth)
                    cand.setHeight(kernelSize + 2*borderWidth)

                im = cand.getImage().getImage()
                max = afwMath.makeStatistics(im, afwMath.MAX).getValue()
                if not numpy.isfinite(max):
                    continue

                psfCellSet.insertCandidate(cand)

                if display:
                    ds9.dot("o", source.getXAstrom() - mi.getX0(), source.getYAstrom() - mi.getY0(),
                            size=4, frame=frame, ctype=ds9.CYAN)
            except Exception, e:
                continue

            source.setFlagForDetection(source.getFlagForDetection() | algorithms.Flags.STAR)
            psfStars += [source]

    return psfStars, psfCellSet

