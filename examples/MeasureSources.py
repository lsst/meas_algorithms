#!/usr/bin/env python
"""
Demonstrate running a simple image-processing pipeline

Run with:
   python MeasureSources.py
or
   python
   >>> import MeasureSources; MeasureSources.run()
"""

import pdb                              # we may want to say pdb.set_trace()
import glob, math, os, sys
from math import *
import eups
import lsst.daf.base as dafBase
import lsst.pex.logging as logging
import lsst.pex.policy as policy
import lsst.afw.detection as afwDetection
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.meas.algorithms as algorithms
import lsst.meas.algorithms.defects as defects
import lsst.meas.algorithms.measureSourceUtils as measureSourceUtils
import lsst.afw.display.ds9 as ds9

try:
    type(verbose)
except NameError:
    verbose = 0
logging.Trace_setVerbosity("meas.algorithms.measure", verbose)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class MO(object):
    """Measure the sources on a frame"""
    def __init__(self, display=False, rhs=None):

        self.display = display
        self.gas = None
        self.savedMask = None

        if rhs:
            try:
                self.exposure = rhs.exposure
                self.gas = rhs.gas
                self.pixscale = rhs.pixscale
                self.psf = rhs.psf
                self.psfImage = rhs.psfImage
                self.savedMask = rhs.savedMask
                self.sourceList = rhs.sourceList
                self.XY0 = rhs.XY0
            except AttributeError:
                raise RuntimeError, ("Unexpected rhs: %s" % rhs)

    def readData(self, fileName=None, subImage=False):
        if not fileName or isinstance(fileName, int):
            if fileName:
                which = fileName
            else:
                which = 1

            fileName = os.path.join(eups.productDir("afwdata"), "CFHT", "D4", "cal-53535-i-797722_%d" % which)
        #
        # We could read into an Exposure, but we're going to want to determine our own WCS
        #
        hdu, metadata = 0, dafBase.PropertySet()
        mi = afwImage.MaskedImageF(fileName, hdu, metadata) # read MaskedImage
        wcs = afwImage.Wcs(metadata)
        self.pixscale = 3600*math.sqrt(wcs.pixArea(wcs.getOriginRaDec()))
        #
        # Just an initial guess
        #
        FWHM = 5
        self.psf = algorithms.createPSF("DGPSF", 15, FWHM/(2*sqrt(2*log(2))))

        if subImage:                           # use sub-image
            self.XY0 = afwImage.PointI(824, 140)
            mi = mi.Factory(mi, afwImage.BBox(self.XY0, 512, 512))
        else:                       # use full image, trimmed to data section
            self.XY0 = afwImage.PointI(32, 2)
            mi = mi.Factory(mi, afwImage.BBox(self.XY0, afwImage.PointI(2079, 4609)))
            mi.setXY0(afwImage.PointI(0, 0))

        mi.getMask().addMaskPlane("DETECTED")
        self.exposure = afwImage.makeExposure(mi, wcs)

        if self.display:
            ds9.mtv(self.exposure)

    def ISR(self, fixCRs=True):
        """Run the ISR stage, removing CRs and patching bad columns"""
        mi = self.exposure.getMaskedImage()

        #
        # Fix defects
        #
        # Mask known bad pixels
        #
        badPixels = defects.policyToBadRegionList(os.path.join(eups.productDir("meas_algorithms"),
                                                               "pipeline/BadPixels.paf"))
        # did someone lie about the origin of the maskedImage?  If so, adjust bad pixel list
        if self.XY0.getX() != mi.getX0() or self.XY0.getY() != mi.getY0():
            dx = self.XY0.getX() - mi.getX0()
            dy = self.XY0.getY() - mi.getY0()
            for bp in badPixels:
                bp.shift(-dx, -dy)

        algorithms.interpolateOverDefects(mi, self.psf, badPixels)
        #
        # Subtract background
        #
        bctrl = afwMath.BackgroundControl(afwMath.NATURAL_SPLINE);
        bctrl.setNxSample(int(mi.getWidth()/256) + 1);
        bctrl.setNySample(int(mi.getHeight()/256) + 1);
	backobj = afwMath.makeBackground(mi.getImage(), bctrl)

        img = mi.getImage(); img -= backobj.getImageF(); del img
        #
        # Remove CRs
        #
        crPolicy = policy.Policy.createPolicy(os.path.join(eups.productDir("meas_algorithms"),
                                                           "pipeline", "CosmicRays.paf"))
        if fixCRs:
            crs = algorithms.findCosmicRays(mi, self.psf, 0, crPolicy)
        #
        # We do a pretty good job of interpolating, so don't propagagate the convolved CR/INTRP bits
        # (we'll keep them for the original CR/INTRP pixels)
        #
        self.savedMask = mi.getMask().Factory(mi.getMask(), True)
        saveBits = self.savedMask.getPlaneBitMask("CR") | \
                   self.savedMask.getPlaneBitMask("BAD") | \
                   self.savedMask.getPlaneBitMask("INTRP") # Bits to not convolve
        self.savedMask &= saveBits

        msk = mi.getMask(); msk &= ~saveBits; del msk # Clear the saved bits

    def measure(self):
        """Detect and measure sources"""
        mi = self.exposure.getMaskedImage()

        #
        # Smooth image
        #
        cnvImage = mi.Factory(mi.getDimensions())
        cnvImage.setXY0(afwImage.PointI(mi.getX0(), mi.getY0()))
        self.psf.convolve(cnvImage, mi, True, cnvImage.getMask().getMaskPlane("EDGE"))

        if self.savedMask:
            msk = cnvImage.getMask(); msk |= self.savedMask; del msk # restore the saved bits

        threshold = afwDetection.Threshold(3, afwDetection.Threshold.STDEV)
        #
        # Only search the part of the frame that was PSF-smoothed
        #        
        llc = afwImage.PointI(self.psf.getKernel().getWidth()/2, self.psf.getKernel().getHeight()/2)
        urc = afwImage.PointI(cnvImage.getWidth() - 1, cnvImage.getHeight() - 1) - llc;
        middle = cnvImage.Factory(cnvImage, afwImage.BBox(llc, urc))
        ds = afwDetection.DetectionSetF(middle, threshold, "DETECTED")
        del middle
        #
        # Reinstate the saved (e.g. BAD) (and also the DETECTED | EDGE) bits in the unsmoothed image
        #
        if self.savedMask:
            self.savedMask <<= cnvImage.getMask()
            msk = mi.getMask(); msk |= self.savedMask; del msk
            del self.savedMask; self.savedMask = None

        if self.display:
            ds9.mtv(mi, frame=0, lowOrderBits=True)
            ds9.mtv(cnvImage, frame=1)

        objects = ds.getFootprints()
        #
        # Time to actually measure
        #
        moPolicy = policy.Policy.createPolicy(os.path.join(eups.productDir("meas_algorithms"),
                                                           "pipeline", "MeasureSources.paf"))
        
        measureSources = algorithms.makeMeasureSources(self.exposure, moPolicy, self.psf)

        self.sourceList = afwDetection.SourceSet()
        for i in range(len(objects)):
            source = afwDetection.Source()
            self.sourceList.append(source)

            source.setId(i)
            source.setFlagForDetection(source.getFlagForDetection() | algorithms.Flags.BINNED1);

            try:
                measureSources.apply(source, objects[i])
            except Exception, e:
                print e
                continue

            if source.getFlagForDetection() & algorithms.Flags.EDGE:
                continue

            if self.display:
                xc, yc = source.getXAstrom() - mi.getX0(), source.getYAstrom() - mi.getY0()
                ds9.dot("+", xc, yc, size=1)
                
                if source.getFlagForDetection() & algorithms.Flags.INTERP_CENTER:
                    continue
                Mxx, Mxy, Myy = source.getFwhmA(), source.getFwhmTheta(), source.getFwhmB()
                ds9.dot("@:%g,%g,%g" % (Mxx, Mxy, Myy), xc, yc)
                
    def getPsfImage(self, fluxLim=1000):
        """Set the Mxx v. Myy image"""
        #
        # OK, we have all the source.  Let's do something with them
        #
        xSize, ySize = 20, 20
        xMax, yMax = 15, 15
        self.psfImage = afwImage.ImageF(xSize, ySize); self.psfImage.set(0)

        for source in self.sourceList:
            if fluxLim != None and source.getPsfMag() < fluxLim: # ignore faint objects
                continue
            #
            # Create an Image of Mxx v. Myy
            #
            i, j = int(source.getFwhmA()*xSize/xMax + 0.5), int(source.getFwhmB()*ySize/yMax + 0.5)
            if i in range(0, xSize) and j in range(0, ySize):
                if i == 0 and j == 0:
                    continue            # ignore the very smallest objects
                
                self.psfImage.set(i, j, self.psfImage.get(i, j) + 1)
        #
        # Embed self.psfImage into a larger image so we can smooth when measuring it
        #
        width, height = self.psfImage.getWidth(), self.psfImage.getHeight()
        
        limg = self.psfImage.Factory(2*width, 2*height)
        limg.set(0)

        sublimg = self.psfImage.Factory(limg, afwImage.BBox(afwImage.PointI(width, height), width, height))
        sublimg <<= self.psfImage
        del sublimg
        #
        # Now measure that image, looking for the highest peak.  Start by building an Exposure
        #
        msk = afwImage.MaskU(limg.getDimensions()); msk.set(0)
        var = afwImage.ImageF(limg.getDimensions()); var.set(1)
        mpsfImage = afwImage.MaskedImageF(limg, msk, var)
        mpsfImage.setXY0(afwImage.PointI(-width, -height))
        del msk; del var
        exposure = afwImage.makeExposure(mpsfImage)
        #
        # Next run an object detector
        #
        stats = afwMath.makeStatistics(self.psfImage, afwMath.STDEV)

        threshold = afwDetection.Threshold(2*stats.getValue(afwMath.STDEV))
        ds = afwDetection.DetectionSetF(mpsfImage, threshold, "DETECTED")
        objects = ds.getFootprints()
        #
        # And measure it
        #
        moPolicy = policy.Policy()
        moPolicy.add("measureObjects.centroidAlgorithm", "NAIVE")
        moPolicy.add("measureObjects.shapeAlgorithm", "SDSS")

        measureSources = algorithms.makeMeasureSources(exposure, moPolicy)

        if self.display:
            frame = 2
            ds9.mtv(mpsfImage.Factory(mpsfImage, afwImage.BBox(afwImage.PointI(width, height), width, height)),
                    frame=frame)

        sourceList = afwDetection.SourceSet()
        for i in range(len(objects)):
            source = afwDetection.Source()
            sourceList.append(source)

            source.setId(i)
            try:
                measureSources.apply(source, objects[i])
            except Exception, e:
                print e
                continue

            Mxx, Mxy, Myy = source.getFwhmA(), source.getFwhmTheta(), source.getFwhmB()

            if False:
                print source.getXAstrom(), source.getYAstrom(), \
                      Mxx, Mxy, Myy, \
                      "flags: ", measureSourceUtils.explainDetectionFlags(source.getFlagForDetection())
                
            if self.display:
                ds9.dot("+", source.getXAstrom(), source.getYAstrom(), size=0.5, ctype=ds9.RED, frame=frame)
                ds9.dot("@:%g,%g,%g" % (Mxx, Mxy, Myy), source.getXAstrom(), source.getYAstrom(), frame=frame)

    def write(self, filename, forFergal=False):
        if filename == "-":
            fd = sys.stdout
        else:
            fd = open(filename, "w")

        for source in self.sourceList:
            if forFergal:               # a format the Fergal used for the meas_astrom tests
                if source.getFlagForDetection() & (algorithms.Flags.EDGE | algorithms.Flags.SATUR_CENTER):
                    continue

                print >> fd, source.getXAstrom(), source.getYAstrom(), source.getPsfMag()
                continue

            print >> fd, "%-3d %7.2f %7.2f  %7.2f %7.2f   %7.3f %7.3f %7.3f   %8.1f" % \
                  (source.getId(),
                   source.getXAstrom(), source.getXAstromErr(),
                   source.getYAstrom(), source.getYAstromErr(),
                   source.getFwhmA(), source.getFwhmTheta(), source.getFwhmB(),
                   source.getPsfMag()),
            if fd == sys.stdout:
                print >> fd, measureSourceUtils.explainDetectionFlags(source.getFlagForDetection())
            else:
                print >> fd, ("0x%x" % source.getFlagForDetection())

    def read(self, filename, pixscale=0.18390):
        fd = open(filename, "r")

        self.pixscale = pixscale

        self.sourceList = afwDetection.SourceSet()
        for line in fd.readlines():
            source = afwDetection.Source()
            self.sourceList.append(source)

            vals = line.split()
            i = 0
            source.setId(int(vals[i])); i += 1
            source.setXAstrom(float(vals[i])); i += 1
            source.setXAstromErr(float(vals[i])); i += 1
            source.setYAstrom(float(vals[i])); i += 1
            source.setYAstromErr(float(vals[i])); i += 1
            source.setFwhmA(float(vals[i])); i += 1
            source.setFwhmTheta(float(vals[i])); i += 1
            source.setFwhmB(float(vals[i])); i += 1
            source.setPsfMag(float(vals[i])); i += 1
            source.setFlagForDetection(int(vals[i], 16)); i += 1
            
    def setWcs(self, fluxLim=None):
        import lsst.meas.astrom.net as astromNet

        if not self.gas:
            self.gas = astromNet.GlobalAstrometrySolution();

            #Read in the indices (i.e the files containing the positions of known asterisms
            #and add them to the gas object
            self.gas.setLogLevel(2)
            print >> sys.stderr, "Reading astrometry_net index files:"
            indices = glob.glob(os.path.join(eups.productDir("astrometry_net_data"), "index-2*.fits"))
            for f in indices:
                print >> sys.stderr, os.path.basename(f),
                self.gas.addIndexFile(f)
            print >> sys.stderr, ""
	
	# Set list of object positions
        starList = afwDetection.SourceSet()
        fluxLim = 10000
        for source in self.sourceList:
            if source.getFlagForDetection() & (algorithms.Flags.EDGE | algorithms.Flags.SATUR_CENTER):
                continue
            
            if fluxLim != None and source.getPsfMag() >= fluxLim: # ignore faint objects
                starList.append(source)

	self.gas.setStarlist(starList)
	self.gas.setNumberStars(len(starList))

        if False:
            self.gas.setImageScaleArcsecPerPixel(self.pixscale)
        else:
            self.gas.setMaximumImageScale(0.9*self.pixscale)
            self.gas.setMaximumImageScale(1.1*self.pixscale)

        if self.gas.blindSolve():
            self.exposure.setWcs(self.gas.getWcs())
        else:
            print "Failed to find WCS solution"

    def kitchenSink(self, subImage=False, fileName=None, fluxLim=3e5, fixCRs=True):
        """Do everything"""

        self.readData(fileName=fileName, subImage=subImage)
        self.ISR(fixCRs=fixCRs)
        self.measure()
        self.getPsfImage()
        if False:
            self.setWcs(fluxLim)

if __name__ == "__main__":
    if False:
        MeasureSources.MO().kitchenSink()
    else:
        try:
            mo = MO(display=1); mo.read("/u/rhl/LSST/meas/algorithms/foo.out"); mo.getPsfImage()
        except Exception, e:
            print e
