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

import math
import lsst.pex.policy        as pexPolicy
import lsst.pex.logging       as pexLog
import lsst.pex.policy        as policy
import lsst.meas.algorithms   as measAlg
import lsst.afw.detection     as afwDetection
import lsst.afw.geom          as afwGeom
import lsst.afw.coord         as afwCoord
import lsst.afw.display.ds9   as ds9


def sourceMeasurement(
    exposure,                 # exposure to analyse
    psf,                      # psf
    footprintLists,           # footprints of the detected objects
    measSourcePolicy,         # measureSources policy
    ):
    """ Source Measurement """

    try:
        import lsstDebug

        display = lsstDebug.Info(__name__).display
    except ImportError, e:
        try:
            display
        except NameError:
            display = False

    if display:
        frame = 0
        ds9.mtv(exposure, title="input", frame=frame)

    # instantiate a measurement object for 
    # - instantiation only involves looking up the algorithms for centroid, shape, and photometry
    #   ... nothing actually gets measured yet.
    exposure.setPsf(psf)
    measureSources = measAlg.makeMeasureSources(exposure, measSourcePolicy)

    # create an empty list to contain the sources we found (as Source objects)
    sourceSet = afwDetection.SourceSet()
    
    for footprintList in footprintLists:
        footprints, isNegative = footprintList

        # loop over all the objects detected
        for i in range(len(footprints)):

            # create a new source, and add it to the list, initialize ...
            source = afwDetection.Source()
            sourceSet.append(source)
            source.setId(i)

            source.setFlagForDetection(source.getFlagForDetection() | measAlg.Flags.BINNED1);

            # actually try to "measure" this object
            # recall: footprints[i] is a footprint for an object, measured values will be recorded in 'source'
            try:
                measureSources.apply(source, footprints[i])
            except Exception, e:
                # logging might be nice here
                #self.log.log(Log.WARN, str(e))
                pass
            #
            # Set the time
            #
            if False and exposure.getDetector():
                pos = afwGeom.Point2I(int(source.getXAstrom()), int(source.getYAstrom()))
                midTime = exposure.getCalib().getMidTime(exposure.getDetector(), pos)
            else:
                midTime = exposure.getCalib().getMidTime()
                
            source.setTaiMidPoint(midTime.get())
            source.setTaiRange(exposure.getCalib().getExptime())

            if display:
                ds9.dot("+", source.getXAstrom(), source.getYAstrom(), size=3, ctype=ds9.RED)
                if display > 1:
                    ds9.dot(("@:%.1f,%.1f,%1f" % (source.getXAstromErr()**2, 0, source.getYAstromErr()**2)),
                            source.getXAstrom(), source.getYAstrom(), size=3, ctype=ds9.RED)

    return sourceSet

def sourceMeasurementPsfFlux(
    exposure,                 # exposure to analyse
    sourceSet                 # the Sources to measure
    ):
    """Measure the psfFluxes for a SourceSet"""

    measSourcePolicy = pexPolicy.Policy(pexPolicy.PolicyString(
        """#<?cfg paf policy?>
        source: {
        psfFlux: "PSF"
        }
        photometry: {
        enabled: true
        }
        """
        ))
    measureSources = measAlg.makeMeasureSources(exposure, measSourcePolicy)

    for source in sourceSet:
        try:
            measureSources.apply(source)
        except Exception, e:
            # logging might be nice here
            #self.log.log(Log.WARN, str(e))
            pass

def computeSkyCoords(wcs, sourceSet):
    log = pexLog.Log(pexLog.getDefaultLog(), 'lsst.meas.utils.sourceMeasurement.computeSkyCoords')
    if sourceSet is None:
        log.log(Log.WARN, "No SourceSet provided" )
        return
    if wcs is None:
        log.log(Log.WARN, "No WCS provided")

    for s in sourceSet:
        (ra, dec, raErr, decErr) = xyToRaDec(
                s.getXFlux(), 
                s.getYFlux(),
                s.getXFluxErr(), 
                s.getYFluxErr(), 
                wcs)
        s.setRaFlux(ra)
        s.setDecFlux(dec)
        s.setRaFluxErr(raErr)
        s.setDecFluxErr(decErr)

        (ra, dec, raErr, decErr) = xyToRaDec(
            s.getXAstrom(), 
            s.getYAstrom(),
            s.getXAstromErr(), 
            s.getYAstromErr(), 
            wcs)
        s.setRaAstrom(ra);
        s.setDecAstrom(dec)
        s.setRaAstromErr(raErr);
        s.setDecAstromErr(decErr)

        # No errors for XPeak, YPeak
        coords = wcs.pixelToSky(s.getXPeak(), s.getYPeak())
        s.setRaPeak(coords.getLongitude(afwCoord.RADIANS))
        s.setDecPeak(coords.getLatitude(afwCoord.RADIANS))

        # Simple RA/decl == Astrom versions
        s.setRa(s.getRaAstrom())
        s.setRaErrForDetection(s.getRaAstromErr())
        s.setDec(s.getDecAstrom())
        s.setDecErrForDetection(s.getDecAstromErr())

def xyToRaDec(x,y, xErr, yErr, wcs, pixToSkyAffineTransform=None):
        """Use wcs to transform pixel coordinates x, y and their errors to 
        sky coordinates ra, dec with errors. If the caller does not provide an
        affine approximation to the pixel->sky WCS transform, an approximation
        is automatically computed (and used to propagate errors). For sources
        from exposures far from the poles, a single approximation can be reused
        without introducing much error.

        Note that the affine transform is expected to take inputs in units of
        pixels to outputs in units of degrees. This is an artifact of WCSLIB
        using degrees as its internal angular unit.

        Sky coordinates and their errors are returned in units of radians.

        """
        sky = wcs.pixelToSky(x, y)
        if pixToSkyAffineTransform is None:
            pixToSkyAffineTransform = wcs.linearizePixelToSky(sky)
        
        t = pixToSkyAffineTransform
        varRa  = t[0]**2 * xErr**2 + t[2]**2 * yErr**2
        varDec = t[1]**2 * xErr**2 + t[3]**2 * yErr**2
        raErr = math.radians(math.sqrt(varRa))
        decErr = math.radians(math.sqrt(varDec))

        return (sky.getLongitude(afwCoord.RADIANS),
                sky.getLatitude(afwCoord.RADIANS),
                raErr,
                decErr)
