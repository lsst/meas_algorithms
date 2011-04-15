#!/usr/bin/env python

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

# -*- lsst-python -*-

import lsst.pex.policy as pexPolicy
import lsst.meas.algorithms as measAlgorithms
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetection
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.afw.display.ds9 as ds9
import math
import unittest
import lsst.utils.tests as utilsTests
import numpy

from apCorrTest import plantSources

try:
    display
except NameError:
    display = False
    displayCoeffs = False

class sincPhotSums(unittest.TestCase):

    def setUp(self):
        self.nx = 64
        self.ny = 64
        self.kwid = 15
        self.sky = 100.0
        self.val  = 10000.0
        self.sigma = 4.0
        coordList = [[self.nx/2, self.ny/2, self.val, self.sigma]]

        # exposure with gaussian
        self.expGaussPsf = plantSources(self.nx, self.ny, self.kwid, self.sky, coordList, addPoissonNoise=False)
        self.mpGaussPsf  = measAlgorithms.makeMeasurePhotometry(self.expGaussPsf)

        # just plain sky (ie. a constant)
        self.mimg = afwImage.MaskedImageF(afwGeom.ExtentI(self.nx, self.ny))
        self.mimg.set(self.sky, 0x0, self.sky)
        self.expSky = afwImage.makeExposure(self.mimg)
        self.mpSky = measAlgorithms.makeMeasurePhotometry(self.expSky)

        if display > 1:
            ds9.mtv(self.expGaussPsf)
        
        for alg in ("NAIVE", "PSF", "SINC",):
            self.mpGaussPsf.addAlgorithm(alg)
            self.mpSky.addAlgorithm(alg)

        self.pol = pexPolicy.Policy(pexPolicy.PolicyString(
            """#<?cfg paf policy?>
            NAIVE.radius: 10.0
            SINC.radius2: 6.0
            SINC.radius1: 0.0
            """
            ))

        
    def tearDown(self):
        del self.mimg
        del self.expGaussPsf
        del self.mpGaussPsf
        del self.expSky
        del self.mpSky
        del self.pol

        
        
    def XXXtestSincPhotSums(self):
        """Verify annular fluxes sum to total aperture flux."""

        # call the Photometry in this measurePhotometry object,
        #  measure the aperture flux
        def measure(mp, r1, r2, posAng, ellipticity):
            self.pol.set("SINC.radius1", r1)
            self.pol.set("SINC.radius2", r2)
            self.pol.set("SINC.angle", (numpy.pi/180.0)*posAng)
            self.pol.set("SINC.ellipticity", ellipticity)
            mp.configure(self.pol)
            peak = afwDetection.Peak(self.nx/2, self.ny/2)
            photom = mp.measure(peak)
            return photom.find("SINC").getFlux()

        # take a list of radii (sorted)
        # - for all possible annuli (ie. rad2 > rad1), measure the flux
        #   by calling the above local function measure()
        def compute(mp, rads, posAng, ellipticity, writeFits=False):
            f = {}
            for rad1 in rads:
                f[rad1] = {}
                for rad2 in rads:
                    if rad2 > rad1:
                        print "running: r1=%.1f r2=%.1f  posAng=%.1f e=%.1f" % \
                              (rad1, rad2, posAng, ellipticity)
                        f[rad1][rad2] = measure(mp, rad1, rad2, posAng, ellipticity)
                        img = measAlgorithms.getCoeffImage(rad1, rad2, posAng, ellipticity);

                        if displayCoeffs:
                            ds9.mtv(img, title="%g %g %g %g.fits" % (rad1, rad2, posAng, ellipticity))

                        if writeFits:
                            img.writeFits("cimg-%.1f-%.1f-%.1f-%.1f.fits" % (rad1, rad2, posAng, ellipticity))
            return f

        # for all the annuli we just obtained fluxes for
        #   print them as a matrix: r1 (rows), and r2 (columns)
        # The trace of the matrix corresponds to the measurements
        #   of consecutive annuli, eg 0-1, 1-2, 2-3, etc
        # So, if the method works, the aperture flux for eg. r=3 should
        #     equal the sum of annular fluxes for r = 0-1, 1-2, 2-3
        def printAndTest(f, rads, reqTolerance):
            
            fmt = "%13s "*(len(rads)+1)
            tup = ("",) + tuple(rads)
            print fmt % tup
            sumtrace = 0.0
            for i1 in range(len(rads)):
                rad1 = rads[i1]
                print "%10.8f  " % (rad1),
                for i2 in range(len(rads)):
                    rad2 = rads[i2]
                    if f[rad1].has_key(rad2):
                        print "%13.8f " % (f[rad1][rad2]),
                    else:
                        print "%13s " % (""),
                    
                    if i2 == i1 + 1:
                        sumtrace += f[rad1][rad2]
                print ""
                        
            n = len(rads)
            f0 = f[rads[0]][rads[n-1]]
            print "%.8f %.8f  %.8f" % (f0, sumtrace, (sumtrace - f0)/f0)

            self.assertTrue( abs(f0-sumtrace)/f0 < reqTolerance)
                    
                

        ######################
        # run the tests
        rads = [0.0, 2.0, 4.0, 6.0]
        posAngs = [0.0, 30.0]
        ellipticities = [0.0, 0.7]

        for i in range(len(ellipticities)):

            # sky
            # - this should be a totally bandlimited 'psf' as it's a constant
            # - there should be no flux error due to power above the nyquist
            #   and we should achieve machine precision when we compare
            #   the aperture flux to the sum of the annuli of which it's composed.
            f = compute(self.mpSky, rads, posAngs[i], ellipticities[i])
            reqTolerance = 1.0e-7 # machine precision (constant is band-limited)
            printAndTest(f, rads, reqTolerance)

            # gaussian
            # - this isn't bandlimited, and we expect to lose a bit of flux
            #   beyond the nyquist.  I selected a fairly broad gaussian (sigma=4 pixels)
            #   and it should be fairly tight in k-space, nonetheless, expect lost flux
            #   and set the tolerance lower.
            f = compute(self.mpGaussPsf, rads, posAngs[i], ellipticities[i])
            reqTolerance = 1.0e-3 # leakage due to truncation at bandlimit
            printAndTest(f, rads, reqTolerance)
        
    def testEllipticalGaussian(self):
        """Test measuring elliptical aperture mags for an elliptical Gaussian"""

        width, height = 200, 200
        xcen, ycen = 0.5*width, 0.5*height
        #
        # Make the object
        #
        gal = afwImage.ImageF(afwGeom.ExtentI(width, height))
        a, b, theta = float(10), float(5), 20
        flux = 1e4
        I0 = flux/(2*math.pi*a*b)

        c, s = math.cos(math.radians(theta)), math.sin(math.radians(theta))
        for y in range(height):
            for x in range(width):
                dx, dy = x - xcen, y - ycen
                u =  c*dx + s*dy
                v = -s*dx + c*dy
                val = I0*math.exp(-0.5*((u/a)**2 + (v/b)**2))
                if val < 0:
                    val = 0
                gal.set(x, y, val)

        objImg = afwImage.makeExposure(afwImage.makeMaskedImage(gal))
        del gal

        if display:
            frame = 0
            ds9.mtv(objImg, frame=frame, title="Elliptical")

        self.assertAlmostEqual(1.0, afwMath.makeStatistics(objImg.getMaskedImage().getImage(),
                                                           afwMath.SUM).getValue()/flux)
        #
        # Now measure some annuli
        #
        mp = measAlgorithms.makeMeasurePhotometry(objImg)
        mp.addAlgorithm("SINC")
    
        policy = pexPolicy.Policy(pexPolicy.PolicyString(
            """#<?cfg paf policy?>
            SINC.radius1: 0.0
            SINC.radius2: 0.0
            SINC.angle: %g
            SINC.ellipticity: %g
            """ % (math.radians(theta), (1 - b/a))
            ))

        peak = afwDetection.Peak(xcen, ycen)
        for r1, r2 in [(0,      0.45*a),
                       (0.45*a, 1.0*a),
                       ( 1.0*a, 2.0*a),
                       ( 2.0*a, 3.0*a),
                       ( 3.0*a, 5.0*a),
                       ( 3.0*a, 10.0*a),
                       ]:
            policy.set("SINC.radius1", r1)
            policy.set("SINC.radius2", r2)

            if display:                 # draw the inner and outer boundaries of the aperture
                Mxx = 1
                Myy = (b/a)**2

                mxx, mxy, myy = c**2*Mxx + s**2*Myy, c*s*(Mxx - Myy), s**2*Mxx + c**2*Myy
                for r in (r1, r2):
                    ds9.dot("@:%g,%g,%g" % (r**2*mxx, r**2*mxy, r**2*myy), xcen, ycen, frame=frame)

            mp.configure(policy)
            photom = mp.measure(peak)

            self.assertAlmostEqual(math.exp(-0.5*(r1/a)**2) - math.exp(-0.5*(r2/a)**2),
                                   photom.find("SINC").getFlux()/flux, 5)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(sincPhotSums)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)


def run(exit = False):
    """Run the tests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
 
