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
import lsst.afw.display.ds9 as ds9
import math
import unittest
import lsst.utils.tests as utilsTests
import numpy

from apCorrTest import plantSources

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
        self.mpGaussPsf  = measAlgorithms.MeasurePhotometryF(self.expGaussPsf)

        # just plain sky (ie. a constant)
        self.mimg = afwImage.MaskedImageF(self.nx, self.ny)
        self.mimg.set(self.sky, 0x0, self.sky)
        self.expSky = afwImage.makeExposure(self.mimg)
        self.mpSky = measAlgorithms.MeasurePhotometryF(self.expSky)

        if False:
            ds9.mtv(self.exposure)
        
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

        
        
    def testSincPhotSums(self):
        """Verify annular fluxes sum to total aperture flux."""

        # call the Photometry in this measurePhotometry object,
        #  measure the aperture flux
        def measure(mp, r1, r2, posAng, ecc):
            self.pol.set("SINC.radius1", r1)
            self.pol.set("SINC.radius2", r2)
            self.pol.set("SINC.positionAngle", (numpy.pi/180.0)*posAng)
            self.pol.set("SINC.ellipticity", ecc)
            mp.configure(self.pol)
            peak = afwDetection.Peak(self.nx/2, self.ny/2)
            photom = mp.measure(peak)
            return photom.find("SINC").getFlux()

        # take a list of radii (sorted)
        # - for all possible annuli (ie. rad2 > rad1), measure the flux
        #   by calling the above local function measure()
        def compute(mp, rads, posAng, ecc, writeFits=False):
            f = {}
            for rad1 in rads:
                f[rad1] = {}
                for rad2 in rads:
                    if rad2 > rad1:
                        print "running: r1=%.1f r2=%.1f  posAng=%.1f e=%.1f" % (rad1, rad2, posAng, ecc)
                        f[rad1][rad2] = measure(mp, rad1, rad2, posAng, ecc)
                        img = measAlgorithms.getCoeffImage(rad1, rad2, posAng, ecc);
                        if writeFits:
                            img.writeFits("cimg-%.1f-%.1f-%.1f-%.1f.fits" % (rad1, rad2, posAng, ecc))
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
        eccs = [0.0, 0.7]

        for i in range(len(eccs)):

            # sky
            # - this should be a totally bandlimited 'psf' as it's a constant
            # - there should be no flux error due to power above the nyquist
            #   and we should achieve machine precision when we compare
            #   the aperture flux to the sum of the annuli of which it's composed.
            f = compute(self.mpSky, rads, posAngs[i], eccs[i])
            reqTolerance = 1.0e-7 # machine precision (constant is band-limited)
            printAndTest(f, rads, reqTolerance)

            # gaussian
            # - this isn't bandlimited, and we expect to lose a bit of flux
            #   beyond the nyquist.  I selected a fairly broad gaussian (sigma=4 pixels)
            #   and it should be fairly tight in k-space, nonetheless, expect lost flux
            #   and set the tolerance lower.
            f = compute(self.mpGaussPsf, rads, posAngs[i], eccs[i])
            reqTolerance = 1.0e-3 # leakage due to truncation at bandlimit
            printAndTest(f, rads, reqTolerance)



        
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
 
