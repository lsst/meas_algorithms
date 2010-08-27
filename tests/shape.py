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

import re, os, sys
import glob
import math
import unittest

import lsst.pex.policy as pexPolicy
import lsst.pex.exceptions as pexExceptions
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.meas.algorithms as algorithms
import lsst.utils.tests as utilsTests
import lsst.afw.detection as afwDetection

import lsst.afw.display.ds9 as ds9

try:
    type(verbose)
except NameError:
    verbose = 0
    display = False

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class ShapeTestCase(unittest.TestCase):
    """A test case for centroiding"""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testCleanup(self):
        """Test that tearDown does"""
        pass

    def testInvalidmeasureShape(self):
        """Test that we cannot instantiate an unknown measureShape"""

        def getInvalid():
            shapeFinder = algorithms.makeMeasureShape(None)
            shapeFinder.addAlgorithm("XXX")

        utilsTests.assertRaisesLsstCpp(self, pexExceptions.NotFoundException, getInvalid)

    def do_testmeasureShape(self, algorithmName):
        """Test that we can instantiate and play with a measureShape"""

        im = afwImage.ImageF(100, 100)

        #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        bkgd = 100; im.set(bkgd)
        x, y = 30, 40
        im.set(x, y, 1000 + bkgd)

        #
        # Add a Gaussian to the image
        #
        sigma_xx, sigma_xy, sigma_yy, ksize = math.pow(1.5, 2), 0, math.pow(2.5, 2), 15

        if False:
            k = afwMath.AnalyticKernel(ksize, ksize,
                                       afwMath.GaussianFunction2D(math.sqrt(sigma_xx), math.sqrt(sigma_yy)))
            cim = im.Factory(im.getDimensions())
            afwMath.convolve(cim, im, k, True)
            im = cim
        else:
            phi = math.radians(0)      # rotate +ve this far
            c, s = math.cos(phi), math.sin(phi)

            sum, sumxx, sumxy, sumyy = 4*[0.0]
            for dx in range(-ksize/2, ksize/2 + 1):
                for dy in range(-ksize/2, ksize/2 + 1):
                    u, v = c*dx + s*dy,  s*dx - c*dy
                    I = 1000*math.exp(-0.5*(u*u/sigma_xx + v*v/sigma_yy))
                    im.set(x + dx, y + dy, bkgd + I)

                    sum += I
                    sumxx += I*dx*dx
                    sumxy += I*dx*dy
                    sumyy += I*dy*dy

        sumxx /= sum; sumxy /= sum; sumyy /= sum

        if False:
            print "RHL %g %g %g" % (sumxx, sumxy, sumyy)

        msk = afwImage.MaskU(im.getDimensions()); msk.set(0)
        var = afwImage.ImageF(im.getDimensions()); var.set(10)
        im = afwImage.MaskedImageF(im, msk, var)
        del msk; del var

        shapeFinder = algorithms.makeMeasureShape(None)
        shapeFinder.addAlgorithm(algorithmName)
        shapeFinder.configure(pexPolicy.Policy(pexPolicy.PolicyString("SDSS.background: %f" % bkgd)))
            
        if display:
            ds9.mtv(im)

        shapeFinder.setImage(afwImage.makeExposure(im))

        s = shapeFinder.measure(afwDetection.Peak(x, y)).find(algorithmName)

        if False:
            Ixx, Iyy, Ixy = s.getIxx(), s.getIyy(), s.getIxy()
            A2 = 0.5*(Ixx + Iyy) + math.sqrt( (0.5*(Ixx - Iyy))**2 + Ixy**2 )
            B2 = 0.5*(Ixx + Iyy) - math.sqrt( (0.5*(Ixx - Iyy))**2 + Ixy**2 )

            print "I_xx:  %.5f %.5f" % (Ixx, sigma_xx)
            print "I_xy:  %.5f %.5f" % (Ixy, sigma_xy)
            print "I_yy:  %.5f %.5f" % (Iyy, sigma_yy)
            print "A2, B2 = %.5f, %.5f" % (A2, B2)            

        self.assertTrue(abs(x - s.getX()) < 1e-4, "%g v. %g" % (x, s.getX()))
        self.assertTrue(abs(y - s.getY()) < 1e-4, "%g v. %g" % (y, s.getY()))
        self.assertTrue(abs(s.getIxx() - sigma_xx) < 1e-3*(1 + sigma_xx))
        self.assertTrue(abs(s.getIxy() - sigma_xy) < 1e-3*(1 + sigma_xy))
        self.assertTrue(abs(s.getIyy() - sigma_yy) < 1e-3*(1 + sigma_yy))

    def testSDSSmeasureShape(self):
        """Test that we can instantiate and play with SDSSmeasureShape"""

        self.do_testmeasureShape("SDSS")

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(ShapeTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the tests"""
    utilsTests.run(suite(), exit)
 
if __name__ == "__main__":
    run(True)
