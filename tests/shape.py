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
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
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

    def do_testmeasureShape(self):
        """Test that we can instantiate and play with a measureShape"""

        algorithmName = "shape.sdss"
        algorithmConfig = algorithms.SdssShapeConfig()

        im = afwImage.ImageF(afwGeom.ExtentI(100))

        #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        bkgd = 100.0; im.set(bkgd)
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
        exp = afwImage.makeExposure(im)

        algorithmConfig.background = bkgd
        schema = afwTable.SourceTable.makeMinimalSchema()
        shapeFinder = algorithms.MeasureSourcesBuilder()\
            .addAlgorithm(algorithmConfig.makeControl())\
            .build(schema)
            
        if display:
            ds9.mtv(im)

        table = afwTable.SourceTable.make(schema)
        table.defineShape(algorithmName)
        table.defineCentroid(algorithmName + ".centroid")
        source = table.makeRecord()
        center = afwGeom.Point2D(x, y)

        shapeFinder.apply(source, exp, center)

        if False:
            Ixx, Iyy, Ixy = source.getIxx(), source.getIyy(), source.getIxy()
            A2 = 0.5*(Ixx + Iyy) + math.sqrt( (0.5*(Ixx - Iyy))**2 + Ixy**2 )
            B2 = 0.5*(Ixx + Iyy) - math.sqrt( (0.5*(Ixx - Iyy))**2 + Ixy**2 )

            print "I_xx:  %.5f %.5f" % (Ixx, sigma_xx)
            print "I_xy:  %.5f %.5f" % (Ixy, sigma_xy)
            print "I_yy:  %.5f %.5f" % (Iyy, sigma_yy)
            print "A2, B2 = %.5f, %.5f" % (A2, B2)            

        self.assertTrue(abs(x - source.getX()) < 1e-4, "%g v. %g" % (x, source.getX()))
        self.assertTrue(abs(y - source.getY()) < 1e-4, "%g v. %g" % (y, source.getY()))
        self.assertTrue(abs(source.getIxx() - sigma_xx) < 1e-3*(1 + sigma_xx),
                        "%g v. %g" % (sigma_xx, source.getIxx()))
        self.assertTrue(abs(source.getIxy() - sigma_xy) < 1e-3*(1 + sigma_xy),
                        "%g v. %g" % (sigma_xy, source.getIxy()))
        self.assertTrue(abs(source.getIyy() - sigma_yy) < 1e-3*(1 + sigma_yy),
                        "%g v. %g" % (sigma_yy, source.getIyy()))

    def testSDSSmeasureShape(self):
        """Test that we can instantiate and play with SDSSmeasureShape"""

        self.do_testmeasureShape()

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
