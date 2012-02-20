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
import lsst.afw.table as afwTable
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
        self.expGaussPsf = plantSources(self.nx, self.ny, self.kwid, self.sky, coordList,
                                        addPoissonNoise=False)

        # just plain sky (ie. a constant)
        self.mimg = afwImage.MaskedImageF(afwGeom.ExtentI(self.nx, self.ny))
        self.mimg.set(self.sky, 0x0, self.sky)
        self.expSky = afwImage.makeExposure(self.mimg)

        if display > 1:
            ds9.mtv(self.expGaussPsf)

    def tearDown(self):
        del self.mimg
        del self.expGaussPsf
        del self.expSky

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
        sincConfig = measAlgorithms.SincFluxConfig(radius1=0.0, radius2=0.0, angle=math.radians(theta),
                                                   ellipticity=(1-b/a))
        for r1, r2 in [(0,      0.45*a),
                       (0.45*a, 1.0*a),
                       ( 1.0*a, 2.0*a),
                       ( 2.0*a, 3.0*a),
                       ( 3.0*a, 5.0*a),
                       ( 3.0*a, 10.0*a),
                       ]:
            sincConfig.radius1 = r1
            sincConfig.radius2 = r2
            schema = afwTable.SourceTable.makeMinimalSchema()
            mp = measAlgorithms.MeasureSourcesBuilder().addAlgorithm(sincConfig.makeControl()).build(schema)

            if display:                 # draw the inner and outer boundaries of the aperture
                Mxx = 1
                Myy = (b/a)**2

                mxx, mxy, myy = c**2*Mxx + s**2*Myy, c*s*(Mxx - Myy), s**2*Mxx + c**2*Myy
                for r in (r1, r2):
                    ds9.dot("@:%g,%g,%g" % (r**2*mxx, r**2*mxy, r**2*myy), xcen, ycen, frame=frame)

            table = afwTable.SourceTable.make(schema)
            source = table.makeRecord()
            center = afwGeom.Point2D(xcen, ycen)

            mp.apply(source, objImg, center)

            self.assertAlmostEqual(math.exp(-0.5*(r1/a)**2) - math.exp(-0.5*(r2/a)**2),
                                   source["flux.sinc"]/flux, 5)

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
 
