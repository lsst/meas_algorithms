#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2014 LSST Corporation.
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

import os, sys

import numpy
import unittest

import lsst.utils.tests
import lsst.pex.exceptions
import lsst.afw.math
import lsst.afw.geom
import lsst.afw.image
import lsst.afw.coord
import lsst.meas.algorithms

numpy.random.seed(50)

class CoaddBoundedFieldTestCase(lsst.utils.tests.TestCase):

    def makeRandomWcs(self, crval, maxOffset=10.0):
        """Make a random TAN Wcs that's complex enough to create an interesting test, by combining
        random offsets and rotations.  We don't include any random scaling, because we don't want
        the Jacobian to be significantly different from one in these tests, as we want to compare
        warped images (which implicitly include the Jacobian) against the behavior of CoaddBoundedField
        (which intentionally does not).
        """
        crpix = lsst.afw.geom.Point2D(*numpy.random.uniform(low=-maxOffset, high=maxOffset, size=2))
        rotate = lsst.afw.geom.LinearTransform.makeRotation(
            numpy.pi*numpy.random.rand()*lsst.afw.geom.radians
            )
        scale = lsst.afw.geom.LinearTransform.makeScaling(
            (0.01*lsst.afw.geom.arcseconds).asDegrees()
            )
        cd = rotate * scale
        return lsst.afw.image.makeWcs(crval, crpix, cd[cd.XX], cd[cd.XY], cd[cd.YX], cd[cd.YY])

    def makeRandomField(self, bbox):
        """Create a random ChebyshevBoundedField"""
        coefficients = numpy.random.randn(3, 3)
        # We add a positive DC offset to make sure our fields more likely to combine constructively;
        # this makes the test more stringent, because we get less accidental small numbers.
        coefficients[0,0] = numpy.random.uniform(low=4, high=6)
        return lsst.afw.math.ChebyshevBoundedField(bbox, coefficients)

    def setUp(self):
        crval = lsst.afw.coord.IcrsCoord(45.0*lsst.afw.geom.degrees, 45.0*lsst.afw.geom.degrees)
        elementBBox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-50, -50), lsst.afw.geom.Point2I(50, 50))
        self.elements = lsst.meas.algorithms.CoaddBoundedField.ElementVector()
        for i in range(10):
            self.elements.append(
                lsst.meas.algorithms.CoaddBoundedField.Element(
                    self.makeRandomField(elementBBox),
                    self.makeRandomWcs(crval),
                    numpy.random.uniform(low=0.5, high=1.5)
                    )
                )
        self.coaddWcs = self.makeRandomWcs(crval, maxOffset=0.0)
        self.bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-75, -75), lsst.afw.geom.Point2I(75, 75))

    def testEvaluate(self):
        """Test the main implementation of CoaddBoundedField::evaluate() by creating images of
        each constituent field, warping them, and coadding them to form a check image.
        """
        field = lsst.meas.algorithms.CoaddBoundedField(self.bbox, self.coaddWcs, self.elements, 0.0)
        coaddImage = lsst.afw.image.ImageF(self.bbox)
        warpCtrl = lsst.afw.math.WarpingControl("bilinear")
        weightMap = lsst.afw.image.ImageF(self.bbox)
        for element in self.elements:
            elementImage = lsst.afw.image.ImageF(element.field.getBBox())
            element.field.fillImage(elementImage)
            warp = lsst.afw.image.ImageF(self.bbox)
            lsst.afw.math.warpImage(warp, self.coaddWcs, elementImage, element.wcs, warpCtrl, 0.0)
            coaddImage.scaledPlus(element.weight, warp)
            warp.getArray()[warp.getArray() != 0.0] = element.weight
            weightMap += warp
        coaddImage /= weightMap
        coaddImage.getArray()[numpy.isnan(coaddImage.getArray())] = 0.0
        fieldImage = lsst.afw.image.ImageF(self.bbox)
        field.fillImage(fieldImage)

        # The coaddImage we're comparing to has artifacts at the edges of all the input exposures,
        # due to the interpolation, so we just check that the number of unequal pixels is small (<10%)
        # This can depend on the random number seed, so if a failure is seen, uncomment the line below
        # to enable a plot of the differences, and verify by eye that they look like edge artifacts.  If
        # they do, just modify the seed (at the top of this file) or change number-of-pixels threshold
        # until the test passes.

        diff = numpy.abs(fieldImage.getArray() - coaddImage.getArray())
        relTo = numpy.abs(fieldImage.getArray())
        rtol = 1E-2
        atol = 1E-7
        bad = numpy.logical_and(diff > rtol*relTo, diff > atol)

        if False:   # enable this to see a plot of the comparison (but it will always fail, since
                    # it doesn't take into account the artifacts in coaddImage)
            self.assertClose(fieldImage.getArray(), coaddImage.getArray(), plotOnFailure=True,
                             rtol=rtol, atol=atol, relTo=relTo)

        self.assertLess(bad.sum(), 0.10*self.bbox.getArea())


    def testPersistence(self):
        """Test that we can round-trip a CoaddBoundedField through FITS."""
        filename = "testCoaddBoundedField.fits"
        field1 = lsst.meas.algorithms.CoaddBoundedField(self.bbox, self.coaddWcs, self.elements, 0.0)
        field1.writeFits(filename)
        field2 = lsst.meas.algorithms.CoaddBoundedField.readFits(filename)
        image1 = lsst.afw.image.ImageD(self.bbox)
        image2 = lsst.afw.image.ImageD(self.bbox)
        field1.fillImage(image1)
        field2.fillImage(image2)
        # use assertClose for array support, not fuzziness; this test should be exact
        self.assertClose(image1.getArray(), image2.getArray(), rtol=0.0, atol=0.0, plotOnFailure=True)
        os.remove(filename)


    def tearDown(self):
        del self.coaddWcs
        del self.bbox
        del self.elements

def suite():
    """Returns a suite containing all the test cases in this module."""
    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(CoaddBoundedFieldTestCase)
    return unittest.TestSuite(suites)

def run(exit = False):
    lsst.utils.tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
