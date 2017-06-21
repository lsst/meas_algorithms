#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2017 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as
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
from __future__ import absolute_import, division, print_function
from builtins import zip
from builtins import range
import os

import numpy as np
import unittest

import lsst.afw.math
import lsst.afw.geom
import lsst.afw.image
import lsst.afw.coord
from lsst.afw.geom.polygon import Polygon
import lsst.meas.algorithms
import lsst.pex.exceptions
import lsst.utils.tests


class CoaddBoundedFieldTestCase(lsst.utils.tests.TestCase):

    def makeRandomWcs(self, crval, maxOffset=10.0):
        """Make a random TAN Wcs that's complex enough to create an interesting test, by combining
        random offsets and rotations.  We don't include any random scaling, because we don't want
        the Jacobian to be significantly different from one in these tests, as we want to compare
        warped images (which implicitly include the Jacobian) against the behavior of CoaddBoundedField
        (which intentionally does not).
        """
        crpix = lsst.afw.geom.Point2D(*np.random.uniform(low=-maxOffset, high=maxOffset, size=2))
        rotate = lsst.afw.geom.LinearTransform.makeRotation(
            np.pi*np.random.rand()*lsst.afw.geom.radians
        )
        scale = lsst.afw.geom.LinearTransform.makeScaling(
            (0.01*lsst.afw.geom.arcseconds).asDegrees()
        )
        cd = rotate * scale
        return lsst.afw.image.makeWcs(crval, crpix, cd[cd.XX], cd[cd.XY], cd[cd.YX], cd[cd.YY])

    def makeRandomField(self, bbox):
        """Create a random ChebyshevBoundedField"""
        coefficients = np.random.randn(3, 3)
        # We add a positive DC offset to make sure our fields more likely to combine constructively;
        # this makes the test more stringent, because we get less accidental small numbers.
        coefficients[0, 0] = np.random.uniform(low=4, high=6)
        return lsst.afw.math.ChebyshevBoundedField(bbox, coefficients)

    def constructElements(self, validBox):
        """Construct the elements of a CoaddBoundedField."""
        np.random.seed(50)
        validPolygon = Polygon(lsst.afw.geom.Box2D(validBox)) if validBox else None
        elements = []
        validBoxes = []

        for i in range(10):
            elements.append(
                lsst.meas.algorithms.CoaddBoundedField.Element(
                    self.makeRandomField(self.elementBBox),
                    self.makeRandomWcs(self.crval),
                    validPolygon,
                    np.random.uniform(low=0.5, high=1.5)
                )
            )
            validBoxes.append(validBox)
        return elements, validBoxes

    def setUp(self):
        self.crval = lsst.afw.coord.IcrsCoord(45.0*lsst.afw.geom.degrees, 45.0*lsst.afw.geom.degrees)
        self.elementBBox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-50, -50), lsst.afw.geom.Point2I(50, 50))
        self.coaddWcs = self.makeRandomWcs(self.crval, maxOffset=0.0)
        self.bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-75, -75), lsst.afw.geom.Point2I(75, 75))
        self.possibleValidBoxes = (None, lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-25, -25),
                                         lsst.afw.geom.Point2I(25, 25)))

    def testEvaluate(self):
        """Test the main implementation of CoaddBoundedField::evaluate() by creating images of
        each constituent field, warping them, and coadding them to form a check image.

        We run this test twice: once with a regular ValidPolygon, and once
        with the polygon undefined (ie, set to None); see DM-11008.
        """
        def runTest(elements, validBoxes):
            field = lsst.meas.algorithms.CoaddBoundedField(self.bbox, self.coaddWcs, elements, 0.0)
            coaddImage = lsst.afw.image.ImageF(self.bbox)
            warpCtrl = lsst.afw.math.WarpingControl("bilinear")
            weightMap = lsst.afw.image.ImageF(self.bbox)
            for element, validBox in zip(elements, validBoxes):
                elementImage = (lsst.afw.image.ImageF(validBox) if validBox
                                else lsst.afw.image.ImageF(element.field.getBBox()))
                element.field.fillImage(elementImage, True)
                warp = lsst.afw.image.ImageF(self.bbox)
                lsst.afw.math.warpImage(warp, self.coaddWcs, elementImage, element.wcs, warpCtrl, 0.0)
                coaddImage.scaledPlus(element.weight, warp)
                warp.getArray()[warp.getArray() != 0.0] = element.weight
                weightMap += warp
            coaddImage /= weightMap
            coaddImage.getArray()[np.isnan(coaddImage.getArray())] = 0.0
            fieldImage = lsst.afw.image.ImageF(self.bbox)
            field.fillImage(fieldImage)

            # The coaddImage we're comparing to has artifacts at the edges of all the input exposures,
            # due to the interpolation, so we just check that the number of unequal pixels is small (<10%)
            # This can depend on the random number seed, so if a failure is seen, uncomment the line below
            # to enable a plot of the differences, and verify by eye that they look like edge artifacts.  If
            # they do, just modify the seed (at the top of this file) or change number-of-pixels threshold
            # until the test passes.

            diff = np.abs(fieldImage.getArray() - coaddImage.getArray())
            relTo = np.abs(fieldImage.getArray())
            rtol = 1E-2
            atol = 1E-7
            bad = np.logical_and(diff > rtol*relTo, diff > atol)

            if False:   # enable this to see a plot of the comparison (but it will always fail, since
                        # it doesn't take into account the artifacts in coaddImage)
                self.assertFloatsAlmostEqual(fieldImage.getArray(), coaddImage.getArray(), plotOnFailure=True,
                                             rtol=rtol, atol=atol, relTo=relTo)

            self.assertLess(bad.sum(), 0.10*self.bbox.getArea())

        for validBox in self.possibleValidBoxes:
            runTest(*self.constructElements(validBox))

    def testEquality(self):
        def runTest(elements, validBoxes):
            field1 = lsst.meas.algorithms.CoaddBoundedField(self.bbox, self.coaddWcs, elements, 0.0)
            field2 = lsst.meas.algorithms.CoaddBoundedField(self.bbox, self.coaddWcs, elements, 0.0)
            self.assertEqual(field1, field2)

            bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-75, -75), lsst.afw.geom.Point2I(75, 75))
            field3 = lsst.meas.algorithms.CoaddBoundedField(bbox, self.coaddWcs, elements, 0.0)
            self.assertEqual(field1, field3)

            coaddWcs = self.coaddWcs.clone()
            field4 = lsst.meas.algorithms.CoaddBoundedField(self.bbox, coaddWcs, elements, 0.0)
            self.assertEqual(field1, field4)

            # NOTE: make a copy of the list; [:] to copy segfaults,
            # copy.copy() doesn't behave nicely on py2 w/our pybind11 objects,
            # and .elements.copy() doesn't exist on py2.
            elementsCopy = list(elements)
            field5 = lsst.meas.algorithms.CoaddBoundedField(self.bbox, self.coaddWcs, elementsCopy, 0.0)
            self.assertEqual(field1, field5)

            # inequality tests below here
            field6 = lsst.meas.algorithms.CoaddBoundedField(self.bbox, self.coaddWcs, elements, 3.0)
            self.assertNotEqual(field1, field6)
            field7 = lsst.meas.algorithms.CoaddBoundedField(self.bbox, self.coaddWcs, [], 0.0)
            self.assertNotEqual(field1, field7)

            bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-74, -75), lsst.afw.geom.Point2I(75, 75))
            field8 = lsst.meas.algorithms.CoaddBoundedField(bbox, self.coaddWcs, elements, 0.0)
            self.assertNotEqual(field1, field8)

            crval = lsst.afw.coord.IcrsCoord(45.0*lsst.afw.geom.degrees, 45.0*lsst.afw.geom.degrees)
            coaddWcs = self.makeRandomWcs(crval, maxOffset=2.0)
            field9 = lsst.meas.algorithms.CoaddBoundedField(self.bbox, coaddWcs, elements, 0.0)
            self.assertNotEqual(field1, field9)

            elementsCopy = list(elements)
            elementsCopy[2].weight = 1000000
            field10 = lsst.meas.algorithms.CoaddBoundedField(self.bbox, self.coaddWcs, elementsCopy, 0.0)
            self.assertNotEqual(field1, field10)
            elementsCopy = list(elements)
            elementsCopy.pop(0)
            field11 = lsst.meas.algorithms.CoaddBoundedField(self.bbox, self.coaddWcs, elementsCopy, 0.0)
            self.assertNotEqual(field1, field11)

        for validBox in self.possibleValidBoxes:
            runTest(*self.constructElements(validBox))

    def testPersistence(self):
        """Test that we can round-trip a CoaddBoundedField through FITS."""
        def runTest(elements, validBoxes):
            filename = "testCoaddBoundedField.fits"
            field1 = lsst.meas.algorithms.CoaddBoundedField(self.bbox, self.coaddWcs, elements, 0.0)
            field1.writeFits(filename)
            field2 = lsst.meas.algorithms.CoaddBoundedField.readFits(filename)
            image1 = lsst.afw.image.ImageD(self.bbox)
            image2 = lsst.afw.image.ImageD(self.bbox)
            field1.fillImage(image1)
            field2.fillImage(image2)
            # use assertFloatsAlmostEqual for array support, not fuzziness; this test should be exact
            self.assertFloatsAlmostEqual(image1.getArray(), image2.getArray(), rtol=0.0, atol=0.0,
                                         plotOnFailure=True)
            os.remove(filename)

        for validBox in self.possibleValidBoxes:
            runTest(*self.constructElements(validBox))

    def tearDown(self):
        del self.coaddWcs
        del self.bbox


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
