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

import re
import os
import glob
import math
import pdb                          # we may want to say pdb.set_trace()
import unittest

import eups
import lsst.pex.exceptions as pexExceptions
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.detection as afwDetection
import lsst.meas.algorithms as algorithms
import lsst.utils.tests as utilsTests

import testLib

try:
    type(verbose)
except NameError:
    verbose = 0

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class CentroidTestCase(unittest.TestCase):
    """A test case for centroiding"""

    def setUp(self):
	pass

    def tearDown(self):
        pass

    def testMeasureCentroid(self):
        """Test that we can instantiate and play with SillyMeasureCentroid"""

        for imageFactory in (afwImage.MaskedImageF,
                             afwImage.MaskedImageD,
                             ):
            im = imageFactory(afwGeom.ExtentI(100, 100))
            exp = afwImage.makeExposure(im)
            centroider =  algorithms.makeMeasureAstrometry(exp)
            centroider.addAlgorithm("SILLY")
            
            x, y = 10, 20
            c = centroider.measure(afwDetection.Source(0), exp, afwGeom.Point2D(x, y)).find()
            self.assertEqual(x, c.getX() - 1)
            self.assertEqual(y, c.getY() - 1)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(CentroidTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)
 
if __name__ == "__main__":
    run(True)
