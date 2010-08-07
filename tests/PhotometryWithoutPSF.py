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
"""
Tests for ticket 1043 - Photometry fails when no PSF is provided
"""

import lsst.pex.policy as pexPolicy
import lsst.meas.algorithms as measAlgorithms
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetection

import math
import unittest
import lsst.utils.tests as utilsTests

# math.isnan() available in 2.6, but not 2.5.2
try:
    math.isnan(1.0)
except AttributeError:
    math.isnan = lambda x: x != x

class ticket1043TestCase(unittest.TestCase):

    def setUp(self):
        self.mi = afwImage.MaskedImageF(100, 100)
        self.mi.set(0, 0x0, 1)

        self.measurePhotom = measAlgorithms.MeasurePhotometryF(afwImage.makeExposure(self.mi))

        for alg in ("NAIVE", "PSF", "SINC",):
            self.measurePhotom.addAlgorithm(alg)

        pol = pexPolicy.Policy(pexPolicy.PolicyString(
            """#<?cfg paf policy?>
            NAIVE.radius: 10.0
            SINC.radius: 3.0
            """
            ))
            
        self.measurePhotom.configure(pol)

    def tearDown(self):
        del self.mi
        del self.measurePhotom

    def testTicket1043(self):
        self.mi.set(50, 50, (1, 0x0, 1))
        peak = afwDetection.Peak(50, 50)

        photom = self.measurePhotom.measure(peak)

        # make sure aperture photometry works
        knownSincApFlux = 1.0111901760101318
        
        self.assertEqual(photom.find("NAIVE").getFlux(), 1.0)
        self.assertEqual(photom.find("SINC").getFlux(),  knownSincApFlux)
        self.assertTrue(math.isnan(photom.find("PSF").getFlux()))

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(ticket1043TestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)


def run(exit = False):
    """Run the tests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
 
