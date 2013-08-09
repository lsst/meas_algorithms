#!/usr/bin/env python

#
# LSST Data Management System
# Copyright 2008-2013 LSST Corporation.
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

import os, os.path
import unittest
import numpy

import lsst.utils.tests as utilsTests

import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord
import lsst.afw.table as afwTable
import lsst.meas.algorithms as measAlg

class Ticket2986Test(unittest.TestCase):
    def test(self):
        schema = afwTable.ExposureTable.makeMinimalSchema()
        schema.addField("ccd", int, doc="CCD number")
        schema.addField("visit", long, doc="Visit number")
        schema.addField("goodpix", int, doc="Number of good pixels")
        schema.addField("weight", float, doc="Weighting for this CCD")
        ccds = afwTable.ExposureCatalog(schema)

        wcs = afwImage.makeWcs(afwCoord.Coord(0.0*afwGeom.degrees, 0.0*afwGeom.degrees),
                               afwGeom.Point2D(0.0, 0.0), 1.0e-4, 0.0, 0.0, 1.0e-4)

        new = ccds.addNew()
        new.set("id", 0)
        new.set("bbox.min", afwGeom.Point2I(0,0))
        new.set("bbox.max", afwGeom.Point2I(1024,1024))

        # The following lines are critical for reproducing the bug, because
        # the code is reading a double starting at the 'ccd' (offset 24), and
        # it sees a zero (from the zero in 'ccd' and the leading zeros in 'visit').
        new.set("ccd", 0)
        new.set("visit", 6789)

        new.set("goodpix", 987654321)
        new.set("weight", 1.0)
        new.setPsf(measAlg.SingleGaussianPsf(23, 23, 2.345))
        new.setWcs(wcs)

        # The following fails with:
        # LsstCppException: 0: lsst::pex::exceptions::RuntimeErrorException thrown at src/CoaddPsf.cc:134 in lsst::afw::geom::Point2D lsst::meas::algorithms::{anonymous}::computeAveragePosition(const ExposureCatalog&, const lsst::afw::image::Wcs&, lsst::afw::table::Key<double>)
        # 0: Message: Could not find a valid average position for CoaddPsf
        measAlg.CoaddPsf(ccds, wcs)

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(Ticket2986Test)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
