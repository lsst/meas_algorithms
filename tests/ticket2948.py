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

class Ticket2948Test(unittest.TestCase):
    def test(self):
        visitSchema = afwTable.ExposureTable.makeMinimalSchema()
        visitGoodPixKey = visitSchema.addField("goodpix", type=int, doc="Number of good pixels in visit")
        visitWeightKey = visitSchema.addField("weight", type=float, doc="Weight for this visit")

        ccdSchema = afwTable.ExposureTable.makeMinimalSchema()
        ccdGoodPixKey = ccdSchema.addField("goodpix", type=int, doc="Number of good pixels in CCD")
        ccdWeightKey = ccdSchema.addField("weight", type=float, doc="Weight for this CCD")

        inputs = afwImage.CoaddInputs(visitSchema, ccdSchema)
        weight = numpy.nan # This is key to the bug
        numGoodPix = 123
        size = 1
        filename = "tests/ticket2948.fits"
        wcs = afwImage.makeWcs(afwCoord.Coord(0.0*afwGeom.degrees, 0.0*afwGeom.degrees),
                               afwGeom.Point2D(0.0, 0.0), 1.0e-4, 0.0, 0.0, 1.0e-4)

        # Add a CCD
        ccdRecord = inputs.ccds.addNew()
        ccdRecord.setId(1)
        ccdRecord.setI(ccdGoodPixKey, numGoodPix)
        ccdRecord.setD(ccdWeightKey, weight)

        ccdRecord.setPsf(measAlg.SingleGaussianPsf(23, 23, 2.345))
        ccdRecord.setWcs(wcs)
        ccdRecord.setBBox(afwGeom.Box2I(afwGeom.Point2I(0,0), afwGeom.Extent2I(size, size)))

        # Add a visit comprising that single CCD
        visitRecord = inputs.visits.addNew()
        visitRecord.setId(2)
        visitRecord.setPsf(ccdRecord.getPsf())
        visitRecord.setWcs(ccdRecord.getWcs())
        visitRecord.setBBox(ccdRecord.getBBox())
        visitRecord.setD(visitWeightKey, weight)

        # Attach to an Exposure
        exp = afwImage.ExposureF(size, size)
        exp.getInfo().setCoaddInputs(inputs)
        exp.setWcs(wcs)
        exp.setPsf(measAlg.CoaddPsf(inputs.ccds, exp.getWcs()))
        exp.writeFits(filename)

        # Bug present if this fails like:
        # LsstCppException: 0: lsst::pex::exceptions::RuntimeErrorException thrown at src/CoaddPsf.cc:133 in lsst::afw::geom::Point2D lsst::meas::algorithms::{anonymous}::computeAveragePosition(const ExposureCatalog&, const lsst::afw::image::Wcs&, lsst::afw::table::Key<double>)
        # 0: Message: Could not find a valid average position for CoaddPsf
        # 1: Rethrown at src/table/io/InputArchive.cc:105 in boost::shared_ptr<lsst::afw::table::io::Persistable> lsst::afw::table::io::InputArchive::Impl::get(int, const lsst::afw::table::io::InputArchive&)
        # 1: Message: loading object with id=4, name='CoaddPsf'
        fromDisk = afwImage.ExposureF(filename)
        fromDisk.getPsf().computeKernelImage()


def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(Ticket2948Test)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
