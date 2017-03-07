#!/usr/bin/env python
#
# LSST Data Management System
#
# Copyright 2008-2016  AURA/LSST.
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
# see <https://www.lsstcorp.org/LegalNotices/>.
#
from __future__ import absolute_import, division, print_function
import unittest


import numpy as np
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord
import lsst.afw.table as afwTable
import lsst.meas.algorithms as measAlg
import lsst.utils.tests


class Ticket2986Test(unittest.TestCase):

    def test(self):
        schema = afwTable.ExposureTable.makeMinimalSchema()
        schema.addField("ccd", np.int32, doc="CCD number")
        schema.addField("visit", np.int32, doc="Visit number")
        schema.addField("goodpix", np.int32, doc="Number of good pixels")
        schema.addField("weight", float, doc="Weighting for this CCD")
        ccds = afwTable.ExposureCatalog(schema)

        wcs = afwImage.makeWcs(afwCoord.Coord(0.0*afwGeom.degrees, 0.0*afwGeom.degrees),
                               afwGeom.Point2D(0.0, 0.0), 1.0e-4, 0.0, 0.0, 1.0e-4)

        new = ccds.addNew()
        new.set("id", 0)
        new.set("bbox_min_x", 0)
        new.set("bbox_min_y", 0)
        new.set("bbox_max_x", 1024)
        new.set("bbox_max_y", 1024)

        # The following lines are critical for reproducing the bug, because
        # the code is reading a double starting at the 'ccd' (offset 24), and
        # it sees a zero (from the zero in 'ccd' and the leading zeros in 'visit').
        new.set("ccd", 0)
        new.set("visit", 6789)

        new.set("goodpix", 987654321)
        new.set("weight", 1.0)
        new.setPsf(measAlg.SingleGaussianPsf(23, 23, 2.345))
        new.setWcs(wcs)

        # In the presence of the bug, the following fails with
        # lsst::pex::exceptions::RuntimeError thrown in src/CoaddPsf.cc
        # with message: "Could not find a valid average position for CoaddPsf"
        measAlg.CoaddPsf(ccds, wcs)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
