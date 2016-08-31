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

"""
Tests for PSF code

Run with:
   python testPsfAttributes.py
or
   python
   >>> import testPsfAttributes; psf.run()
"""
from __future__ import absolute_import, division, print_function
import math
import unittest

import lsst.afw.geom as afwGeom
import lsst.meas.algorithms as measAlg
import lsst.utils.tests

try:
    type(verbose)
    import lsst.afw.display.ds9 as ds9
except NameError:
    verbose = 0
    display = False

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


class PsfAttributesTestCase(lsst.utils.tests.TestCase):

    def testGaussian(self):
        """Check that we can measure a single Gaussian's attributes."""
        sigma0 = 5.0
        aEff0 = 4.0*math.pi*sigma0**2

        xwid = int(12*sigma0)
        ywid = xwid

        # set the peak of the outer guassian to 0 so this is really a single gaussian.
        psf = measAlg.SingleGaussianPsf(xwid, ywid, sigma0)

        if False and display:
            im = psf.computeImage(afwGeom.PointD(xwid//2, ywid//2))
            ds9.mtv(im, title="N(%g) psf" % sigma0, frame=0)

        psfAttrib = measAlg.PsfAttributes(psf, xwid//2, ywid//2)
        sigma = psfAttrib.computeGaussianWidth(psfAttrib.ADAPTIVE_MOMENT)
        m1 = psfAttrib.computeGaussianWidth(psfAttrib.FIRST_MOMENT)
        m2 = psfAttrib.computeGaussianWidth(psfAttrib.SECOND_MOMENT)
        noise = psfAttrib.computeGaussianWidth(psfAttrib.NOISE_EQUIVALENT)
        bick = psfAttrib.computeGaussianWidth(psfAttrib.BICKERTON)
        aEff = psfAttrib.computeEffectiveArea()

        if verbose:
            print("Adaptive            %g v %g" % (sigma0, sigma))
            print("First moment        %g v %g" % (sigma0, m1))
            print("Second moment       %g v %g" % (sigma0, m2))
            print("Noise Equivalent    %g v %g" % (sigma0, sigma))
            print("Bickerton           %g v %g" % (sigma0, bick))
            print("Effective area      %g v %f" % (aEff0, aEff))

        self.assertLessEqual(abs(sigma0 - sigma), 1.0e-2)
        self.assertLessEqual(abs(sigma0 - m1), 3.0e-2)
        self.assertLessEqual(abs(sigma0 - m2), 1.0e-2)
        self.assertLessEqual(abs(sigma0 - noise), 1.0e-2)
        self.assertLessEqual(abs(sigma0 - bick), 1.0e-2)
        self.assertLessEqual(abs(aEff0 - aEff), 1.0e-2)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
