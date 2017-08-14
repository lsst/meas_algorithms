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


import lsst.afw.detection as afwDet
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
import lsst.meas.algorithms as measAlg
import lsst.utils.tests

try:
    type(display)
    import lsst.afw.display.ds9 as ds9
except NameError:
    display = False

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


class CandidateMaskingTestCase(lsst.utils.tests.TestCase):
    """Testing masking around PSF candidates."""

    def setUp(self):
        self.x, self.y = 123, 45
        self.exp = afwImage.ExposureF(256, 256)
        self.exp.getMaskedImage().getImage()[self.x, self.y] = 1.0
        self.exp.getMaskedImage().getVariance().set(0.01)

        schema = afwTable.SourceTable.makeMinimalSchema()
        self.catalog = afwTable.SourceCatalog(schema)

    def tearDown(self):
        del self.exp
        del self.catalog

    def createCandidate(self, threshold=0.1):
        """Create a PSF candidate from self.exp

        @param threshold: Threshold for creating footprints on image
        """

        source = self.catalog.addNew()
        fpSet = afwDet.FootprintSet(self.exp.getMaskedImage(), afwDet.Threshold(threshold), "DETECTED")
        if display:
            ds9.mtv(self.exp, frame=1)
            for fp in fpSet.getFootprints():
                for peak in fp.getPeaks():
                    ds9.dot("x", peak.getIx(), peak.getIy(), frame=1)

        # There might be multiple footprints; only the one around self.x,self.y should go in the source
        found = False
        for fp in fpSet.getFootprints():
            if fp.contains(afwGeom.Point2I(self.x, self.y)):
                found = True
                break
        self.assertTrue(found, "Unable to find central peak in footprint: faulty test")

        source.setFootprint(fp)
        return measAlg.PsfCandidateF(source, self.exp, self.x, self.y)

    def checkCandidateMasking(self, badPixels, extraPixels=[], size=25, threshold=0.1, pixelThreshold=0.0):
        """Check that candidates are masked properly

        We add various pixels to the image and investigate the masking.

        @param badPixels: (x,y,flux) triplet of pixels that should be masked
        @param extraPixels: (x,y,flux) triplet of additional pixels to add to image
        @param size: Size of candidate
        @param threshold: Threshold for creating footprints on image
        @param pixelThreshold: Threshold for masking pixels on candidate
        """
        image = self.exp.getMaskedImage().getImage()
        for x, y, f in badPixels + extraPixels:
            image[x, y] = f
        cand = self.createCandidate(threshold=threshold)
        oldPixelThreshold = cand.getPixelThreshold()
        try:
            cand.setPixelThreshold(pixelThreshold)
            candImage = cand.getMaskedImage(size, size)
            mask = candImage.getMask()
            if display:
                ds9.mtv(candImage, frame=2)
                ds9.mtv(candImage.getMask().convertU(), frame=3)

            detected = mask.getMaskPlane("DETECTED")
            intrp = mask.getMaskPlane("INTRP")
            for x, y, f in badPixels:
                x -= self.x - size//2
                y -= self.y - size//2
                self.assertTrue(mask.get(x, y, intrp))
                self.assertFalse(mask.get(x, y, detected))
        finally:
            # Ensure this static variable is reset
            cand.setPixelThreshold(oldPixelThreshold)

    def testBlends(self):
        """Test that blended objects are masked."""
        """
        We create another object next to the one of interest,
        joined by a bridge so that they're part of the same
        footprint.  The extra object should be masked.
        """
        self.checkCandidateMasking([(self.x+2, self.y, 1.0)], [(self.x+1, self.y, 0.5)])

    def testNeighborMasking(self):
        """Test that neighbours are masked."""
        """
        We create another object separated from the one of
        interest, which should be masked.
        """
        self.checkCandidateMasking([(self.x+5, self.y, 1.0)])

    def testFaintNeighborMasking(self):
        """Test that faint neighbours are masked."""
        """
        We create another faint (i.e., undetected) object separated
        from the one of interest, which should be masked.
        """
        self.checkCandidateMasking([(self.x+5, self.y, 0.5)], threshold=0.9, pixelThreshold=1.0)


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
