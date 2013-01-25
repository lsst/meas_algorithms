
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
Tests for CoaddPsf code

Run with:
   python CoaddPsf.py
"""

import os, sys
from math import *
import numpy
import unittest
import eups
import lsst.utils.tests as utilsTests
import lsst.pex.exceptions as pexExceptions
import lsst.pex.logging as logging

import math
import pdb
import numpy

import lsst.pex.config as pexConfig
import lsst.afw.geom as afwGeom
import lsst.afw.detection as afwDetection
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.coord as afwCoord
import lsst.pipe.base as pipeBase
import lsst.afw.cameraGeom as cameraGeom
import lsst.meas.algorithms as measAlg

try:
    type(verbose)
except NameError:
    verbose = 0
    logging.Trace.setVerbosity("meas.algorithms.Interp", verbose)
    logging.Trace.setVerbosity("afw.detection.Psf", verbose)
    display = False

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
def getPsfAttributes(psf, point):
        psfAttrib = measAlg.PsfAttributes(psf, point)
        sigma = psfAttrib.computeGaussianWidth(psfAttrib.ADAPTIVE_MOMENT)
        m1    = psfAttrib.computeGaussianWidth(psfAttrib.FIRST_MOMENT)
        m2    = psfAttrib.computeGaussianWidth(psfAttrib.SECOND_MOMENT)
        aEff  = psfAttrib.computeEffectiveArea()
        return sigma,m1, m2

def testRelDiff(A,B,delta):
    retval = abs((A-B)/(.5*(A+B)))
    if (retval > delta):
        return False
    return True

# Check that the coaddpsf passed in is really using the correct components and weighting them properly
# The components in this case are all single gaussians, and we will just add the moments

def getCoaddPsfAttributes(coaddpsf, point):
    count = coaddpsf.getComponentCount()
    coaddWcs = coaddpsf.getCoaddWcs()
    weight_sum = 0.0
    m1_sum = 0.0
    m2_sum = 0.0
    components = []
    for i in range(count):
        wcs = coaddpsf.getWcs(i)
        psf = coaddpsf.getPsf(i)
        bbox = afwGeom.Box2D(coaddpsf.getBBox(i))
        point_rel = wcs.skyToPixel(coaddWcs.pixelToSky(afwGeom.Point2D(point)))
        if bbox.contains(point_rel):
            #print i, "contains", point_rel   
            weight = coaddpsf.getWeight(i)
            attribs = measAlg.PsfAttributes(psf, afwGeom.Point2I(0,0))
            m1 = attribs.computeGaussianWidth(attribs.FIRST_MOMENT)
            m2 = attribs.computeGaussianWidth(attribs.SECOND_MOMENT)
            m1_sum += m1*weight
            m2_sum += m2*weight
            weight_sum += weight
            components.append((coaddpsf.getId(i), m1, m2, weight))
        else: 
            #print i, "doesn't contain", point_rel    
            pass
    if weight_sum == 0.0:
        return 0,0,0
    else:
        return m1_sum/weight_sum, m2_sum/weight_sum, components

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#   This is a test which checks to see that all of the ExposureCatalog rows are correctly
#   ingested by the CoaddPsf constructor, and that they can be read back in the right order
#   and with the right values
#   The weightname mechanism is also tested.  Whatever input column name is used should be
#   mapped to "weight"

class CreatePsfTest(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test(self):
        """Check that we can create a CoaddPsf with 10 elements"""
        # this is the coadd Wcs we want
        cd11 = 5.55555555e-05
        cd12 = 0.0
        cd21 = 0.0
        cd22 = 5.55555555e-05
        crval1 = 0.0
        crval2 = 0.0
        crpix = afwGeom.PointD(1000, 1000)
        crval = afwCoord.Coord(afwGeom.Point2D(crval1, crval2))
        wcsref = afwImage.makeWcs(crval,crpix,cd11,cd12,cd21,cd22)


        schema = afwTable.ExposureTable.makeMinimalSchema()
        schema.addField("customweightname", type="D", doc="Coadd weight")
        mycatalog = afwTable.ExposureCatalog(schema)

#       Imagine a ccd in each of positions +-1000 pixels from the center
        for i in range(1,10,1):
            record = mycatalog.getTable().makeRecord()
            psf = afwDetection.createPsf("DoubleGaussian", 100, 100, i, 1.00, 0.0);
            record.setPsf(psf)
            crpix = afwGeom.PointD(i*1000.0, i*1000.0)
            wcs = afwImage.makeWcs(crval,crpix,cd11,cd12,cd21,cd22)

            record.setWcs(wcs)
            record['customweightname'] = 1.0 * (i+1)
            record['id'] = i
            bbox = afwGeom.Box2I(afwGeom.Point2I(0,0), afwGeom.Extent2I(i*1000, i*1000))
            record.setBBox(bbox)
            mycatalog.append(record)

        mypsf = measAlg.CoaddPsf(mycatalog, wcsref, 'customweightname')
        mypsf.setDefaultImageSize(afwGeom.Extent2I(100, 100))

        print "OK, starting create test"
        
        self.assertTrue(mypsf.getComponentCount() == 9)
        for i in range(1,10,1):
            wcs = mypsf.getWcs(i-1)
            psf = mypsf.getPsf(i-1)
            bbox = mypsf.getBBox(i-1)
            weight = mypsf.getWeight(i-1)
            id = mypsf.getId(i-1)
            self.assertTrue(i == id)
            self.assertTrue(weight == 1.0*(i+1))
            self.assertTrue(bbox.getEndX() == 1000* i)
            self.assertTrue(wcs.getPixelOrigin().getX() == (1000.0 * i))
            sigma,m1,m1 = getPsfAttributes(psf, afwGeom.Point2I(0,0))
            self.assertTrue(testRelDiff(i,sigma,.01))

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#  This test is no longer functional because it depends on a Biaxil Gaussian Psf
#  which I need to implement in Python instead of C++
class RotatePsfTest(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test(self):
        """Check that we can create a CoaddPsf with 10 elements"""
        # this is the coadd Wcs we want
        cd11 = 5.55555555e-05
        cd12 = 0.0
        cd21 = 0.0
        cd22 = 5.55555555e-05
        crval1 = 0.0
        crval2 = 0.0
        crpix = afwGeom.PointD(1000, 1000)
        crval = afwCoord.Coord(afwGeom.Point2D(crval1, crval2))
        wcsref = afwImage.makeWcs(crval,crpix,cd11,cd12,cd21,cd22)

        cd21 = 5.55555555e-05
        cd12 = 5.55555555e-05
        cd11 = 0.0
        cd22 = 0.0

        wcs = afwImage.makeWcs(crval, crpix, cd11, cd12, cd21,cd22)
        schema = afwTable.ExposureTable.makeMinimalSchema()
        schema.addField("weight", type="D", doc="Coadd weight")
        mycatalog = afwTable.ExposureCatalog(schema)
        # make a single record with an oblong Psf
        record = mycatalog.getTable().makeRecord()
        psf = measAlg.BiaxialGaussianPsf(100,100,1.0,6.0,0.0)
        record.setPsf(psf)
        img = psf.computeImage(afwGeom.PointD(0,0))
        img.writeFits("img1.fits")
        record.setWcs(wcs)
        record['weight'] = 1.0
        record['id'] = 1
        bbox = afwGeom.Box2I(afwGeom.Point2I(0,0), afwGeom.Extent2I(2000, 2000))
        record.setBBox(bbox)
        mycatalog.append(record) 
        mypsf = measAlg.CoaddPsf(mycatalog, wcsref)
        mypsf.setDefaultImageSize(afwGeom.Extent2I(100, 100))
        img = mypsf.computeImage(afwGeom.Point2D(1000,1000), afwGeom.Extent2I(100,100), True, False)
        img.writeFits("img2.fits")
        psfAttrib = measAlg.PsfAttributes(psf, afwGeom.Point2I(1000,1001))
        m1    = psfAttrib.computeGaussianWidth(psfAttrib.FIRST_MOMENT)
        m2    = psfAttrib.computeGaussianWidth(psfAttrib.SECOND_MOMENT)
        coaddpsfAttrib = measAlg.PsfAttributes(mypsf, afwGeom.Point2I(1000,1001))
        coaddm1    = coaddpsfAttrib.computeGaussianWidth(psfAttrib.FIRST_MOMENT)
        coaddm2    = coaddpsfAttrib.computeGaussianWidth(psfAttrib.SECOND_MOMENT)
        print "m1: ",m1,coaddm1,"m2: ",m2,coaddm2
        self.assertTrue(testRelDiff(m1, coaddm1, .01))
        self.assertTrue(testRelDiff(m2, coaddm2, .01))

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
class SimpleGaussianTest(unittest.TestCase):
    """Simple Gaussian Psf test case for CoaddPsf"""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test(self):
        """Check that we can measure a single Gaussian's attributes"""
        sigma0 = 5;
        # set the peak of the outer guassian to 0 so this is really a single gaussian.

        psf = afwDetection.createPsf("DoubleGaussian", 60, 60, 1.5*sigma0, 1, 0.0)

        if False and display:
            im = psf.computeImage(afwGeom.PointD(xwid/2, ywid/2))
            ds9.mtv(im, title="N(%g) psf" % sigma0, frame=0)

        # this is the coadd Wcs we want
        cd11 = 5.55555555e-05
        cd12 = 0.0
        cd21 = 0.0
        cd22 = 5.55555555e-05
        crval1 = 0.0
        crval2 = 0.0
        crpix = afwGeom.PointD(1000, 1000)
        crval = afwCoord.Coord(afwGeom.Point2D(crval1, crval2))
        wcsref = afwImage.makeWcs(crval,crpix,cd11,cd12,cd21,cd22)


        schema = afwTable.ExposureTable.makeMinimalSchema()
        schema.addField("weight", type="D", doc="Coadd weight")
        mycatalog = afwTable.ExposureCatalog(schema)

        sigma = [5,6,7,8]  # 5 pixels is the same as a sigma of 1 arcsec.

        # lay down a simple pattern of four ccds, set in a pattern of 1000 pixels around the center
        offsets = [(1999,1999), (1999,0), (0, 0), (0,1999)]

#       Imagine a ccd in each of positions +-1000 pixels from the center 
        for i in range(4):
            record = mycatalog.getTable().makeRecord()
            psf = afwDetection.createPsf("DoubleGaussian", 100, 100, sigma[i], 1.00, 1.0);
            record.setPsf(psf)
            crpix = afwGeom.PointD(offsets[i][0], offsets[i][1])
            wcs = afwImage.makeWcs(crval,crpix,cd11,cd12,cd21,cd22)

            # print out the coorinates of this supposed 2000x2000 ccd in wcsref coordinates
            beginCoord = wcs.pixelToSky(0,0)
            endCoord = wcs.pixelToSky(2000, 2000)
            #print "\n", "ccd %d corner coords: "%i, beginCoord.getPosition(), endCoord.getPosition(), endCoord.getEpoch()
            beginPix = wcsref.skyToPixel(beginCoord)
            endPix = wcsref.skyToPixel(endCoord)
            #print "\n", "ccd %d corner pixs: "%i, beginPix, endPix, ", in wcsref system"
            record.setWcs(wcs)
            record['weight'] = 1.0
            record['id'] = i
            bbox = afwGeom.Box2I(afwGeom.Point2I(0,0), afwGeom.Extent2I(2000, 2000))
            record.setBBox(bbox)
            mycatalog.append(record) 
            #img = psf.computeImage(afwGeom.Point2D(1000,1000), afwGeom.Extent2I(100,100), True, False)
            #img.writeFits("img%d.fits"%i)
        print "OK, starting gaussian test"
        mypsf = measAlg.CoaddPsf(mycatalog, wcsref) #, 'weight')
        mypsf.setDefaultImageSize(afwGeom.Extent2I(100, 100))
        i,m1,m2 = getCoaddPsfAttributes(mypsf, afwGeom.Point2I(1000,1000))
        #print "At 1000,1000: count = ",i,m1,m2
        i,m1,m2 = getCoaddPsfAttributes(mypsf, afwGeom.Point2I(1001,1001))
        #print "At 1001,1001: count = ",i,m1,m2
        i,m1,m2 = getCoaddPsfAttributes(mypsf, afwGeom.Point2I(1000,1001))
        #print "At 1000,1001: count = ",i,m1,m2
        psfAttrib = measAlg.PsfAttributes(mypsf, afwGeom.Point2I(1000,1001))
        m1    = psfAttrib.computeGaussianWidth(psfAttrib.FIRST_MOMENT)
        m2    = psfAttrib.computeGaussianWidth(psfAttrib.SECOND_MOMENT)
        #print "Coadd itself:",m1,m2
        #self.assertTrue(abs(sigma0 - sigma) <= 1.0e-2)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#   This test checks to be sure that the weights are being applied correctly in doComputeImage
#   Since the 2nd moments are linear in the function value, we can simply weight the moments
#   and be sure that the resulting moments are correct
class WeightTest(unittest.TestCase):
    """Simple Gaussian Psf test case for CoaddPsf with Psf weighting"""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test(self):
        """Check that we can measure a single Gaussian's attributes"""
        sigma0 = 5;
        # set the peak of the outer guassian to 0 so this is really a single gaussian.

        psf = afwDetection.createPsf("DoubleGaussian", 60, 60, 1.5*sigma0, 1, 0.0)

        if False and display:
            im = psf.computeImage(afwGeom.PointD(xwid/2, ywid/2))
            ds9.mtv(im, title="N(%g) psf" % sigma0, frame=0)

        # this is the coadd Wcs we want
        cd11 = 5.55555555e-05
        cd12 = 0.0
        cd21 = 0.0
        cd22 = 5.55555555e-05
        crval1 = 0.0
        crval2 = 0.0
        crpix = afwGeom.PointD(1000, 1000)
        crval = afwCoord.Coord(afwGeom.Point2D(crval1, crval2))
        wcsref = afwImage.makeWcs(crval,crpix,cd11,cd12,cd21,cd22)


        schema = afwTable.ExposureTable.makeMinimalSchema()
        schema.addField("weight", type="D", doc="Coadd weight")
        mycatalog = afwTable.ExposureCatalog(schema)

        sigma = [5,6,7,8]  # 5 pixels is the same as a sigma of 1 arcsec.

        # lay down a simple pattern of four ccds, set in a pattern of 1000 pixels around the center
        offsets = [(1999,1999), (1999,0), (0, 0), (0,1999)]

#       Imagine a ccd in each of positions +-1000 pixels from the center 
        for i in range(4):
            record = mycatalog.getTable().makeRecord()
            psf = afwDetection.createPsf("DoubleGaussian", 100, 100, sigma[i], 1.00, 0.0);
            record.setPsf(psf)
            crpix = afwGeom.PointD(offsets[i][0], offsets[i][1])
            wcs = afwImage.makeWcs(crval,crpix,cd11,cd12,cd21,cd22)

            # print out the coorinates of this supposed 2000x2000 ccd in wcsref coordinates
            record.setWcs(wcs)
            record['weight'] = 1.0 * (i+1)
            record['id'] = i
            bbox = afwGeom.Box2I(afwGeom.Point2I(0,0), afwGeom.Extent2I(2000, 2000))
            record.setBBox(bbox)
            mycatalog.append(record) 
        print "OK, starting weight test"
        mypsf = measAlg.CoaddPsf(mycatalog, wcsref) #, 'weight')
        mypsf.setDefaultImageSize(afwGeom.Extent2I(100, 100))
        m1,m2,comps = getCoaddPsfAttributes(mypsf, afwGeom.Point2I(1000,1000))
        #print "At 1000,1000: count = %d,"%len(comps),m1,m2
        psfAttrib = measAlg.PsfAttributes(mypsf, afwGeom.Point2I(1000,1000))
        m1coadd    = psfAttrib.computeGaussianWidth(psfAttrib.FIRST_MOMENT)
        m2coadd    = psfAttrib.computeGaussianWidth(psfAttrib.SECOND_MOMENT)
        #print "Coadd itself:",m1coadd,m2coadd
        #print comps
        self.assertTrue(testRelDiff(m1,m1coadd,.01))
        m1,m2,comps = getCoaddPsfAttributes(mypsf, afwGeom.Point2I(1001,1001))
        #print "At 1001,1001: count = %d"%len(comps),m1,m2
        psfAttrib = measAlg.PsfAttributes(mypsf, afwGeom.Point2I(1001,1001))
        m1coadd    = psfAttrib.computeGaussianWidth(psfAttrib.FIRST_MOMENT)
        m2coadd    = psfAttrib.computeGaussianWidth(psfAttrib.SECOND_MOMENT)
        #print "Coadd itself:",m1coadd,m2coadd
        #print comps
        self.assertTrue(testRelDiff(m1,m1coadd,.01))
        m1,m2,comps = getCoaddPsfAttributes(mypsf, afwGeom.Point2I(1001,1000))
        #print "At 1001,1001: count = %d"%len(comps),m1,m2
        psfAttrib = measAlg.PsfAttributes(mypsf, afwGeom.Point2I(1001,1000))
        m1coadd    = psfAttrib.computeGaussianWidth(psfAttrib.FIRST_MOMENT)
        m2coadd    = psfAttrib.computeGaussianWidth(psfAttrib.SECOND_MOMENT)
        #print "Coadd itself:",m1coadd,m2coadd
        #print comps
        self.assertTrue(testRelDiff(m1,m1coadd,.01))
        m1,m2,comps = getCoaddPsfAttributes(mypsf, afwGeom.Point2I(1000,1001))
        #print "At 1000,1001: count = %d"%len(comps),m1,m2
        psfAttrib = measAlg.PsfAttributes(mypsf, afwGeom.Point2I(1000,1001))
        m1coadd    = psfAttrib.computeGaussianWidth(psfAttrib.FIRST_MOMENT)
        m2coadd    = psfAttrib.computeGaussianWidth(psfAttrib.SECOND_MOMENT)
        #print "Coadd itself:",m1coadd,m2coadd
        self.assertTrue(testRelDiff(m1,m1coadd,.01))
        #print comps

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(CreatePsfTest)
    #suites += unittest.makeSuite(RotatePsfTest)
    suites += unittest.makeSuite(SimpleGaussianTest)
    suites += unittest.makeSuite(WeightTest)
    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the utilsTests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
