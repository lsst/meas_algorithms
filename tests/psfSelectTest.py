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

# todo:
# - growth curves
# - 

import math
import pdb                          # we may want to say pdb.set_trace()
import unittest

import numpy
import eups

import lsst.daf.base            as dafBase
import lsst.afw.math            as afwMath
import lsst.pex.exceptions      as pexEx
import lsst.pex.policy          as pexPolicy
import lsst.pex.logging         as pexLog
import lsst.afw.image           as afwImage
import lsst.afw.detection       as afwDet
import lsst.afw.geom            as afwGeom
import lsst.afw.geom.ellipses   as geomEllip
import lsst.meas.algorithms     as measAlg

import lsst.afw.cameraGeom      as cameraGeom
import lsst.afw.cameraGeom.utils as cameraUtils

import lsst.utils.tests         as utilsTests
import lsst.sdqa                as sdqa

import lsst.afw.display.ds9     as ds9

import sourceDetectionBickTmp   as srcDet
import sourceMeasurementBickTmp as srcMeas

try:
    type(verbose)
except NameError:
    verbose = 0

try:
    display
except NameError:
    display = False


def plantSources(x0, y0, nx, ny, sky, nObj, wid, distorter, useRandom=False):

    img0 = afwImage.ImageF(afwGeom.ExtentI(nx, ny))
    img = afwImage.ImageF(afwGeom.ExtentI(nx, ny))
    
    ixx0, iyy0, ixy0 = wid*wid, wid*wid, 0.0
    m0 = geomEllip.Quadrupole(ixx0, iyy0, ixy0)

    edgeBuffer = 40.0*wid
    
    flux = 1.0e4
    nkx, nky = int(10*wid) + 1, int(10*wid) + 1
    xhwid,yhwid = nkx/2, nky/2

    nRow = int(math.sqrt(nObj))
    xstep = (nx - 0*edgeBuffer)/(nRow+1)
    ystep = (ny - 0*edgeBuffer)/(nRow+1)

    if useRandom:
	nObj = nRow*nRow

    goodAdded0 = 0
    goodAdded = 0
    
    for i in range(nObj):

	# get our position
	if useRandom:
	    xcen0, ycen0 = numpy.random.uniform(nx), numpy.random.uniform(ny)
	else:
	    xcen0, ycen0 = xstep*((i%nRow) + 1), ystep*(int(i/nRow) + 1)
	ixcen0, iycen0 = int(xcen0), int(ycen0)

	# distort position and shape
	p = distorter.distort(afwGeom.Point2D(xcen0, ycen0))
	m = distorter.distort(afwGeom.Point2D(x0+xcen0, y0+ycen0), m0)
	xcen, ycen = xcen0, ycen0 #p.getX(), p.getY()
	if (xcen < 1.5*edgeBuffer or (nx - xcen) < 1.5*edgeBuffer or
	    ycen < 1.5*edgeBuffer or (ny - ycen) < 1.5*edgeBuffer):
	    continue
	ixcen, iycen = int(xcen), int(ycen)
	ixx, iyy, ixy = m.getIXX(), m.getIYY(), m.getIXY()

	# plant the object
	tmp = 0.25*(ixx-iyy)**2 + ixy**2
	a2 = 0.5*(ixx+iyy) + numpy.sqrt(tmp)
	b2 = 0.5*(ixx+iyy) - numpy.sqrt(tmp)
	#ellip = 1.0 - numpy.sqrt(b2/a2)
	theta = 0.5*numpy.arctan2(2.0*ixy, ixx-iyy)
	a = numpy.sqrt(a2)
	b = numpy.sqrt(b2)

	c, s = math.cos(theta), math.sin(theta)
	good0, good = True, True
	for y in range(nky):
	    iy = iycen + y - yhwid
	    iy0 = iycen0 + y - yhwid
	    
	    for x in range(nkx):
		ix = ixcen + x - xhwid
		ix0 = ixcen0 + x - xhwid
		
		if ix >= 0 and ix < nx and iy >= 0 and iy < ny:
		    dx, dy = ix - xcen, iy - ycen
		    u =  c*dx + s*dy
		    v = -s*dx + c*dy
		    I0 = flux/(2*math.pi*a*b)		
		    val = I0*math.exp(-0.5*((u/a)**2 + (v/b)**2))
		    if val < 0:
			val = 0
		    prevVal = img.get(ix, iy)
		    img.set(ix, iy, val+prevVal)
		else:
		    good = False

		if ix0 >=0 and ix0 < nx and iy0 >= 0 and iy0 < ny:
		    dx, dy = ix - xcen, iy - ycen
		    I0 = flux/(2*math.pi*wid*wid)		
		    val = I0*math.exp(-0.5*((dx/wid)**2 + (dy/wid)**2))
		    if val < 0:
			val = 0
		    prevVal = img0.get(ix0, iy0)
		    img0.set(ix0, iy0, val+prevVal)
		else:
		    good0 = False
		    
	goodAdded0 += 1 if good0 else 0
	goodAdded += 1 if good else 0

    # add sky and noise
    img += sky
    img0 += sky
    noise = afwImage.ImageF(afwGeom.ExtentI(nx, ny))
    noise0 = afwImage.ImageF(afwGeom.ExtentI(nx, ny))
    for i in range(nx):
	for j in range(ny):
	    noise.set(i, j, numpy.random.poisson(img.get(i,j) ))
	    noise0.set(i, j, numpy.random.poisson(img0.get(i,j) ))


    edgeWidth = int(0.5*edgeBuffer)
    mask = afwImage.MaskU(afwGeom.ExtentI(nx, ny))
    left   = afwGeom.Box2I(afwGeom.Point2I(0,0), afwGeom.ExtentI(edgeWidth, ny))
    right  = afwGeom.Box2I(afwGeom.Point2I(nx - edgeWidth,0), afwGeom.ExtentI(edgeWidth, ny))
    top    = afwGeom.Box2I(afwGeom.Point2I(0,ny - edgeWidth), afwGeom.ExtentI(nx, edgeWidth))
    bottom = afwGeom.Box2I(afwGeom.Point2I(0,0), afwGeom.ExtentI(nx, edgeWidth))
    
    for pos in [left, right, top, bottom]:
	msk = afwImage.MaskU(mask, pos, False)
	msk.set(msk.getPlaneBitMask('EDGE'))
	
    expos = afwImage.makeExposure(afwImage.makeMaskedImage(noise, mask, afwImage.ImageF(noise, True)))
    expos0 = afwImage.makeExposure(afwImage.makeMaskedImage(noise0, mask, afwImage.ImageF(noise0, True)))

    im = expos.getMaskedImage().getImage()
    im0 = expos0.getMaskedImage().getImage()
    im -= sky
    im0 -= sky

    
    return expos, goodAdded, expos0, goodAdded0


#################################################################
# quick and dirty detection (note: we already subtracted background)
def detectAndMeasure(exposure, detPolicy, measPolicy):
    
    # detect
    dsPos, dsNeg   = srcDet.detectSources(exposure, exposure.getPsf(), detPolicy)
    footprintLists = [[dsPos.getFootprints(),[]]]
    # ... and measure
    sourceList     = srcMeas.sourceMeasurement(exposure, exposure.getPsf(),
                                                   footprintLists, measPolicy)
    return sourceList

    
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
display = False
class PsfSelectionTestCase(unittest.TestCase):
    """Test the aperture correction."""

    def setUp(self):
	self.x0, self.y0 = 0, 0
        self.nx, self.ny = 512, 512 #2048, 4096
	self.sky = 100.0
	self.nObj = 100

	self.sCamCoeffs = [0.0, 1.0, 7.16417e-08, 3.03146e-10, 5.69338e-14, -6.61572e-18]
	self.sCamDistorter = cameraGeom.RadialPolyDistortion(self.sCamCoeffs)
	
	self.sCamCoeffsExag = [0.0, 1.0, 3.0e-04, 3.03146e-7, 5.69338e-11, -6.61572e-15]
	self.sCamDistorterExag = cameraGeom.RadialPolyDistortion(self.sCamCoeffsExag)

        self.detPolicy = pexPolicy.Policy.createPolicy(pexPolicy.DefaultPolicyFile("meas_algorithms",
                                                                             "detectionDictionaryBickTmp.paf",
                                                                             "tests"))

        self.measPolicy = pexPolicy.Policy.createPolicy(pexPolicy.DefaultPolicyFile("meas_algorithms",
									      "MeasureSourcesDictionary.paf",
									      "policy"))
        self.secondMomentStarSelectorPolicy = pexPolicy.Policy.createPolicy(
            pexPolicy.DefaultPolicyFile("meas_algorithms", "policy/secondMomentStarSelectorDictionary.paf"))
	self.secondMomentStarSelectorPolicy.set('fluxLim', 5000.0)

	self.starSelector = measAlg.makeStarSelector("secondMomentStarSelector",
						     self.secondMomentStarSelectorPolicy)



        pcaPsfDeterminerPolicy = pexPolicy.Policy.createPolicy(
            pexPolicy.DefaultPolicyFile("meas_algorithms", "policy/pcaPsfDeterminerDictionary.paf"))
        width, height = self.nx, self.ny
	nEigenComponents = 3
        pcaPsfDeterminerPolicy.set("sizeCellX", width//3)
        pcaPsfDeterminerPolicy.set("sizeCellY", height//3)
        pcaPsfDeterminerPolicy.set("nEigenComponents", nEigenComponents)
        pcaPsfDeterminerPolicy.set("spatialOrder", 1)
        pcaPsfDeterminerPolicy.set("kernelSizeMin", 31)
        pcaPsfDeterminerPolicy.set("nStarPerCell", 0)
        pcaPsfDeterminerPolicy.set("nStarPerCellSpatialFit", 0) # unlimited

        self.psfDeterminer = measAlg.makePsfDeterminer("pcaPsfDeterminer", pcaPsfDeterminerPolicy)


        rawPolicy = pexPolicy.Policy.createPolicy(pexPolicy.DefaultPolicyFile("afw",
                                                                              "singleCcdGeom.paf",
                                                                              "tests"))
        self.cameraPolicy = cameraUtils.getGeomPolicy(rawPolicy)

        
    def tearDown(self):
	del self.detPolicy
	del self.measPolicy
	del self.secondMomentStarSelectorPolicy
	del self.sCamDistorter
	del self.sCamDistorterExag
	del self.starSelector
	del self.psfDeterminer
        del self.cameraPolicy


    def testPsfCandidate(self):

	distorter = self.sCamDistorterExag

	# make an exposure
	psfSigma = 1.5
	exposDist, nGoodDist, expos0, nGood0 = plantSources(self.x0, self.y0,
							    self.nx, self.ny,
							    self.sky, self.nObj, psfSigma, distorter)

	# set the psf
	kwid = 21 #int(10*psfSigma) + 1
	psf = afwDet.createPsf("SingleGaussian", kwid, kwid, psfSigma)
	exposDist.setPsf(psf)
	
	# Distorter lives in Detector in an Exposure
	camera = cameraUtils.makeCamera(self.cameraPolicy)
        raft = cameraGeom.cast_Raft(camera[0])
        detector = cameraGeom.cast_Ccd(raft[0])
	
	exposDist.setDetector(detector)

	detector.setDistortion(distorter)
        print "detect"
	sourceList       = detectAndMeasure(exposDist, self.detPolicy, self.measPolicy)

        metadata = dafBase.PropertyList()
        print "select"
        psfCandidateList = self.starSelector.selectStars(exposDist, sourceList)
        print "determine"
        psf, cellSet = self.psfDeterminer.determinePsf(exposDist, psfCandidateList, metadata)

        
	settings = {'scale': 'minmax', 'zoom':"to fit", 'mask':'transparency 80'}
	ds9.mtv(exposDist, frame=1, title="full", settings=settings)
	for s in sourceList:
	    x, y = s.getXAstrom(), s.getYAstrom()
            
	    if x < 400 or y < 400:
		continue

	    cand = measAlg.makePsfCandidate(s, exposDist)
	    img = cand.getImage()	    
	    uimg = cand.getUndistImage(img.getWidth(), img.getHeight())

	    buffer = 3
	    warpAlg = "lanczos3"

	    oimg = cand.getOffsetImage(warpAlg, buffer)
	    uoimg = cand.getUndistOffsetImage(warpAlg, buffer)

	    ds9.mtv(img, frame=2, title="image", settings=settings)
	    ds9.mtv(oimg, frame=3, title="offset", settings=settings)
	    ds9.mtv(uimg, frame=4, title="undistorted", settings=settings)
	    ds9.mtv(uoimg, frame=5, title="undistorted offset", settings=settings)

	    psfImg = psf.computeImage(afwGeom.Point2D(x,y), True, False)
	    ds9.mtv(psfImg, frame=6, title="psf", settings=settings)
            
	    psfImgDist = psf.computeImage(afwGeom.Point2D(x, y))
	    ds9.mtv(psfImgDist, frame=7, title="psfDist", settings=settings)


    def testDistortedImage(self):

	distorter = self.sCamDistorterExag
	
	psfSigma = 1.5
	exposDist, nGoodDist, expos0, nGood0 = plantSources(self.x0, self.y0,
							    self.nx, self.ny,
							    self.sky, self.nObj, psfSigma, distorter)
	
	kwid = int(10*psfSigma) + 1
	psf = afwDet.createPsf("SingleGaussian", kwid, kwid, psfSigma)
	exposDist.setPsf(psf)
	expos0.setPsf(psf)

	expos = exposDist
	nGood = nGoodDist

	# Distorter lives in Detector in an Exposure
	camera = cameraUtils.makeCamera(self.cameraPolicy)
        raft = cameraGeom.cast_Raft(camera[0])
        detector = cameraGeom.cast_Ccd(raft[0])

	expos.setDetector(detector)

	########################
	# try without distorter
	detector.setDistortion(None)
	sourceList       = detectAndMeasure(expos, self.detPolicy, self.measPolicy)
	psfCandidateListx = self.starSelector.selectStars(expos, sourceList)

	########################
	# try with distorter
	detector.setDistortion(distorter)
	sourceList       = detectAndMeasure(expos, self.detPolicy, self.measPolicy)
	psfCandidateList = self.starSelector.selectStars(expos, sourceList)

	print "uncorrected nAdded,nCand: ", len(psfCandidateListx), nGood
	print "dist-corrected nAdded,nCand: ", len(psfCandidateList), nGood
	
        sys.exit()
	########################
	# display
	if True: #display:
	    iDisp = 1
	    ds9.mtv(exposDist, frame=iDisp);	    iDisp += 1
	    ds9.mtv(expos0, frame=iDisp);	    iDisp += 1
	    ds9.mtv(expos, frame=iDisp)
	    for c in psfCandidateList:
		s = c.getSource()
		ds9.dot("o", s.getXAstrom(), s.getYAstrom(), frame=iDisp)

	# we shouldn't expect to get all available stars
	self.assertTrue(len(psfCandidateListx) < nGood)
	# we should get all of them
	self.assertEqual(len(psfCandidateList), nGood)

		
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(PsfSelectionTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)
 
if __name__ == "__main__":
    run(True)


