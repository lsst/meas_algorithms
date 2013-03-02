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
import unittest

import numpy
import eups

import lsst.daf.base as dafBase
import lsst.afw.math            as afwMath
import lsst.pex.exceptions      as pexEx
import lsst.pex.policy          as policy
import lsst.pex.config          as pexConf
import lsst.pex.logging         as pexLog
import lsst.afw.image           as afwImage
import lsst.afw.detection       as afwDet
import lsst.afw.geom            as afwGeom
import lsst.afw.table           as afwTable
import lsst.meas.algorithms     as measAlg

import lsst.utils.tests         as utilsTests

import lsst.afw.display.ds9       as ds9

try:
    type(verbose)
except NameError:
    verbose = 0

def plantSources(bbox, kwid, sky, coordList, addPoissonNoise=True):
    """Make an exposure with stars (modelled as Gaussians)
    
    @param bbox: parent bbox of exposure
    @param kwid: kernel width (and height; kernel is square)
    @param sky: amount of sky background (counts)
    @param coordList: a list of [x, y, counts, sigma], where:
        * x,y are relative to exposure origin
        * counts is the integrated counts for the star
        * sigma is the Gaussian sigma in pixels
    @param addPoissonNoise: add Poisson noise to the exposure?
    """
    # make an image with sources
    img = afwImage.ImageD(bbox)
    meanSigma = 0.0
    for coord in coordList:
        x, y, counts, sigma = coord
        meanSigma += sigma

        # make a single gaussian psf
        psf = afwDet.createPsf("SingleGaussian", kwid, kwid, sigma)

        # make an image of it and scale to the desired number of counts
        normPeak = False
        thisPsfImg = psf.computeImage(afwGeom.PointD(int(x), int(y)), normPeak)
        thisPsfImg *= counts

        # bbox a window in our image and add the fake star image
        imgSeg = img.Factory(img, thisPsfImg.getBBox(afwImage.PARENT), afwImage.PARENT)
        imgSeg += thisPsfImg
    meanSigma /= len(coordList)

    img += sky

    # add Poisson noise
    if (addPoissonNoise):
        numpy.random.seed(seed=1) # make results reproducible
        imgArr = img.getArray()
        imgArr[:] = numpy.random.poisson(imgArr)

    # bundle into a maskedimage and an exposure
    mask = afwImage.MaskU(bbox)
    var = img.convertFloat()
    img -= sky
    mimg = afwImage.MaskedImageF(img.convertFloat(), mask, var)
    exposure = afwImage.makeExposure(mimg)

    # insert an approximate psf
    psf = afwDet.createPsf("SingleGaussian", kwid, kwid, meanSigma)
    exposure.setPsf(psf)

    return exposure


    
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
display = False
class ApertureCorrectionTestCase(unittest.TestCase):
    """Test the aperture correction."""

    def setUp(self):
        self.nx, self.ny = 128, 128
        self.bbox = afwGeom.Box2I(afwGeom.Point2I(0,0), afwGeom.Extent2I(self.nx, self.ny))
        self.ngrid        = 5
        self.sigma0      = 1.5
        self.val         = 40000.0
        self.sky         = 100.0
        self.kwid        = int(self.sigma0*7)
        if not self.kwid%2: self.kwid += 1

        # enable/disable the Assert statements
        # - diabling allows all tests in a method to run and print output
        # rather than throw an exception before remaining tests can run
        self.doTest      = True #False

        # how big do we allow the error to be
        #self.maxErrorFrac = 0.005       # half a percent fine for NAIVE
        self.maxErrorFrac = 0.003        # 0.3 percent needed for SINC
        # how many sigma can the measured value be from theoretical?
        # note: we're checking all candidate stars
        #self.nSigmaErrorLimit = 1.1     # 1.1 good for NAIVE
        self.nSigmaErrorLimit = 7.0      # 3.2 needed for SINC to past constant psf test

        # Note: SINC is not exactly as expected because of tapering of the aperture.
        
        self.metadata = dafBase.PropertyList()

        # detection configuration
        self.detConfig = measAlg.SourceDetectionConfig()
        self.detConfig.reEstimateBackground = False

        # measurement configuration
        self.measConfig = measAlg.SourceMeasurementConfig()
        
        # psf algorithms and policies
        # apcorr policies
        self.apCorrConfig = measAlg.ApertureCorrection.ConfigClass()
        self.apCorrConfig.polyStyle = "standard" # this does better than cheby ??
        self.apCorrConfig.order     = 2
        self.apCorrConfig.alg1.name = "flux.psf"
        self.apCorrConfig.alg2.name = "flux.sinc"
        self.apCorrConfig.alg2["flux.sinc"].radius = 3.0
        self.alg1 = self.apCorrConfig.alg1
        self.alg2 = self.apCorrConfig.alg2

        # logs
        self.log = pexLog.getDefaultLog()
        self.log.setThreshold(self.log.WARN)

        self.nDisp = 1
        
    def tearDown(self):
        del self.bbox
        del self.detConfig
        del self.measConfig
        del self.apCorrConfig
        del self.metadata
        del self.log
        pass
    
    #################################################################
    # quick and dirty detection (note: we already subtracted background)
    def detectAndMeasure(self, exposure):
        schema = afwTable.SourceTable.makeMinimalSchema()
        self.detConfig.validate()
        self.measConfig.validate()
        detTask = measAlg.SourceDetectionTask(config=self.detConfig, schema=schema)
        measTask = measAlg.SourceMeasurementTask(config=self.measConfig, schema=schema)
        # detect
        table = afwTable.SourceTable.make(schema)
        sources = detTask.makeSourceCatalog(table, exposure).sources
        # ... and measure
        measTask.run(exposure, sources)
        return sources

    ###################################################
    # Compute the theoretical fraction of flux inside
    #   a radius, r, for a Gaussian.
    # Solution in 2D is analytic: integral of r*exp(-r**2)
    #   is returned (solve by parts)
    ###################################################
    def apCorrTheory(self, sigma, r):
        return 1.0 - math.exp(-r**2/(2.0*sigma**2))
        

    ###################################################
    # Compute the flux+err expected for the different types
    #   of photometry used: flux.psf, flux.sinc, flux.naive
    ###################################################
    def getKnownFluxes(self, psfImg, radius, counts, sigma):

        flux = {"flux.psf": 0.0, "flux.sinc": 0.0, "flux.naive": 0.0 }
        fluxErr = {"flux.psf": 0.0, "flux.sinc": 0.0, "flux.naive": 0.0 }
        measErr = {"flux.psf": 0.0, "flux.sinc": 0.0, "flux.naive": 0.0 }

        xw, yw = psfImg.getWidth(), psfImg.getHeight()
        x0, y0 = psfImg.getX0(), psfImg.getY0()
        ix, iy = xw/2, yw/2
        #if not xw % 2: ix -= 1
        #if not yw % 2: iy -= 1
            
        psfSum, psfSumSqrd = 0.0, 0.0
        for j in range(xw):
            for i in range(yw):
                w = psfImg.get(i, j)
                psfSum += w
                psfSumSqrd += w*w
                f = w*counts

                # add up the psf flux
                fluxErr["flux.psf"] += w*w*f
                flux["flux.psf"] += w*f

                # add up the naive fluxes
                dx, dy = i-ix, j-iy
                if (dx*dx + dy*dy <= radius*radius):
                    flux["flux.naive"] += f
                    fluxErr["flux.naive"] += f
                    # use smallest value as error ... ad-hoc
                    if f < measErr["flux.naive"] or measErr["flux.naive"] == 0:
                        measErr["flux.naive"] = f

        # renormalize the psf fluxes
        flux["flux.psf"] *= psfSum/psfSumSqrd
        fluxErr["flux.psf"] = math.sqrt(fluxErr["flux.psf"])*psfSum/psfSumSqrd
        measErr["flux.psf"] = 0.0 #fluxErr["flux.psf"]/math.sqrt(psfSum)
        
        fluxErr["flux.naive"] = math.sqrt(fluxErr["flux.naive"])
        measErr["flux.naive"] = math.sqrt(measErr["flux.naive"])

        # use the analytic form for the integral of a single gaussian for the sinc
        # - it's not quite right because of the cos tapering
        frac = self.apCorrTheory(sigma, self.alg2.active.radius)
        flux["flux.sinc"] = counts*frac
        fluxErr["flux.sinc"] = math.sqrt(flux["flux.sinc"])
        measErr["flux.sinc"] = 0.0

        return flux, fluxErr, measErr


    #######################################################
    # Get the known aperture correction based on the values
    #   from getKnownFluxes()
    # Errors are propegrated, but theoretical errors are
    #   ill defined for the knownFluxes
    #######################################################
    def getKnownApCorr(self, fluxKnown, fluxKnownErr, measKnownErr):
        apCorr    = fluxKnown[self.alg2.name]/fluxKnown[self.alg1.name]
        apCorrErr = apCorr*(measKnownErr[self.alg1.name]/fluxKnown[self.alg1.name] +
                            measKnownErr[self.alg2.name]/fluxKnown[self.alg2.name])
        return apCorr, apCorrErr


    ########################################################
    # print a summary of what we measured
    ########################################################
    def printSummary(self, psfImg, fluxKnown, fluxKnownErr, measKnownErr, ac):
    
        # print diagnostics on the star selection
        print "Used %d apCorr stars (%d good)" % (self.metadata.get("numAvailStars"),
                                                  self.metadata.get("numGoodStars"))
        
        # have a look at the know values
        print "Flux known (%s): %.2f +/- %.2f" % (self.alg1.name, fluxKnown[self.alg1.name],
                                                  fluxKnownErr[self.alg1.name])
        print "Flux known (%s): %.2f +/- %.2f" % (self.alg2.name, fluxKnown[self.alg2.name],
                                                  fluxKnownErr[self.alg2.name])
        apCorr, apCorrErr    = self.getKnownApCorr(fluxKnown, fluxKnownErr, measKnownErr)
        print "Aperture Corr'n Known: %.4f +/- %.4f" % (apCorr, apCorrErr)
        apcorr, apcorrErr = ac.computeAt(self.nx/2, self.ny/2)
        print "Aperture Corr'n meas: %.4f +/- %.4f" % (apcorr, apcorrErr)

            

    #########################################################
    # The main workhorse of code
    # - plant a list of given objects
    # - get a psf
    # - get the apCorr
    # - test the results
    #########################################################
    def plantAndTest(self, coordList):

        # plant them in the image, and measure them
        exposure   = plantSources(self.bbox, self.kwid, self.sky, coordList)
        sourceList = self.detectAndMeasure(exposure)
        mimg = exposure.getMaskedImage()
        img = mimg.getImage()

        if display:
            ds9.mtv(img,      frame=self.nDisp, title="Delta functions")
            self.nDisp += 1
            ds9.mtv(mimg,     frame=self.nDisp, title="convolved image")
            self.nDisp += 1
        
        # try getPsf()
        starSelectorFactory = measAlg.starSelectorRegistry["secondMoment"]
        starSelectorConfig = starSelectorFactory.ConfigClass()
        starSelectorConfig.badFlags = []
        starSelector = starSelectorFactory(starSelectorConfig)
        
        psfDeterminerFactory = measAlg.psfDeterminerRegistry["pca"]
        psfDeterminerConfig = psfDeterminerFactory.ConfigClass()
        psfDeterminerConfig.sizeCellX = self.nx/4
        psfDeterminerConfig.sizeCellY = self.ny/4
        psfDeterminer = psfDeterminerFactory(psfDeterminerConfig)
        
        psfCandidateList = starSelector.selectStars(exposure, sourceList)
        psf, cellSet = psfDeterminer.determinePsf(exposure, psfCandidateList, self.metadata)
        exposure.setPsf(psf)

        ##########################################
        # try the aperture correction
        self.log.setThreshold(self.log.INFO)
        ac = measAlg.ApertureCorrection(exposure, cellSet, self.metadata, self.apCorrConfig, log=self.log)
        
        if display:

            # show the apCorr and error as images
            acImg = afwImage.ImageF(self.nx, self.ny)
            acErrImg = afwImage.ImageF(self.nx, self.ny)
            for j in range(self.ny):
                for i in range(self.nx):
                    apCo, apCoErr = ac.computeAt(i, j)
                    acImg.set(i, j, apCo)
                    acErrImg.set(i, j, apCoErr)

            ds9.mtv(acImg,    frame=self.nDisp, title="Apcorr Image")
            self.nDisp += 1            
            ds9.mtv(acErrImg, frame=self.nDisp, title="Apcorr Error Image")
            self.nDisp += 1
        

        # print info for the middle object
        xmid, ymid, valid, sigmid = coordList[len(coordList)/2]
        normPeak = False
        psfImg = psf.computeImage(afwGeom.PointD(int(xmid), int(ymid)), normPeak)
        fluxKnown, fluxKnownErr, measKnownErr = self.getKnownFluxes(psfImg, self.alg2.active.radius, 
                                                                    self.val, sigmid)
        self.printSummary(psfImg, fluxKnown, fluxKnownErr, measKnownErr, ac)

        if display:
            ds9.mtv(psfImg,   frame=self.nDisp, title="Psf Image")
            self.nDisp += 1            

            
        ############################################
        # for each thing we planted ... check it
        iCoord = -1
        everyNth = 2
        for coord in coordList:
            iCoord += 1

            # ok ... not *every* planted object
            if iCoord % everyNth:
                continue
            
            x, y, val, sigma = coord
        
            normPeak = False
            psfImg = psf.computeImage(afwGeom.PointD(int(x), int(y)), normPeak)
            fluxKnown, fluxKnownErr, measKnownErr = self.getKnownFluxes(psfImg, self.alg2.active.radius,
                                                                        self.val, sigma)

            corrKnown, corrErrKnown           = self.getKnownApCorr(fluxKnown, fluxKnownErr, measKnownErr)
            corrMeasMiddle, corrErrMeasMiddle = ac.computeAt(x, y)
            
            print "%3d %3d %5.3f %6.4f %6.4f  %5.3f" % (x, y, sigma, corrMeasMiddle, corrKnown,
                                                       corrMeasMiddle/corrKnown),

            
            ###################
            # Tests
            ###################
            
            # verify we're within error (1 stdev)
            discrep = abs(corrKnown - corrMeasMiddle)
            error = self.nSigmaErrorLimit*(corrErrMeasMiddle)   # ie. +/-  ~nSigErrLim*sigma
            print "discrep: %6.4f %6.4f" % (discrep, error),
            if (discrep < error):
                print "pass",
            else:
                print "FAIL",
            if self.doTest:
                self.assertTrue(discrep < error)

            # and that error is small
            errFrac = corrErrMeasMiddle/corrMeasMiddle
            print "errFrac: %5.3f" % (errFrac),
            if (errFrac < self.maxErrorFrac):
                print "pass"
            else:
                print "FAIL"
            if self.doTest:
                self.assertTrue(errFrac < self.maxErrorFrac)

                

    #####################################################
    # Test for Constant Psf
    #####################################################
    def testApCorrConstantPsf(self):
        """ Verify that we can recover the known aperture correction for a *constant* psf."""

        dx = self.nx/(self.ngrid + 1)
        dy = self.ny/(self.ngrid + 1)

        # decide where to put fake psfs on a grid
        coordList = []
        for i in range(self.ngrid):
            for j in range(self.ngrid):
                x, y = (1+i)*dx, (1+j)*dy
                coordList.append([x, y, self.val, self.sigma0])

        self.plantAndTest(coordList)
        
            
    #####################################################
    # Test for Linearly varying aperture correction
    #####################################################
    def testApCorrLinearVaryingPsf(self):
        """Verify that we can recover the known aperture correction for *linearly varying* psf."""

        dx = self.nx/(self.ngrid + 1)
        dy = self.ny/(self.ngrid + 1)

        # vary apCorr by dApCorr linearly across the image
        apCorr = self.apCorrTheory(self.sigma0, self.alg2.active.radius)
        # want aperture correction to vary by this much across the field ...
        dApCorr = 0.05*apCorr
        # ... and that means changing sigma by this much:
        # (integrate r*exp(-r**2/sigma**2), and solve sigma)
        sig2   = self.alg2.active.radius*(-2.0*math.log(1.0 - (apCorr+dApCorr)))**-0.5

        # deriv to scale sigma by
        dsigmaDx   = (sig2 - self.sigma0)/self.nx
        
        # decide where to put fake psfs on a grid
        coordList = []
        for i in range(self.ngrid):
            for j in range(self.ngrid):
                x, y = (1+i)*dx, (1+j)*dy
                coordList.append([x, y, self.val, self.sigma0+dsigmaDx*x])

        self.plantAndTest(coordList)


    #####################################################
    # Test for Quadratic varying aperture correction
    #####################################################
    def testApCorrQuadraticVaryingPsf(self):
        """Verify that we can recover the known aperture correction for *quadraticly varying* psf."""

        dx = self.nx/(self.ngrid + 1)
        dy = self.ny/(self.ngrid + 1)
        xmid, ymid = self.nx/2, self.ny/2
        
        # vary apCorr by dApCorr quadratically across the image
        apCorr = self.apCorrTheory(self.sigma0, self.alg2.active.radius)
        # want aperture correction to vary by this much across the field ...
        dApCorr = -0.05*apCorr
        # ... and that means changing sigma by this much:
        # (integrate r*exp(-r**2/sigma**2), and solve sigma)
        sig2   = self.alg2.active.radius*(-2.0*math.log(1.0 - (apCorr+dApCorr)))**-0.5

        # define 2nd derivs for a parabola
        dsigmaDx2   = (sig2 - self.sigma0)/((0.5*self.nx)**2)
        dsigmaDy2   = (sig2 - self.sigma0)/((0.5*self.ny)**2)
        
        # decide where to put fake psfs on a grid
        coordList = []
        for i in range(self.ngrid):
            for j in range(self.ngrid):
                x, y = (1+i)*dx, (1+j)*dy
                # center the parabola in the middle of the frame
                xp, yp = x-xmid, y-ymid 
                coordList.append([x, y, self.val, self.sigma0 + dsigmaDx2*xp*xp + dsigmaDy2*yp*yp])

        self.plantAndTest(coordList)
        
        
        
        
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(ApertureCorrectionTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)
 
if __name__ == "__main__":
    run(True)


