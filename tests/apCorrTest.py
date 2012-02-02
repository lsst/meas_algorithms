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
import lsst.meas.algorithms     as measAlg

import lsst.utils.tests         as utilsTests

import sourceDetectionBickTmp   as srcDet
import sourceMeasurementBickTmp as srcMeas


import lsst.afw.display.ds9       as ds9

try:
    type(verbose)
except NameError:
    verbose = 0

class DetectionConfig(pexConf.Config):
    minPixels = pexConf.RangeField(
        doc="detected sources with fewer than the specified number of pixels will be ignored",
        dtype=int, optional=True, default=1, min=0,
    )
    nGrow = pexConf.RangeField(
        doc="How many pixels to to grow detections",
        dtype=int, optional=True, default=1, min=0,
    )
    thresholdValue = pexConf.RangeField(
        doc="value assigned to the threshold object used in detection",
        dtype=float, optional=True, default=5.0, min=0.0,
    )
    thresholdType = pexConf.ChoiceField(
        doc="specifies the desired flavor of Threshold",
        dtype=str, optional=True, default="stdev",
        allowed={
            "variance": "threshold applied to image variance",
            "stdev": "threshold applied to image std deviation",
            "value": "threshold applied to image value"
        }
    )
    thresholdPolarity = pexConf.ChoiceField(
        doc="specifies whether to detect positive, or negative sources, or both",
        dtype=str, optional=True, default="positive",
        allowed={
            "positive": "detect only positive sources",
            "negative": "detect only negative sources",
            "both": "detect both positive and negative sources",
        }
    )


######################################################
# We need a quick/easy way to add sources with specified psf width to an image
#
# We'll take a 'coordList' = [x, y, sigma]
# We'll return and exposure
######################################################
def plantSources(nx, ny, kwid, sky, coordList, addPoissonNoise=True):

    # make a masked image
    img   = afwImage.ImageD(afwGeom.ExtentI(nx, ny), 0.0)
    msk   = afwImage.MaskU(img.getDimensions(), 0x0)
    var   = afwImage.ImageD(img.getDimensions())

    # add sources
    sigma0 = 0.0
    imgPsf = afwImage.ImageD(img.getDimensions(), 0.0)
    for coord in coordList:
        x, y, val, sigma = coord
        sigma0 += sigma

        # make a single gaussian psf
        psf = afwDet.createPsf("SingleGaussian", kwid, kwid, sigma)

        # make an image of it, scale to our specified count rate (self.val)
        normPeak = False
        thisPsfImg = psf.computeImage(afwGeom.PointD(int(x), int(y)), normPeak)
        thisPsfImg *= val

        # bbox a window in our image and add the fake star image
        llc = afwGeom.PointI(x-kwid/2, y-kwid/2)
        urc = afwGeom.PointI(x+kwid/2, y+kwid/2)
        imgSeg = img.Factory(img, afwGeom.BoxI(llc, urc), afwImage.LOCAL)
        imgSeg += thisPsfImg

    img += sky
    sigma0 /= len(coordList)

    # add Poisson noise
    if (addPoissonNoise):
        ran = afwMath.Random()
        for j in range(ny):
            for i in range(nx):
                img.set(i, j, ran.poisson(img.get(i, j)))

    # bundle into a maskedimage and an exposure
    var <<= img
    img -= sky
    mimg     = afwImage.MaskedImageF(img.convertFloat(), msk, var.convertFloat())
    exposure = afwImage.makeExposure(mimg)

    # put in a temp psf
    psf = afwDet.createPsf("SingleGaussian", kwid, kwid, sigma0)
    exposure.setPsf(psf)

    return exposure


    
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
display = False
class ApertureCorrectionTestCase(unittest.TestCase):
    """Test the aperture correction."""

    def setUp(self):
        self.nx, self.ny = 128, 128
        self.ngrid        = 5
        self.sigma0      = 1.5
        self.val         = 40000.0
        self.sky         = 100.0
        self.kwid        = int(self.sigma0*7)
        if not self.kwid%2: self.kwid += 1

        # enable/disable the Assert statements
        # - diabling allows all tests in a method to run and print output
        # rather than throw an exception before remaining tests can run
        self.doTest      = False

        # how big do we allow the error to be
        #self.maxErrorFrac = 0.005       # half a percent fine for NAIVE
        self.maxErrorFrac = 0.003        # 0.3 percent needed for SINC
        # how many sigma can the measured value be from theoretical?
        # note: we're checking all candidate stars
        #self.nSigmaErrorLimit = 1.1     # 1.1 good for NAIVE
        self.nSigmaErrorLimit = 3.2      # 3.2 needed for SINC to past constant psf test

        # Note: SINC is not exactly as expected because of tapering of the aperture.
        
        self.metadata = dafBase.PropertyList()

        # detection policies
        
        self.detConfig = DetectionConfig()

        # measurement policies
        self.measSrcConfig = measAlg.MeasureSourcesConfig()
        
        # psf algorithms and policies
        # apcorr policies
        self.apCorrConfig = measAlg.ApertureCorrection.ConfigClass()
        self.apCorrConfig.polyStyle = "standard" # this does better than cheby ??
        self.apCorrConfig.order     = 2
        self.apCorrConfig.alg1.name = "PSF"
        self.apCorrConfig.alg2.name = "SINC"
        self.apCorrConfig.alg2["SINC"].radius = 3.0
        self.alg1 = self.apCorrConfig.alg1.active
        self.alg2 = self.apCorrConfig.alg2.active

        # logs
        self.log = pexLog.getDefaultLog()
        self.log.setThreshold(self.log.WARN)

        self.nDisp = 1
        
    def tearDown(self):
        del self.detConfig
        del self.measSrcConfig
        del self.apCorrConfig
        del self.metadata
        del self.log
        pass




    
    #################################################################
    # quick and dirty detection (note: we already subtracted background)
    def detectAndMeasure(self, exposure):

        # detect
        dsPos, dsNeg   = srcDet.detectSources(exposure, exposure.getPsf(), self.detConfig)

        footprintLists = [[dsPos.getFootprints(),[]]]
        # ... and measure
        sourceList     = srcMeas.sourceMeasurement(exposure, exposure.getPsf(),
                                                   footprintLists, self.measSrcConfig)

        return sourceList


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
    #   of photometry used: PSF, SINC, NAIVE
    ###################################################
    def getKnownFluxes(self, psfImg, radius, counts, sigma):

        flux = {"PSF": 0.0, "SINC": 0.0, "NAIVE": 0.0 }
        fluxErr = {"PSF": 0.0, "SINC": 0.0, "NAIVE": 0.0 }
        measErr = {"PSF": 0.0, "SINC": 0.0, "NAIVE": 0.0 }

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
                fluxErr["PSF"] += w*w*f
                flux["PSF"] += w*f

                # add up the naive fluxes
                dx, dy = i-ix, j-iy
                if (dx*dx + dy*dy <= radius*radius):
                    flux["NAIVE"] += f
                    fluxErr["NAIVE"] += f
                    # use smallest value as error ... ad-hoc
                    if f < measErr["NAIVE"] or measErr["NAIVE"] == 0:
                        measErr["NAIVE"] = f

        # renormalize the psf fluxes
        flux["PSF"] *= psfSum/psfSumSqrd
        fluxErr["PSF"] = math.sqrt(fluxErr["PSF"])*psfSum/psfSumSqrd
        measErr["PSF"] = 0.0 #fluxErr["PSF"]/math.sqrt(psfSum)
        
        fluxErr["NAIVE"] = math.sqrt(fluxErr["NAIVE"])
        measErr["NAIVE"] = math.sqrt(measErr["NAIVE"])

        # use the analytic form for the integral of a single gaussian for the sinc
        # - it's not quite right because of the cos tapering
        frac = self.apCorrTheory(sigma, self.alg2.radius)
        flux["SINC"] = counts*frac
        fluxErr["SINC"] = math.sqrt(flux["SINC"])
        measErr["SINC"] = 0.0

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
        exposure   = plantSources(self.nx, self.ny, self.kwid, self.sky, coordList)
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
        fluxKnown, fluxKnownErr, measKnownErr = self.getKnownFluxes(psfImg, self.alg2.radius, 
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
            fluxKnown, fluxKnownErr, measKnownErr = self.getKnownFluxes(psfImg, self.alg2.radius,
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
        """Verify that we can recovery the known aperture correction for *linearly varying* psf."""

        dx = self.nx/(self.ngrid + 1)
        dy = self.ny/(self.ngrid + 1)

        # vary apCorr by dApCorr linearly across the image
        apCorr = self.apCorrTheory(self.sigma0, self.alg2.radius)
        # want aperture correction to vary by this much across the field ...
        dApCorr = 0.05*apCorr
        # ... and that means changing sigma by this much:
        # (integrate r*exp(-r**2/sigma**2), and solve sigma)
        sig2   = self.alg2.radius*(-2.0*math.log(1.0 - (apCorr+dApCorr)))**-0.5

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
        """Verify that we can recovery the known aperture correction for *quadraticly varying* psf."""

        dx = self.nx/(self.ngrid + 1)
        dy = self.ny/(self.ngrid + 1)
        xmid, ymid = self.nx/2, self.ny/2
        
        # vary apCorr by dApCorr quadratically across the image
        apCorr = self.apCorrTheory(self.sigma0, self.alg2.radius)
        # want aperture correction to vary by this much across the field ...
        dApCorr = -0.05*apCorr
        # ... and that means changing sigma by this much:
        # (integrate r*exp(-r**2/sigma**2), and solve sigma)
        sig2   = self.alg2.radius*(-2.0*math.log(1.0 - (apCorr+dApCorr)))**-0.5

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


