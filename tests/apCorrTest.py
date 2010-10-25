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
# - clean up plantFindSources - no need to add 1 to whole image and stack, use shift
# - test obvious orders: 0, 1, 2
# - try other more stable polyinterps ... cheby?
# - growth curves

import re
import os
import glob
import math
import pdb                          # we may want to say pdb.set_trace()
import unittest

import eups
import lsst.pex.exceptions      as pexEx
import lsst.pex.policy          as policy
import lsst.pex.logging         as pexLog
import lsst.afw.image           as afwImage
import lsst.afw.detection       as afwDet
import lsst.afw.geom            as afwGeom
import lsst.meas.algorithms     as algorithms
import lsst.utils.tests         as utilsTests
import lsst.sdqa                as sdqa

import numpy
import lsst.afw.math            as afwMath
import lsst.meas.algorithms.ApertureCorrection as apCorr
import lsst.meas.algorithms.Psf as Psf

import testLib

import lsst.afw.display.ds9       as ds9

try:
    type(verbose)
except NameError:
    verbose = 0

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
display = True
class ApertureCorrectionTestCase(unittest.TestCase):
    """Test the aperture correction."""

    def setUp(self):
        self.nx, self.ny = 128, 128
        self.ngrid        = 5
        self.sigma0      = 1.5
        self.val         = 20000.0
        self.sky         = 100.0
        self.alg1        = "PSF"
        self.alg2        = "NAIVE"
        self.rad1        = 0.0
        self.rad2        = 3.0
        self.kwid        = int(self.sigma0*7)
        if not self.kwid%2: self.kwid += 1

        # sdqa
        self.sdqaRatings = sdqa.SdqaRatingSet() # do I really need to make my own?

        # psf policies
        self.psfPolicy = policy.Policy.createPolicy(policy.DefaultPolicyFile("meas_algorithms", 
                                                                        "PsfDeterminationDictionary.paf",
                                                                        "policy"))
        self.psfAlgPolicy    = self.psfPolicy.get("psfPolicy")
        self.psfSelectPolicy = self.psfPolicy.get("selectionPolicy")
        self.psfSelectPolicy.set("sizeCellX", self.nx/4)
        self.psfSelectPolicy.set("sizeCellY", self.ny/4)


        # apcorr policies
        self.apCorrPolicy = policy.Policy.createPolicy(policy.DefaultPolicyFile("meas_algorithms", 
                                                                           "ApertureCorrectionDictionary.paf",
                                                                           "policy"))
        self.selectPolicy = self.apCorrPolicy.get("selectionPolicy")
        self.apCorrPolicy.set("order", 2)
        self.apCorrPolicy.set("algorithm1", self.alg1)
        self.apCorrPolicy.set("algorithm2", self.alg2)
        self.apCorrPolicy.set("radius1", self.rad1)
        self.apCorrPolicy.set("radius2", self.rad2)


        # logs
        self.log = pexLog.getDefaultLog()
        self.log.setThreshold(self.log.INFO)
        
        
    def tearDown(self):
        del self.psfAlgPolicy
        del self.psfSelectPolicy
        del self.psfPolicy
        del self.apCorrPolicy
        del self.selectPolicy
        del self.log
        pass



    def plantFindSources(self, coordList):

        # make an image and add fake stars
        img   = afwImage.ImageD(self.nx, self.ny, 0.0)
        msk   = afwImage.MaskU(img.getDimensions(), 0x0)
        msk.addMaskPlane("DETECTED")
        var   = afwImage.ImageD(self.nx, self.ny)

        # put delta functions in the image
        sigma0 = 0.0
        for coord in coordList:
            x, y, sigma = coord
            sigma0 += sigma

            # add a delta function
            imgDF = afwImage.ImageD(self.nx, self.ny, 0.0)
            imgDF.set(x, y, self.sky+self.val)

            # make a kernel
            gauss = afwMath.GaussianFunction2D(sigma, sigma)
            kernel = afwMath.AnalyticKernel(self.kwid, self.kwid, gauss)

            # convolve and add the final image
            imgPsf = afwImage.ImageD(self.nx, self.ny, 0.0)
            afwMath.convolve(imgPsf, imgDF, kernel)
            img += imgPsf
            
        img += self.sky
        sigma0 /= len(coordList)

        # add Poisson noise and mask the edge
        edgeBit = msk.getPlaneBitMask("EDGE")
        if True:
            ran = afwMath.Random()
            for j in range(self.ny):
                for i in range(self.nx):
                    img.set(i, j, ran.poisson(img.get(i, j)))

                    if (i < self.kwid or
                        i > self.nx - self.kwid or
                        j < self.kwid or
                        j > self.ny - self.kwid):
                        msk.set(i, j, edgeBit)
                    
        # make a maskedimage and an exposure
        var <<= img
        img -= self.sky
        mimg   = afwImage.MaskedImageF(img.convertFloat(),
                                       msk,
                                       var.convertFloat())
        exposure = afwImage.makeExposure(mimg)
        
        # put in a temp psf
        psf = afwDet.createPsf("SingleGaussian", self.kwid, self.kwid, sigma0) #FWHM/(2*sqrt(2*log(2))))
        exposure.setPsf(psf)

        
        ####
        # quick and dirty detection
        cnvImage = mimg.Factory(mimg.getDimensions())
        afwMath.convolve(cnvImage, mimg, kernel, afwMath.ConvolutionControl())
        llc = afwImage.PointI(kernel.getWidth()/2, kernel.getHeight()/2)
        urc = afwImage.PointI(cnvImage.getWidth() - 1, cnvImage.getHeight() - 1) - llc;
        middle = cnvImage.Factory(cnvImage, afwImage.BBox(llc, urc))

        threshold = afwDet.Threshold(3, afwDet.Threshold.STDEV)
        ds = afwDet.FootprintSetF(middle, threshold, "DETECTED")
        ds.setMask(mimg.getMask(), "DETECTED")
        del middle
        objects = ds.getFootprints()

        ####
        # quick and dirty measurement
        moPolicy = policy.Policy.createPolicy(os.path.join(eups.productDir("meas_algorithms"),
                                                           "examples", "MeasureSources.paf"))
        moPolicy = moPolicy.getPolicy("measureObjects")
        measureSources = algorithms.makeMeasureSources(exposure, moPolicy)

        sourceList = afwDet.SourceSet()
        for i in range(len(objects)):
            source = afwDet.Source()
            sourceList.append(source)

            source.setId(i)
            source.setFlagForDetection(source.getFlagForDetection() | algorithms.Flags.BINNED1);
            measureSources.apply(source, objects[i])
        
        
        return exposure, sourceList, kernel



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
        frac = 1.0 - math.exp(-radius**2/(2.0*sigma**2))
        flux["SINC"] = counts*frac
        fluxErr["SINC"] = math.sqrt(flux["SINC"])
        measErr["SINC"] = 0.0

        return flux, fluxErr, measErr

    
    def getKnownApCorr(self, fluxKnown, fluxKnownErr, measKnownErr):
        apCorr    = fluxKnown[self.alg2]/fluxKnown[self.alg1]
        apCorrErr = apCorr*(measKnownErr[self.alg1]/fluxKnown[self.alg1] +
                            measKnownErr[self.alg2]/fluxKnown[self.alg2])
        return apCorr, apCorrErr

    
    def printSummary(self, psfImg, fluxKnown, fluxKnownErr, measKnownErr):
    
        # print diagnostics on the star selection
        sdqaRatings = dict(zip([r.getName() for r in self.sdqaRatings], [r for r in self.sdqaRatings]))
        print "Used %d apCorr stars (%d good)" % (sdqaRatings["phot.apCorr.numAvailStars"].getValue(),
                                                  sdqaRatings["phot.apCorr.numGoodStars"].getValue())
        
        # have a look at the know values
        print "Flux known (%s): %.2f +/- %.2f" % (self.alg1, fluxKnown[self.alg1], fluxKnownErr[self.alg1])
        print "Flux known (%s): %.2f +/- %.2f" % (self.alg2, fluxKnown[self.alg2], fluxKnownErr[self.alg2])
        apCorr, apCorrErr    = self.getKnownApCorr(fluxKnown, fluxKnownErr, measKnownErr)
        print "Aperture Corr'n Known: %.3f +/- %.3f" % (apCorr, apCorrErr)

        
    def testApCorrConstantPsf(self):
        """Test that we can model the corrections for fake objects"""

        dx = self.nx/(self.ngrid + 1)
        dy = self.ny/(self.ngrid + 1)

        # decide where to put fake psfs on a grid
        coordList = []
        for i in range(self.ngrid):
            for j in range(self.ngrid):
                x, y = (1+i)*dx, (1+j)*dy
                coordList.append([x, y, self.sigma0])

        # plant them in the image, and measure them
        exposure, sourceList, kernel = self.plantFindSources(coordList)
        mimg = exposure.getMaskedImage()
        img = mimg.getImage()


        # try getPsf()
        psf, cellSet, psfSourceSet = Psf.getPsf(exposure, sourceList, self.psfPolicy, self.sdqaRatings)
        exposure.setPsf(psf)

        
        # try apCorr()
        ac = apCorr.ApertureCorrection(exposure, psfSourceSet,
                                       self.apCorrPolicy, self.selectPolicy, self.sdqaRatings,
                                       self.log, useAll=True)

        normPeak = False
        psfImg = psf.computeImage(afwGeom.makePointD(mimg.getWidth()/2, mimg.getHeight()/2), normPeak)
        fluxKnown, fluxKnownErr, measKnownErr = self.getKnownFluxes(psfImg, self.rad2, self.val, self.sigma0)
        self.printSummary(psfImg, fluxKnown, fluxKnownErr, measKnownErr)

        corrKnown, corrErrKnown           = self.getKnownApCorr(fluxKnown, fluxKnownErr, measKnownErr)
        corrMeasMiddle, corrErrMeasMiddle = ac.computeCorrectionAt(self.nx/2, self.ny/2)

        
        ###################
        # Tests
        ###################
        self.assertAlmostEqual(corrKnown, corrMeasMiddle, 4)
        #self.assertAlmostEqual(corrErrKnown, corrErrMeasMiddle, 4)
        
        if display:

            # show the apCorr and error as images
            acImg = afwImage.ImageF(self.nx, self.ny)
            acErrImg = afwImage.ImageF(self.nx, self.ny)
            for j in range(self.ny):
                for i in range(self.nx):
                    apCo, apCoErr = ac.computeCorrectionAt(i, j)
                    acImg.set(i, j, apCo)
                    acErrImg.set(i, j, apCoErr)

            ds9.mtv(img,      frame=1, title="Delta functions")
            ds9.mtv(mimg,     frame=2, title="convolved image")
            ds9.mtv(acImg,    frame=3, title="Apcorr Image")
            ds9.mtv(acErrImg, frame=4, title="Apcorr Error Image")
            ds9.mtv(psfImg,   frame=5, title="Psf Image")

        

            
        
        
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


