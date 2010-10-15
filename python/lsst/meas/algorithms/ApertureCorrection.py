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

import glob, math, os, sys, re
from math import *
import numpy
import eups
import lsst.daf.base                   as dafBase
import lsst.pex.logging                as pexLog
import lsst.pex.policy                 as pexPolicy
import lsst.afw.detection              as afwDet
import lsst.afw.image                  as afwImage
import lsst.afw.math                   as afwMath
import lsst.meas.algorithms            as measAlg
import lsst.meas.algorithms.defects    as defects
import lsst.meas.algorithms.utils      as maUtils
import lsst.sdqa                       as sdqa

import lsst.afw.display.ds9            as ds9

import numpy.linalg                    as linalg
import lsst.meas.algorithms.Psf        as Psf


# to do:
# - allow instantiation with a psf, correction then based on direct measure of psf image.
# - allow spatialCells
# - add warnings about singular fit matrix



###################################################################
#
# Handle polynomial fits in 2d
#
###################################################################
class PolyFit2D(object):

    ##############################
    # constructor
    ##############################
    def __init__(self, x, y, z, order):
        
        self.order = order

        # get the exponents for this order
        self.exp = []
        for i in range(0, order+1):
            for j in range(0, i+1):
                xExp, yExp = j, i-j
                self.exp.append([xExp, yExp])

        # build a numpy array of the terms
        n = len(x)
        terms = []
        for i in range(len(self.exp)):
            xExp, yExp = self.exp[i]
            terms.append(x**xExp * y**yExp)

        # compute the least-squares fit
        # - note: the .T attribute is the transpose
        self.coeff, self.resid, self.rank, self.singval = linalg.lstsq(numpy.array(terms).T, z)


    ###################################
    # accessor to get the polyfit values at requested points
    ##############################
    def get(self, x, y):

        result = 0.0
        for i in range(len(self.exp)):
            xExp, yExp = self.exp[i]
            term = self.coeff[i] * x**xExp * y**yExp
            result += term
        return result



    
    

######################################################
#
# Class to manage aperture corrections
#
######################################################
class ApertureCorrection(object):

    #################
    # Constructor
    #################
    def __init__(self, exposure, sourceList, apCorrPolicy, psfPolicy, sdqaRatings, log=None, useAll=False):

        self.exposure     = exposure
        self.sourceList   = sourceList
        self.apCorrPolicy = apCorrPolicy
        self.psfPolicy    = psfPolicy
        self.sdqaRatings  = sdqaRatings
        self.log          = log

        # use a default log if we didn't get one
        if self.log is None:
            self.log = pexLog.getDefaultLog()
            self.log.setThreshold(pexLog.Log.WARN)
        self.log = pexLog.Log(self.log, "ApertureCorrection")

        # if we're not using all sources, select some suitable ones
        if not useAll:
            # for now, use the PSF star selection routine
            self.sourceList, self.cellSet = Psf.selectPsfSources(self.exposure,
                                                                 self.sourceList,
                                                                 self.psfPolicy)
        else:
            sizePsfCellX = psfPolicy.getInt("sizeCellX")
            sizePsfCellY = psfPolicy.getInt("sizeCellY")
            p0 = afwImage.PointI(exposure.getX0(), exposure.getY0())
            self.cellSet = afwMath.SpatialCellSet(afwImage.BBox(p0,
                                                                exposure.getWidth(),
                                                                exposure.getHeight()),
                                                  sizePsfCellX, sizePsfCellY)
            for s in self.sourceList:
                cand = measAlg.makePsfCandidate(s, exposure.getMaskedImage())
                self.cellSet.insertCandidate(cand)
                

        # unpack the policy    
        alg = [apCorrPolicy.get("algorithm1"), apCorrPolicy.get("algorithm2")]
        rad = [apCorrPolicy.get("radius1"),    apCorrPolicy.get("radius2")]
        self.order = apCorrPolicy.get("order")
        
        ###########
        # get the photometry for the requested algorithms
        self.mp = measAlg.makeMeasurePhotometry(self.exposure)
        for i in range(len(alg)):
            self.mp.addAlgorithm(alg[i])

            if rad[i] > 0.0:
                param = "%s.radius: %.1f" % (alg[i], rad[i])
            else:
                param = "%s.enabled: true" % (alg[i])
            policyStr = "#<?cfg paf policy?>\n%s\n" % (param)
            pol = pexPolicy.Policy.createPolicy(pexPolicy.PolicyString(policyStr))
            self.mp.configure(pol)


        ###########
        # test suitability of sources
        # - each source is ok
        # - distribution is complete

            
        ###########
        # get the point-to-point aperture corrections
        xList = numpy.array([])
        yList = numpy.array([])
        fluxList = [[],[]]
        self.apCorrList = numpy.array([])
        for s in self.sourceList:
            x, y = s.getXAstrom(), s.getYAstrom()
            
            try:
                p = self.mp.measure(afwDet.Peak(x, y))
            except Exception, e:
                self.log.log(log.WARN, "Failed to measure source at %.2f, %.2f." % (x, y))
                continue
                
            fluxes = []
            fluxErrs = []
            for a in alg:
                n = p.find(a)
                flux  = n.getFlux()
                fluxErr = n.getFluxErr()
                fluxes.append(flux)
                fluxErrs.append(fluxErr)

            apCorr = fluxes[1]/fluxes[0]
            self.log.log(self.log.INFO, "Using source: %7.2f %7.2f  %9.2f+/-%5.2f / %9.2f+/-%5.2f = %5.3f" %
                         (x, y, fluxes[0], fluxErrs[0], fluxes[1], fluxErrs[1], apCorr))
            
            fluxList[0].append(fluxes[0])
            fluxList[1].append(fluxes[1])

            xList = numpy.append(xList, x)
            yList = numpy.append(yList, y)
            self.apCorrList = numpy.append(self.apCorrList, apCorr)

        mean = numpy.mean(numpy.array(fluxList), axis=1)
        stdev = numpy.std(numpy.array(fluxList), axis=1)
        print "ap1: %.2f +/- %.2f" % (mean[0], stdev[0])
        print "ap2: %.2f +/- %.2f" % (mean[1], stdev[1])
        print "apCorr: %.3f +/- %.3f" % (numpy.mean(self.apCorrList), numpy.std(self.apCorrList))
        
            
        ###########
        # fit a polynomial to the aperture corrections
        self.fit = PolyFit2D(xList, yList, self.apCorrList, self.order)

        # if len(resid) == 0, the solution has too high an order for the number of sources
        if len(self.fit.resid) < 1:
            self.log.log(self.log.WARN, "Note enough stars for requested polyn. order in Aperture Correction.")


        ###########
        # Generate some stuff for SDQA
        numGoodStars  = 0
        numAvailStars = 0

        for cell in self.cellSet.getCellList():
            numGoodStars += cell.size()

        for cell in self.cellSet.getCellList():
            for cand in cell.begin(False):  # don't ignore BAD stars
                numAvailStars += 1
            
        sdqaRatings.append(sdqa.SdqaRating("phot.apCorr.numGoodStars", numGoodStars,
            0, sdqa.SdqaRating.AMP))
        sdqaRatings.append(sdqa.SdqaRating("phot.apCorr.numAvailStars",
            numAvailStars,  0, sdqa.SdqaRating.AMP))
        sdqaRatings.append(sdqa.SdqaRating("phot.apCorr.spatialLowOrdFlag", 0,  0,
            sdqa.SdqaRating.AMP))

        
    ###########################################
    # Accessor to get the apCorr at this x,y
    ###########################################
    def computeCorrectionAt(self, x, y):
        return self.fit.get(x, y)
