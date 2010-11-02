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

# to do:
# - allow instantiation with a psf, correction then based on direct measure of psf image.
# - figure out why 'standard' polynomials are better than 'cheby'

##################################################################
#
# Generate the first 'order' terms in a requested polynomial style
#
##################################################################
class PolyGenerator(object):
    """
    Generate the first n terms for a requested polynomial style (standard or cheby)
    x     = numpy.ndarray(myPythonArray)
    polyS = PolyGenerator(3, 'standard')
    terms = polyS.getTerms(x)
    """
    
    def __init__(self, order, style="standard"):
        self.order = order
        self.style = style
        
    def getTerms(self, x):
        """Return a numpy array containing the first n terms in the polynomial expansion of x"""
        
        if isinstance(x, numpy.ndarray):
            terms = [numpy.ones(len(x)), x]
        else:
            terms = [1, x]

        if re.search("cheby", self.style, re.IGNORECASE):
            for i in range(2, self.order+1):
                terms.append(2.0*x*terms[i-1] - terms[i-2])
        else:
            for i in range(2, self.order+1):
                terms.append(x*terms[i-1])

        return terms

    
    
###################################################################
#
# Handle polynomial fits in 2d
#
###################################################################
class PolyFit2D(object):
    """
    Handle polynomial fitting in 2d
    
    strategy is to fit:
    a0 *poly0(x)*poly0(y) +
    a1x*poly1(x)*poly0(y) + a1y*poly0(x)*poly1(y) + 
    a2xx*poly2(x)*poly0(y) + a2xy*poly1(x)*poly1(y) + a2yy*poly0(x)*poly2(y) + ...

    where polyN() is the nth-order polynomial of type 'poly'
    eg. for standard polynomical poly2(x) = x**2, but
        for chebyshev            poly2(x) = 2*x**2 - 1
    """

    
    ##############################
    # constructor
    ##############################
    def __init__(self, x, y, z, poly):

        self.poly  = poly
        self.order = self.poly.order
        self.errOrder = 0

        # compute the polynomical terms for x and y
        xTerms, yTerms = self.poly.getTerms(x), self.poly.getTerms(y)
                
        ################
        # values
        # get the order-pairs for all terms (including) cross terms
        self.orderPairs = self._computeOrderPairs(self.order)
        # compute the least-squares fit
        self.coeff, self.resid, self.rank, self.singval = self._fit(z, xTerms, yTerms, self.orderPairs)
        
        ################
        # errors
        # - done separately as errOrder may differ
        # - put the squared residuals in a vector to fit
        residuals = numpy.array([])
        for i in range(len(z)):
            dz = z[i] - self.getVal(x[i], y[i])
            residuals = numpy.append(residuals, dz*dz)
            
        self.errOrderPairs = self._computeOrderPairs(self.errOrder)
        self.errCoeff, self.errResid, self.errRank, self.errSingval = self._fit(residuals, xTerms, yTerms,
                                                                                self.errOrderPairs)


    ##################################
    # get the order-pairs for all terms (including) cross terms
    # these are the x,y polynomial orders for each term, including cross terms:
    # ie. 0th order: [0, 0],
    #     1st order: [1, 0], [0, 1]
    #     2nd order: [2, 0], [1, 1], [0, 2]
    #     etc.
    ##################################
    def _computeOrderPairs(self, order):
    
        # get the order-pairs for all terms (including) cross terms
        orderPairs = []
        for i in range(0, order+1):
            for j in range(0, i+1):
                xOrd, yOrd = j, i-j
                orderPairs.append([xOrd, yOrd])
        return orderPairs

    
    ##################################
    # compute the polynomical terms for x and y
    # eg. for 2nd order x, get: [1, x, x**2]
    # pass everything to lstsq() and return what we get
    ##################################
    def _fit(self, z, xTerms, yTerms, orderPairs):
        
        # take the products of the x,y terms to build a numpy array for the linear fit
        terms = []
        for i in range(len(orderPairs)):
            xOrd, yOrd = orderPairs[i]
            terms.append(xTerms[xOrd] * yTerms[yOrd])
        # - note: the .T attribute is the transpose
        terms = numpy.array(terms).T
        coeff, resid, rank, singval = linalg.lstsq(terms, z)
        return coeff, resid, rank, singval
        
        
    ###################################
    # accessor to get the polyfit values at requested points
    ##############################
    def getVal(self, x, y, truncOrder=None):
        """Get the value at x, y"""
        
        if truncOrder is None:
            truncOrder = self.order
        if truncOrder > self.order:
            raise AttributeError, ("truncOrder must be <= original order of polynomical.")
        
        result = 0.0
        xTerms = self.poly.getTerms(x)
        yTerms = self.poly.getTerms(y)
        for i in range(len(self.orderPairs)):
            xOrd, yOrd = self.orderPairs[i]
            if xOrd+yOrd > truncOrder:
                continue
            term = self.coeff[i] * xTerms[xOrd]*yTerms[yOrd]
            result += term
        return result

    
    ###################################
    # accessor to get the polyfit values at requested points
    # This could be done as part of get, but we might want errOrder to be different from order
    # ... in which case, it must be separate.
    ##############################
    def getErr(self, x, y, truncOrder=None):
        """Get the error at x, y"""

        if truncOrder is None:
            truncOrder = self.order
        if truncOrder > self.order:
            raise AttributeError, ("truncOrder must be <= original order of polynomical.")
        
        err = 0.0
        xTerms = self.poly.getTerms(x)
        yTerms = self.poly.getTerms(y)
        for i in range(len(self.errOrderPairs)):
            xOrd, yOrd = self.errOrderPairs[i]
            if xOrd+yOrd > truncOrder:
                continue
            term = self.errCoeff[i] * xTerms[xOrd]*yTerms[yOrd]
            err += term
        return numpy.sqrt(err)


    
#####################################################
# Control Object for aperture correction
#####################################################
class ApertureCorrectionControl(object):
    """
    Handle input parameters for Aperture Control

    This is a thin replacement for a policy.  The constructor accepts only a policy.
    """
    
    # construct
    def __init__(self, policy):
        self.alg1      = policy.get("algorithm1")
        self.alg2      = policy.get("algorithm2")
        self.rad1      = policy.get("radius1")
        self.rad2      = policy.get("radius2")
        self.polyStyle = policy.get("polyStyle")
        self.order     = policy.get("order")

        
######################################################
#
# Class to manage aperture corrections
#
######################################################
class ApertureCorrection(object):
    """Class to manage aperture corrections.

    exposure    = an afw.Exposure containing the sources to use to compute the aperture correction
    sources     = either a afw.detection.sourceSet, or an afw.math.spatialCellSet
    sdqaRatings = self-explanatory
    apCorrCtrl  = An ApertureControl object (created with an aperture control policy)
    selectPolicy = a policy for selecting candidate stars (not required if a spatialCellSet is provided)
    log         = a pex.logging log
    doSelect    = True = select suitable stars, False = use sources provided.
        
    If a spatialCellSet is provided, it is assumed that no further selection is required,
    as a cellSet does not contain sufficient information to select candidates (namely fluxes).
    
    If a sourceSet is provided, selectPolicy must be provided as it contains sizeCellX/Y
    needed to create a spatialCellSet from a sourceSet.
        
    """
    
    #################
    # Constructor
    #################
    def __init__(self, exposure, sources, sdqaRatings, apCorrCtrl, selectPolicy=None,
                 log=None, doSelect=True):

        self.exposure     = exposure
        self.sources      = sources
        self.apCorrCtrl   = apCorrCtrl
        self.selectPolicy = selectPolicy
        self.sdqaRatings  = sdqaRatings
        self.log          = log

        self.xwid, self.ywid = self.exposure.getWidth(), self.exposure.getHeight()

        # use a default log if we didn't get one
        if self.log is None:
            self.log = pexLog.getDefaultLog()
            self.log.setThreshold(pexLog.Log.WARN)
        self.log = pexLog.Log(self.log, "ApertureCorrection")

        # unpack the control object
        alg = [apCorrCtrl.alg1, apCorrCtrl.alg2]
        rad = [apCorrCtrl.rad1, apCorrCtrl.rad2]
        self.order     = apCorrCtrl.order
        self.polyStyle = apCorrCtrl.polyStyle

        
        ###########
        # if we're not using all sources, select some suitable ones

        # if they've given us a cellSet, we're good to go
        if isinstance(sources, afwMath.SpatialCellSet):
            self.cellSet = sources

        # if they've given us a sourceSet, make sure we got a policy to let us
        # either select stars in a cellSet, or just put stars in a cellSet
        elif isinstance(sources, afwDet.SourceSet):

            if selectPolicy is None:
                raise RuntimeError, (
                    "Must provide a selectionPolicy when 'sources' input is an afw.sourceSet")

            # if they want us to pick the stars, call a selection routine.
            if doSelect:
                pkg = selectPolicy.get("package")
                __import__(pkg)
                selectPsf = sys.modules[pkg]
                self.sources, self.cellSet = selectPsf.selectPsfSources(self.exposure,
                                                                        self.sources,
                                                                        self.selectPolicy)

            # if they want us to just use what we were given, insert sources into a cellSet
            else:
                sizePsfCellX = self.selectPolicy.getInt("sizeCellX")
                sizePsfCellY = self.selectPolicy.getInt("sizeCellY")
                p0 = afwImage.PointI(exposure.getX0(), exposure.getY0())
                self.cellSet = afwMath.SpatialCellSet(afwImage.BBox(p0,
                                                                    exposure.getWidth(),
                                                                    exposure.getHeight()),
                                                      sizePsfCellX, sizePsfCellY)
                for s in self.sources:
                    cand = measAlg.makePsfCandidate(s, exposure.getMaskedImage())
                    self.cellSet.insertCandidate(cand)
                    
        # if they gave us neither a sourceSet, nor a cellSet ... bark
        else:
            raise AttributeError, ("'sources' must be an afwDetection.SourceSet or an afwMath.SpatialCellSet")

        
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
        # get the point-to-point aperture corrections
        xList = numpy.array([])
        yList = numpy.array([])
        fluxList = [[],[]]
        self.apCorrList = numpy.array([])
        for cell in self.cellSet.getCellList():
            for cand in cell.begin(True): # ignore bad candidates
                x, y = cand.getXCenter(), cand.getYCenter()
                
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
                self.log.log(self.log.INFO,
                             "Using source: %7.2f %7.2f  %9.2f+/-%5.2f / %9.2f+/-%5.2f = %5.3f" %
                             (x, y, fluxes[0], fluxErrs[0], fluxes[1], fluxErrs[1], apCorr))

                fluxList[0].append(fluxes[0])
                fluxList[1].append(fluxes[1])

                xList = numpy.append(xList, x)
                yList = numpy.append(yList, y)
                self.apCorrList = numpy.append(self.apCorrList, apCorr)


                
        ###########
        # fit a polynomial to the aperture corrections
        self.fitOrder = self.order
        # if cheby, we'll overfit and truncate
        if re.search("cheby", self.polyStyle, re.IGNORECASE):
            self.fitOrder += 1
            
        poly = PolyGenerator(self.fitOrder, style=self.polyStyle)
        self.fit = PolyFit2D(xList, yList, self.apCorrList, poly)


        ###########
        # check sanity
        
        # if len(resid) == 0, the solution has too high an order for the number of sources
        if len(self.fit.resid) < 1:
            self.log.log(self.log.WARN,
                         "Not enough stars for requested polyn. order in Aperture Correction.")
        
        # warn if singular
        # press et al says (in SVD section): thresh = 0.5*sqrt(m+n+1.)*w[0]*epsilon
        mDim     = len(self.fit.coeff)
        nDim     = len(xList)
        epsilon  = 1.0e-15 #... machine precision
        svThresh = 0.5*math.sqrt(mDim+nDim+1.0)*self.fit.singval[0]*epsilon

        if self.fit.singval[-1] < svThresh:
            self.log.log(self.log.WARN, "Singular value below threshold in apCorr fit (%.14f < %.14f)" %
                         (sv, svThresh))

            
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


        self.log.log(self.log.INFO, "numGoodStars: %d" % (numGoodStars))
        self.log.log(self.log.INFO, "numAvailStars: %d" % (numAvailStars))
        mean = numpy.mean(numpy.array(fluxList), axis=1)
        stdev = numpy.std(numpy.array(fluxList), axis=1)
        self.log.log(self.log.INFO, "mean ap1: %.2f +/- %.2f" % (mean[0], stdev[0]))
        self.log.log(self.log.INFO, "mean ap2: %.2f +/- %.2f" % (mean[1], stdev[1]))
        self.log.log(self.log.INFO, "mean apCorr: %.3f +/- %.3f" %
                     (numpy.mean(self.apCorrList), numpy.std(self.apCorrList)))
        x, y = self.xwid/2, self.ywid/2
        self.log.log(self.log.INFO, "apCorr(%d,%d): %.3f +/- %.3f" %
                     (x, y, self.fit.getVal(x,y,self.order), self.fit.getErr(x,y,self.order)))
        
        
        
    ###########################################
    # Accessor to get the apCorr at this x,y
    ###########################################
    def computeAt(self, x, y):
        """
        Compute the aperture correction and its error at x, y
        Return [value, error]
        """
        return self.fit.getVal(x, y, self.order), self.fit.getErr(x, y, self.order)
