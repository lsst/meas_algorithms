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

import glob
import math
import os
import sys
import re

import numpy
import numpy.linalg as linalg

import eups
import lsst.daf.base as dafBase
import lsst.pex.logging as pexLog
import lsst.pex.policy as pexPolicy
import lsst.pex.config as pexConfig
import lsst.afw.detection as afwDet
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom

from . import algorithmsLib as measAlg
from . import measurement

import utils as maUtils

import lsst.afw.display.ds9 as ds9

# to do:
# - allow instantiation with a psf, correction then based on direct measure of psf image.
# - figure out why 'standard' polynomials are better than 'cheby'

##################################################################
#
# Generate the first 'order' terms in a requested polynomial style
#
##################################################################
class Poly1D(object):
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
    def __init__(self, x, y, z, w, poly):

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
        self.coeff, self.resid, self.rank, self.singval = self._fit(z, w, xTerms, yTerms, self.orderPairs)
        
        ################
        # errors
        # - done separately as errOrder may differ
        # - put the squared residuals in a vector to fit
        self.residuals = numpy.array([])
        for i in range(len(z)):
            dz = z[i] - self.getVal(x[i], y[i])
            self.residuals = numpy.append(self.residuals, dz*dz)
            
        self.errOrderPairs = self._computeOrderPairs(self.errOrder)
        self.errCoeff, self.errResid, self.errRank, self.errSingval = self._fit(self.residuals, w*w,
                                                                                xTerms, yTerms,
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
    def _fit(self, z, w, xTerms, yTerms, orderPairs):
        
        # take the products of the x,y terms to build a numpy array for the linear fit
        terms = []
        for i in range(len(orderPairs)):
            xOrd, yOrd = orderPairs[i]
            terms.append(w * xTerms[xOrd] * yTerms[yOrd])
        # - note: the .T attribute is the transpose
        terms = numpy.array(terms).T
        coeff, resid, rank, singval = linalg.lstsq(terms, w*z)
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


class ApertureCorrectionConfig(pexConfig.Config):
    """Config for ApertureCorrection
    """
    polyStyle = pexConfig.Field(
        doc = "Style of polynomial fit",
        dtype = str,
        default = "standard",
    )
    order = pexConfig.Field(
        doc = "Polynomial interpolation order across the chip.",
        dtype = int,
        default = 2,
    )
    alg1 = measurement.AlgorithmRegistry.filter(measurement.FluxConfig).makeField(
        doc = "Photometric algorithm 1 (aperture correct _from_ this algorithm).",
        multi = False,
        default = "flux.psf",
    )
    alg2 = measurement.AlgorithmRegistry.filter(measurement.FluxConfig).makeField(
        doc = "Photometric algorithm 2 (aperture correct _to_ this algorithm).",
        default = "flux.sinc",
    )
        
######################################################
#
# Class to manage aperture corrections
#
######################################################
class ApertureCorrection(object):
    """Class to manage aperture corrections.

    @param exposure    an afw.Exposure containing the sources to use to compute the aperture correction
    @param cellSet     an afw.math.spatialCellSet containing coords to use
    @param metadata    somewhere to put interesting information 
    @param config      configuration; an instance of self.ConfigClass
    @param log         a pex.logging log
        
    If a spatialCellSet is provided, it is assumed that no further selection is required,
    as a cellSet does not contain sufficient information to select candidates (namely fluxes).
    
    If a sourceSet is provided, selectPolicy must be provided as it contains sizeCellX/Y
    needed to create a spatialCellSet from a sourceSet.
        
    """
    ConfigClass = ApertureCorrectionConfig

    #################
    # Constructor
    #################
    def __init__(self, exposure, cellSet, metadata, config, log=None):

        import lsstDebug
        display = lsstDebug.Info(__name__).display

        if display:
            frame = 0
            ds9.mtv(exposure, frame=frame, title="Exposure for calibration")

        self.xwid, self.ywid = exposure.getWidth(), exposure.getHeight()

        # use a default log if we didn't get one
        if log is None:
            log = pexLog.getDefaultLog()
            log.setThreshold(pexLog.Log.WARN)
        log = pexLog.Log(log, "ApertureCorrection")

        # unpack the control object
        self.order     = config.order
        self.polyStyle = config.polyStyle

        ###########
        # get the photometry for the requested algorithms
        mp = measAlg.MeasureSources()
        mp.addAlgorithm(config.alg1.apply())
        mp.addAlgorithm(config.alg2.apply())
        table = afwTable.SourceTable.make(mp.getSchema())

        ###########
        # get the point-to-point aperture corrections
        xList = numpy.array([])
        yList = numpy.array([])
        fluxList = [[],[]]
        self.apCorrList = numpy.array([])
        self.apCorrErrList = numpy.array([])
        for cell in cellSet.getCellList():
            for cand in cell.begin(True): # ignore bad candidates
                cand = measAlg.cast_PsfCandidateF(cand)
                s = cand.getSource()

                if s.get("flags.pixel.interpolated.center"):
                    continue

                x, y = cand.getXCenter(), cand.getYCenter()
                
                source = table.makeRecord()
                source.setFootprint(s.getFootprint())
                center = afwGeom.Point2D(x, y)

                mp.apply(source, exposure, center)

                fluxes = [source[config.alg1.name], source[config.alg2.name]]
                fluxErrs = [source[config.alg1.name + ".err"], source[config.alg2.name + ".err"]]

                if fluxes[0] <= 0.0 or fluxes[1] <= 0.0:
                    log.log(log.WARN, "Non-positive flux for source at %.2f,%.2f (%f,%f)" %
                            (x, y, fluxes[0], fluxes[1]))
                    continue

                apCorr = fluxes[1]/fluxes[0]
                apCorrErr = apCorr*math.sqrt( (fluxErrs[0]/fluxes[0])**2 + (fluxErrs[1]/fluxes[1])**2 )
                log.log(log.DEBUG,
                             "Using source: %7.2f %7.2f  %9.2f+/-%5.2f / %9.2f+/-%5.2f = %5.3f+/-%5.3f" %
                             (x, y, fluxes[0], fluxErrs[0], fluxes[1], fluxErrs[1], apCorr, apCorrErr))
                if numpy.isnan(apCorr) or numpy.isnan(apCorrErr):
                    continue

                if display:
                    size = rad[0]
                    if size == 0:
                        size = rad[1]
                    ds9.dot("o", x, y, size=size, ctype=ds9.WHITE, frame=frame)
                    ds9.dot("%.3f" % (apCorr), x, y - size - 10, ctype=ds9.WHITE, frame=frame)
                    ds9.dot("%d" % (cand.getId()), x, y + size + 10, ctype=ds9.WHITE, frame=frame)

                fluxList[0].append(fluxes[0])
                fluxList[1].append(fluxes[1])

                xList = numpy.append(xList, x)
                yList = numpy.append(yList, y)
                self.apCorrList = numpy.append(self.apCorrList, apCorr)
                self.apCorrErrList = numpy.append(self.apCorrErrList, apCorrErr)

        if len(self.apCorrList) == 0:
            raise RuntimeError("No good aperture correction measurements.")                

        ###########
        # fit a polynomial to the aperture corrections
        self.fitOrder = self.order
        # if cheby, we'll overfit and truncate
        if re.search("cheby", self.polyStyle, re.IGNORECASE):
            self.fitOrder += 1
            
        poly     = Poly1D(self.fitOrder, style=self.polyStyle)
        weights  = (1.0/self.apCorrErrList)**2
        self.fit = PolyFit2D(xList, yList, self.apCorrList, weights, poly)


        # do a 3-sigma clip and refit
        if False:
            # This currently only works with fake data and is disabled
            # with real data, it clips out almost all measurements
            apCorrFit    = self.fit.getVal(xList, yList)
            apCorrErrFit = self.fit.getErr(xList, yList)
            nSig = 3.0
            igood = numpy.where( numpy.abs(self.apCorrList - apCorrFit) < nSig*apCorrErrFit )[0]
            self.fit = PolyFit2D(xList[igood], yList[igood], self.apCorrList[igood], w[igood], poly)
        
            
        ###########
        # check sanity
        
        # if len(resid) == 0, the solution has too high an order for the number of sources
        if len(self.fit.resid) < 1:
            log.log(log.WARN,
                         "Not enough stars for requested polyn. order in Aperture Correction.")
        
        # warn if singular
        # press et al says (in SVD section): thresh = 0.5*sqrt(m+n+1.)*w[0]*epsilon
        mDim     = len(self.fit.coeff)
        nDim     = len(xList)
        epsilon  = 1.0e-15 #... machine precision
        svThresh = 0.5*math.sqrt(mDim+nDim+1.0)*self.fit.singval[0]*epsilon

        if self.fit.singval[-1] < svThresh:
            log.log(log.WARN, "Singular value below threshold in apCorr fit (%.14f < %.14f)" %
                         (sv, svThresh))

            
        ###########
        # Generate some stuff for QA
        numGoodStars  = len(self.apCorrList)
        numAvailStars = 0

        for cell in cellSet.getCellList():
            for cand in cell.begin(True):  # ignore BAD stars ... they're not available really
                numAvailStars += 1
            
        metadata.set("numGoodStars", numGoodStars)
        metadata.set("numAvailStars", numAvailStars)

        log.log(log.INFO, "%s to %s" % (config.alg1.name, config.alg2.name))
        log.log(log.INFO, "numGoodStars: %d" % (numGoodStars))
        log.log(log.INFO, "numAvailStars: %d" % (numAvailStars))
        #mean = numpy.mean(numpy.array(fluxList), axis=1)
        #stdev = numpy.std(numpy.array(fluxList), axis=1)
        #log.log(log.INFO, "mean ap1: %.2f +/- %.2f" % (mean[0], stdev[0]))
        #log.log(log.INFO, "mean ap2: %.2f +/- %.2f" % (mean[1], stdev[1]))
        log.log(log.INFO, "mean apCorr: %.4f +/- %.4f" %
                     (numpy.mean(self.apCorrList), numpy.std(self.apCorrList)))
        x, y = self.xwid/2, self.ywid/2
        log.log(log.INFO, "apCorr(%d,%d): %.4f +/- %.4f" %
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
