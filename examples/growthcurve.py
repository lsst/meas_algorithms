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

#
# Original filename: examples/growthcurve.py
#
# Author: Steve Bickerton
# Email: bick@astro.princeton.edu
# Date: Mon 2009-10-26 13:42:37
# 
# Summary:
#
# python version growthcurve.cc example
# 
"""
%prog [options] arg
"""

import sys
import re
import optparse
import os
import datetime
import numpy
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.detection as afwDet
import lsst.meas.algorithms as algorithms

# =====================================================================
# a functor for the PSF
#
class Gaussian(): # public std::binary_function<double, double, double> {
    def __init__(self, xcen, ycen, sigma, a) :
        self.xcen = xcen
        self.ycen = ycen
        self.sigma = sigma
        self.a  = a
    def __call__(self, x, y):
        xx = x - self.xcen
        yy = y - self.ycen
        ss = self.sigma*self.sigma
        coeff = self.a * (1.0/(2.0*numpy.pi*ss))
        expon = numpy.exp(-(xx*xx + yy*yy) / (2.0*ss))
        return coeff*expon


# =====================================================================
# a radial functor for the PSF
#
# This functor isn't currently used in the routine
# I'll leave it here in case I (someday) figure out how to integrate a python functor
class RGaussian(): #public std::unary_function<double, double> {

    def __init__(self, sigma, a, apradius, aptaper):
        self.sigma = sigma
        self.a = a
        self.apradius = apradius
        self.aptaper = aptaper

    def __call__ (self, r):
        ss = self.sigma*self.sigma
        gauss = self.a * (1.0/(2.0*numpy.pi*ss)) * numpy.exp(-(r*r)/(2.0*ss));
        aperture = 0.0
        if ( r <= apradius ):
            aperture = 1.0
        elif ( r > apradius and r < apradius + aptaper ):
            aperture = 0.5*(1.0 + cos(numpy.pi*(r - apradius)/aptaper))
        return aperture*gauss*(r*2.0*numpy.pi)


#############################################################
#
# Main body of code
#
#############################################################

def main():

    date = datetime.datetime.now().strftime("%a %Y-%m-%d %H:%M:%S")

    ########################################################################
    # command line arguments and options
    ########################################################################
    
    parser = optparse.OptionParser(usage = __doc__)
    #parser.add_option("-a", "--aa", dest="aa", type=float,
    #                  default=1.0, help="default=%default")
    
    opts, args = parser.parse_args()

    if len(args) == 0:
        r1, r2, dr = 3.0, 3.0, 0.5
    elif len(args) == 3:
        r1, r2, dr = map(float, args)
    else:
        parser.print_help()
        sys.exit(1)

        
    # make a list of radii to compute the growthcurve points
    radius = []
    nR = int( (r2 - r1)/dr + 1 )
    for iR in range(nR):
        radius.append(r1 + iR*dr)


    # make an image big enough to hold the largest requested aperture
    xwidth = 2*(0 + 128)
    ywidth = xwidth

    # initializations
    sigmas = [1.5, 2.5]  # the Gaussian widths of the psfs we'll use
    nS = len(sigmas)
    a = 100.0
    aptaper = 2.0
    xcen = xwidth/2
    ycen = ywidth/2

    
    print "# sig rad  Naive Sinc Psf"
    for iS in range(nS):
        sigma = sigmas[iS];

        # make a Gaussian star to measure
        gauss  = afwMath.GaussianFunction2D(sigma, sigma)
        kernel = afwMath.AnalyticKernel(xwidth, ywidth, gauss)
        kimg   = afwImage.ImageD(kernel.getDimensions())
        kernel.computeImage(kimg, False)
        kimg  *= 100.0
        mimg   = afwImage.MaskedImageF(kimg.convertFloat(),
                                       afwImage.MaskU(kimg.getDimensions(), 0x0),
                                       afwImage.ImageF(kimg.getDimensions(), 0.0))

        # loop over all the radii in the growthcurve
        for iR in range(nR):

            psfH = int(2.0*(r2 + 2.0))
            psfW = int(2.0*(r2 + 2.0))

            psf = afwDet.DoubleGaussianPsf(psfW, psfH, sigma)

            # get the aperture fluxes for Naive and Sinc methods
            mpNaive   = algorithms.makeMeasurePhotometry("NAIVE", radius[iR])
            photNaive = mpNaive.apply(mimg, xcen, ycen, psf, 0.0)
            mpSinc    = algorithms.makeMeasurePhotometry("SINC", radius[iR])
            photSinc  = mpSinc.apply(mimg, xcen, ycen, psf, 0.0)

            fluxNaive = photNaive.getApFlux()
            fluxSinc  = photSinc.getApFlux()
            fluxPsf   = photSinc.getPsfFlux()
            
            # get the exact flux for the theoretical smooth PSF
            # rpsf = RGaussian(sigma, a, radius[iR], aptaper)
            # *** not sure how to integrate a python functor ***
            
            print "%.2f %.2f  %.3f %.3f %.3f" % (sigma, radius[iR], fluxNaive, fluxSinc, fluxPsf)
        

#############################################################
# end
#############################################################

if __name__ == '__main__':
    main()
