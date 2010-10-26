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

# This is not a minimal set of imports
import glob, math, os, sys
from math import *
import numpy
import eups
import lsst.daf.base as dafBase
import lsst.pex.logging as logging
import lsst.pex.policy as policy
import lsst.afw.detection as afwDetection
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.meas.algorithms as algorithms
import lsst.meas.algorithms.defects as defects
import lsst.meas.algorithms.utils as maUtils
import lsst.sdqa as sdqa

import lsst.afw.display.ds9 as ds9


    
    
def getPsf(exposure, sourceList, psfPolicy, sdqaRatings):
    """Return the PSF for the given Exposure and set of Sources, given a Policy

The policy is documented in ip/pipeline/policy/CrRejectDictionary.paf    
    """

    #############################################
    # select the stars
    psfSelectPolicy = psfPolicy.get("selectionPolicy")
    
    selectPackage   = psfSelectPolicy.get("package")
    __import__(selectPackage)
    psfSel    = sys.modules[selectPackage]
    psfStars, psfCellSet = psfSel.selectPsfSources(exposure, sourceList, psfSelectPolicy)
    

    #############################################
    # get the psf with the chosen stars
    psfAlgPolicy    = psfPolicy.get("psfPolicy")
    
    algPackage      = psfAlgPolicy.get("package")
    __import__(algPackage)
    psfAlg    = sys.modules[algPackage]
    psf, cellSet, psfSourceSet = psfAlg.getPsf(exposure, psfStars, psfCellSet, psfAlgPolicy, sdqaRatings)


    return psf, cellSet, psfSourceSet
