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
'''A registry for star selectors

A star selector should be a class with the following API:

def __init__(self, policy):
    """Construct a star selector
    
    @param[in] policy: star selection policy; see <path to dictionary>
    """

def selectStars(self, exposure, sourceList):
    """Return a list of PSF candidates that represent likely stars
    
    The list of PSF candidates may be used by a PSF fitter to construct a PSF.
    
    @param[in] exposure: the exposure containing the sources (lsst.afw.image.Exposure)
    @param[in] sourceList: a list of sources that may be stars (lsst.afw.detection.SourceSet)
    
    @return psfCandidateList: a list of PSF candidates (each an lsst.meas.algorithms.PsfCandidate)
    """
'''

import secondMomentStarSelector
import algorithmsLib

_starSelectorRegistry = dict(
    secondMomentStarSelector = secondMomentStarSelector.SecondMomentStarSelector,
    sizeMagnitudeStarSelector = algorithmsLib.SizeMagnitudeStarSelector,
)

def getStarSelectorNames():
    """Return the names of all registered star selectors
    """
    global _starSelectorRegistry
    return _starSelectorRegistry.keys()

def makeStarSelector(name, policy):
    """Construct a star selector given its name and policy
    
    @param[in] name: name of star selector
    @param[in] policy: policy for star selector
    @raise KeyError if the name is unrecognized
    """
    global _starSelectorRegistry
    return _starSelectorRegistry[name](policy)

def registerStarSelector(name, starSelector):
    '''Register a new star selector.
    
    @param[in] name: name of star selector. The name must be globally unique. To avoid conficts,
        please use the package name as a prefix for selectors that aren't built in e.g. my_pkg.myStarSelector
    @param[in] starSelector: a star selector (see file doc string).
    
    @raise KeyError if the name is a duplicate
    '''
    global _starSelectorRegistry
    if name in _starSelectorRegistry:
        raise KeyError("A star selector already exists with name %r" % (name,))
    _starSelectorRegistry[name] = starSelector
