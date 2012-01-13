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
from .algorithmRegistry import AlgorithmRegistry
from .secondMomentStarSelector import SecondMomentStarSelector
from .algorithmsLib import SizeMagnitudeStarSelector

__all__ = ["starSelectorRegistry"]

class StarSelectorRegistry(AlgorithmRegistry):
    '''A registry of star selectors
    
    A star selector is be a class with the following API
    (Warning: this will eventually change to support prior knowledge of objects that are, or are not, stars):
    
    ConfigClass = # a Config class for this star selector
    
    def __init__(self, config):
        """Construct a star selector
        
        @param[in] config: an instance of self.ConfigClass to configure this class
        """
    
    def selectStars(self, exposure, sourceList):
        """Return a list of PSF candidates that represent likely stars
        
        The list of PSF candidates may be used by a PSF fitter to construct a PSF.
        
        @param[in] exposure: the exposure containing the sources (lsst.afw.image.Exposure)
        @param[in] sourceList: a list of sources that may be stars (lsst.afw.detection.SourceSet)
        
        @return psfCandidateList: a list of PSF candidates (each an lsst.meas.algorithms.PsfCandidate)
        """
    '''
    _requiredAttributes = ("selectStars",)

starSelectorRegistry = StarSelectorRegistry()
starSelectorRegistry.register("secondMoment", SecondMomentStarSelector)
# cannot register SizeMagnitudeStarSelector until it either has a ConfigClass attribute
# or there is a shell class (which could be defined here) that has one,
# and the constructor takes a Config.
#starSelectorRegistry.register("sizeMagnitude", SizeMagnitudeStarSelector)
