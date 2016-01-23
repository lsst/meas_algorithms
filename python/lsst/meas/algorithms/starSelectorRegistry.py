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
from lsst.pex.config import makeRegistry
from .secondMomentStarSelector import SecondMomentStarSelector
from .sizeMagnitudeStarSelectorFactory import sizeMagnitudeStarSelectorFactory

__all__ = ["starSelectorRegistry"]

starSelectorRegistry = makeRegistry(
    '''A registry of star selector factories
    
        A star selector factory makes a class with the following API:
        
        def __init__(self, config):
            """Construct a star selector
            
            @param[in] config: an instance of pexConfig.Config that configures this algorithm
            """
        
        def selectStars(self, exposure, sourceCat, matches=None):
            """Return a list of PSF candidates that represent likely stars
            
            The list of PSF candidates may be used by a PSF fitter to construct a PSF.
            
            @param[in] exposure  the exposure containing the sources (an lsst.afw.image.Exposure)
            @param[in] sourceCat catalog of sources that may be stars (an lsst.afw.table.SourceCatalog)
            @param[in] matches  list of reference object/source matches
                (an lsst.afw.table.ReferenceMatchVector), or None. Some star selectors
                will ignore this argument, others may require it.
            
            @return psfCandidateList: a list of PSF candidates (each an lsst.meas.algorithms.PsfCandidate)
            """
    '''
)

starSelectorRegistry.register("secondMoment", SecondMomentStarSelector)
starSelectorRegistry.register("sizeMagnitude", sizeMagnitudeStarSelectorFactory)
