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
import lsst.pex.config as pexConfig
from .secondMomentStarSelector import SecondMomentStarSelectorConfig
from .sizeMagnitudeStarSelectorConfig import SizeMagnitudeStarSelectorConfig

__all__ = ["starSelectorRegistry"]

starSelectorRegistry = pexConfig.makeConfigRegistry(
    doc = '''A registry of star selector configs
    
        A star selector config is a subclass of pexConfig.Config that has the following API:
        
        def makeAlgorithm(self):
            """Make a star selector using the current config
            """
        
        A star selector is a class with the following API:
        
        def __init__(self, config):
            """Construct a star selector
            
            @param[in] config: an instance of pexConfig.Config that configures this algorithm
            """
        
        def selectStars(self, exposure, sourceList):
            """Return a list of PSF candidates that represent likely stars
            
            The list of PSF candidates may be used by a PSF fitter to construct a PSF.
            
            @param[in] exposure: the exposure containing the sources (lsst.afw.image.Exposure)
            @param[in] sourceList: a list of sources that may be stars (lsst.afw.detection.SourceSet)
            
            @return psfCandidateList: a list of PSF candidates (each an lsst.meas.algorithms.PsfCandidate)
            """
        ''',
)

starSelectorRegistry["secondMoment"] = SecondMomentStarSelectorConfig
starSelectorRegistry["sizeMagnitude"] = SizeMagnitudeStarSelectorConfig
