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
from .shapeletPsfDeterminer import ShapeletPsfDeterminerConfig
from .pcaPsfDeterminer import PcaPsfDeterminerConfig

__all__ = ["psfDeterminerRegistry"]

psfDeterminerRegistry = pexConfig.makeConfigRegistry(
    doc = '''A registry of PSF determiner configs
    
        A PSF determiner config is a subclass of pexConfig.Config with the following API:
        
        def makeAlgorithm():
            """Construct a PSF determiner using the current config
            """
    
        A PSF determiner is a class with the following API:
        
        def __init__(self, config):
            """Construct a PSF Determiner
        
            @param[in] config: an instance of pexConfig.Config that configures this algorithm
            """
        
        def determinePsf(exposure, psfCandidateList, metadata=None):
            """Determine a PSF model
            
            @param[in] exposure: exposure containing the psf candidates (lsst.afw.image.Exposure)
            @param[in] psfCandidateList: a sequence of PSF candidates (each an lsst.meas.algorithms.PsfCandidate);
                typically obtained by detecting sources and then running them through a star selector
            @param[in,out] metadata: a place to save interesting items
        
            @return
            - psf: the fit PSF; a subclass of lsst.afw.detection.Psf
            - cellSet: the spatial cell set used to determine the PSF (lsst.afw.math.SpatialCellSet)
            """
        ''',
)

psfDeterminerRegistry["pca"] = PcaPsfDeterminerConfig
psfDeterminerRegistry["shapelet"] = ShapeletPsfDeterminerConfig
