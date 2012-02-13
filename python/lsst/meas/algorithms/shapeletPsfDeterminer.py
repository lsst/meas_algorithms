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
from .algorithmsLib import ShapeletPsf

class ShapeletPsfDeterminerConfig(pexConfig.Config):
    sizeCellX = pexConfig.Field(
        doc = "Width of SpatialCell (pixels)",
        dtype = int,
#        minValue = 10,
        check = lambda x: x >- 10,
        default = 256,
    )
    sizeCellY = pexConfig.Field(
        doc = "Height of SpatialCell (pixels)",
        dtype = int,
#        minValue = 10,
        check = lambda x: x >- 10,
        default = 256,
    )
    shapeletOrder = pexConfig.Field(
        doc = "The order of the shapelet measurements",
        dtype = int,
        default = 10,
    )
    shapeletSigma = pexConfig.Field(
        doc = "The sigma to use; if <= 0 then determine from the data",
        dtype = float,
        default = -1.0,
    )
    psfAperture = pexConfig.Field(
        doc = "The aperture radius to use (arcsec)",
        dtype = float,
        default = 5.0,
    )
    nStarsPerCell = pexConfig.Field(
        doc = "The maximum number of stars to use per cell",
        dtype = int,
        default = 5,
    )
    interpOrder = pexConfig.Field(
        doc = "The order of the polynomial fit in (x,y)",
        dtype = int,
        default = 2,
    )
    interpNSigmaClip = pexConfig.Field(
        doc = "The number of sigma to use for outlier rejection",
        dtype = float,
        default = 3.0,
    )
    pcaThresh = pexConfig.Field(
        doc = "The theshold value for which principal components to keep.",
        dtype = float,
        default = 1.0e-5,
    )
    colorTerm = pexConfig.Field(
        doc = "** Warning: not implemented! ** Need some way to define what color to use.",
        dtype = str,
        default = "r-i",
    )

    def makeAlgorithm(self):
        """Make a ShapeletPsfDeterminer using the current configuration
        """
        return PcaPsfDeterminer(self)

class ShapeletPsfDeterminer(object):
    ConfigClass = ShapeletPsfDeterminerConfig
    
    def __init__(self, config):
        """Construct a PCA PSF Fitter

        @param config: an instance of ShapeletPsfDeterminerConfig
        """
        self._policy = pexConfig.makePolicy(config)

    def determinePsf(exposure, psfCandidateList, metadata=None):
        """Determine a Shapelet PSF model to an exposure given a list of PSF candidates
        
        @param[in] exposure: exposure containing the psf candidates (lsst.afw.image.Exposure)
        @param[in] psfCandidateList: a sequence of PSF candidates (each an lsst.meas.algorithms.PsfCandidate);
            typically obtained by detecting sources and then running them through a star selector
        @param[in,out] metadata: somewhere to store interesting things about the processing;
            in this case it is ignored
    
        @return
        - psf: a shapelete PSF (lsst.meas.algorithms.ShapeletPsf)
        - psfCellSet: the spacial cell set used to determine the PSF (lsst.afw.math.SpatialCellSet)
        """
        psf = ShapeletPsf(exposure, psfCandidateList, self._policy)
    
        return psf, psf.getCellSet()
