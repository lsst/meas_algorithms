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
import algorithmsLib

# this is the sort of thing wanted, but I don't know how to get it down into the C++ code
# and C++ has no support for constructing Config in C++ itself
# class ShapeletPsfDeterminerConfig(pexConfig.Config):
#     sizeCellX = pexConfig.Field(
#         dtype = int,
#         doc = "Width of SpatialCell (pixels)",
#         minValue = 10,
#         default = 256,
#     )
#     sizeCellY = pexConfig.Field(
#         dtype = int,
#         doc = "Height of SpatialCell (pixels)",
#         minValue = 10,
#         default = 256,
#     )
#     shapeletOrder = pexConfig.Field(
#         dtype = int,
#         doc = "The order of the shapelet measurements",
#         default = 10,
#     )
#     shapeletSigma = pexConfig.Field(
#         dtype = double,
#         doc = "The sigma to use; if <= 0 then determine from the data",
#         default = -1.0,
#     )
#     psfAperture = pexConfig.Field(
#         dtype = double,
#         doc = "The aperture radius to use (arcsec)",
#         default = 5.0,
#     )
#     nStarsPerCell = pexConfig.Field(
#         dtype = int,
#         doc = "The maximum number of stars to use per cell",
#         default = 5,
#     )
#     interpOrder = pexConfig.Field(
#         dtype = int,
#         doc = "The order of the polynomial fit in (x,y)",
#         default = 2,
#     )
#     interpNSigmaClip = pexConfig.Field(
#         dtype = double,
#         doc = "The number of sigma to use for outlier rejection",
#         default = 3.0,
#     )
#     pcaThresh = pexConfig.Field(
#         dtype = double,
#         doc = "The theshold value for which principal components to keep.",
#         default = 1.0e-5,
#     )
#     colorTerm = pexConfig.Field(
#         dtype = string,
#         doc = "** Warning: not implemented! ** Need some way to define what color to use.",
#         default = "r-i",
#     )


class ShapeletPsfDeterminer(object):
    def __init__(self, policy):
        """Construct a PCA PSF Fitter

        @param policy: see policy/ShapeletPsfDeterminerDictionary.paf
        """
        self._policy = policy

    def determinePsf(exposure, psfCandidateList, metadata=None):
        """Determine a Shapelet PSF model to an exposure given a list of PSF candidates
        
        @param[in] exposure: exposure containing the psf candidates (lsst.afw.image.Exposure)
        @param[in] psfCandidateList: a sequence of PSF candidates (each an lsst.meas.algorithms.PsfCandidate);
            typically obtained by detecting sources and then running them through a star selector
        @param[in,out] metadata: somewhere to store interesting things about the processing
    
        @return psf: a shapelete PSF (lsst.meas.algorithms.ShapeletPsf)
        """
        return algorithmsLib.ShapeletPsf(exposure, psfCandidateList, self._policy)
    
        return psf, psf.getCellSet()
