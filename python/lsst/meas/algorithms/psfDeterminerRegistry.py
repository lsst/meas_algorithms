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
'''
A registry of PSF determiners

A PSF determiner should be a class with the following API:

def __init__(self, policy):
    """Construct a PSF Determiner

    @param[in] policy: see <path to policy dictionary>
    """

def determinePsf(exposure, psfCandidateList, sdqaRatings=None):
    """Determine a PSF model
    
    @param[in] exposure: exposure containing the psf candidates (lsst.afw.image.Exposure)
    @param[in] psfCandidateList: a sequence of PSF candidates (each an lsst.meas.algorithms.PsfCandidate);
        typically obtained by detecting sources and then running them through a star selector
    @param[in,out] sdqaRatingSet: an lsst.sdqa.SdqaRatingSet(); if not None it may gain relevant SDQA ratings

    @return
    - psf: the fit PSF; a subclass of lsst.afw.detection.Psf
    - cellSet: the spatial cell set used to determine the PSF (lsst.afw.math.SpatialCellSet)
    """

'''
import shapeletPsfDeterminer
import pcaPsfDeterminer

_psfDeterminerRegistry = dict(
    pcaPsfDeterminer = pcaPsfDeterminer.PcaPsfDeterminer,
    shapeletPsfDeterminer = shapeletPsfDeterminer.ShapeletPsfDeterminer,
)

def getPsfDeterminerNames():
    """Return a list of names of PSF determiners
    """
    global _psfDeterminerRegistry
    return _psfDeterminerRegistry.keys()

def makePsfDeterminer(name, policy):
    """Construct a PSF determiner from its name and policy

    @param[in] name: name of PSF determiner
    @param[in] policy: policy for PSF determiner
    @raise KeyError if the name is unrecognized
    """
    global _psfDeterminerRegistry
    return _psfDeterminerRegistry[name](policy)

def registerPsfDeterminer(name, psfDeterminer):
    """Register a new PSF determiner. The name must be globally unique.
    
    To avoid issues, you may wish to include the name of your package in the determiner name, e.g.:
    "my_pkg.myPsfDeterminer"
    
    @raise KeyError if the name is a duplicate
    """
    global _psfDeterminerRegistry
    if name in _psfDeterminerRegistry:
        raise KeyError("A PSF determiner already exists with name %r" % (name,))
    _psfDeterminerRegistry[name] = psfDeterminer
