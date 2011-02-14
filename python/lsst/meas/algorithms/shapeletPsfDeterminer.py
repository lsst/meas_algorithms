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

class ShapeletPsfDeterminer(object):
    def __init__(self, policy):
        """Construct a PCA PSF Fitter

        @param policy: see policy/ShapeletPsfDeterminerDictionary.paf
        """
        self._policy = policy

    def determinePsf(exposure, psfCandidateList, sdqaRatingSet=None):
        """Determine a Shapelet PSF model to an exposure given a list of PSF candidates
        
        @param[in] exposure: exposure containing the psf candidates (lsst.afw.image.Exposure)
        @param[in] psfCandidateList: a sequence of PSF candidates (each an lsst.meas.algorithms.PsfCandidate);
            typically obtained by detecting sources and then running them through a star selector
        @param[in,out] sdqaRatingSet: an lsst.sdqa.SdqaRatingSet()
    
        @return psf: a shapelete PSF (lsst.meas.algorithms.ShapeletPsf)
        """
        return algorithmsLib.ShapeletPsf(exposure, psfCandidateList, self._policy)
    
        return psf, psf.getCellSet()
