#
# LSST Data Management System
#
# Copyright 2008-2017  AURA/LSST.
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
# see <https://www.lsstcorp.org/LegalNotices/>.
#
__all__ = ["MakePsfCandidatesConfig", "MakePsfCandidatesTask"]

import numpy as np

from lsst.afw.table import SourceCatalog
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pex.exceptions
import lsst.pipe.base as pipeBase
from . import makePsfCandidate


class MakePsfCandidatesConfig(pexConfig.Config):
    kernelSize = pexConfig.Field(
        doc="size of the kernel to create",
        dtype=int,
        default=21,
    )
    borderWidth = pexConfig.Field(
        doc="number of pixels to ignore around the edge of PSF candidate postage stamps",
        dtype=int,
        default=0,
    )


class MakePsfCandidatesTask(pipeBase.Task):
    """Create PSF candidates given an input catalog.
    """
    ConfigClass = MakePsfCandidatesConfig
    _DefaultName = "makePsfCandidates"

    def __init__(self, **kwds):
        pipeBase.Task.__init__(self, **kwds)

    def run(self, starCat, exposure, psfCandidateField=None):
        """Make a list of PSF candidates from a star catalog.

        Parameters
        ----------
        starCat : `lsst.afw.table.SourceCatalog`
            Catalog of stars, as returned by
            ``lsst.meas.algorithms.starSelector.run()``.
        exposure : `lsst.afw.image.Exposure`
            The exposure containing the sources.
        psfCandidateField : `str` or None
            Name of flag field to set True for PSF candidates, or None to not
            set a field; the field is left unchanged for non-candidates.

        Returns
        -------
        struct : `lsst.pipe.base.Struct`
            Results struct containing:

            - ``psfCandidates`` : List of PSF candidates
                (`list` of `lsst.meas.algorithms.PsfCandidate`).
            - ``goodStarCat`` : Subset of ``starCat`` that was successfully made
                into PSF candidates (`lsst.afw.table.SourceCatalog`).

        """
        psfResult = self.makePsfCandidates(starCat, exposure)

        if psfCandidateField is not None:
            isStarKey = starCat.schema[psfCandidateField].asKey()
            for star in psfResult.goodStarCat:
                star.set(isStarKey, True)

        return psfResult

    def makePsfCandidates(self, starCat, exposure):
        """Make a list of PSF candidates from a star catalog.

        Parameters
        ----------
        starCat : `lsst.afw.table.SourceCatalog`
            Catalog of stars, as returned by
            ``lsst.meas.algorithms.starSelector.run()``.
        exposure : `lsst.afw.image.Exposure`
            The exposure containing the sources.

        Returns
        -------
        struct : `lsst.pipe.base.Struct`
            Results struct containing:

            - ``psfCandidates`` : List of PSF candidates
                (`list` of `lsst.meas.algorithms.PsfCandidate`).
            - ``goodStarCat`` : Subset of ``starCat`` that was successfully made
                into PSF candidates (`lsst.afw.table.SourceCatalog`).
        """
        goodStarCat = SourceCatalog(starCat.schema)

        psfCandidateList = []
        didSetSize = False
        for star in starCat:
            try:
                psfCandidate = makePsfCandidate(star, exposure)

                # The setXXX methods are class static, but it's convenient to call them on
                # an instance as we don't know exposures's pixel type
                # (and hence psfCandidate's exact type)
                if not didSetSize:
                    psfCandidate.setBorderWidth(self.config.borderWidth)
                    psfCandidate.setWidth(self.config.kernelSize + 2*self.config.borderWidth)
                    psfCandidate.setHeight(self.config.kernelSize + 2*self.config.borderWidth)
                    didSetSize = True

                im = psfCandidate.getMaskedImage().getImage()
            except lsst.pex.exceptions.LengthError as err:
                self.log.debug("Failed to make a psfCandidate from star %d: %s", star.getId(), err)
                continue
            except lsst.pex.exceptions.Exception as err:
                self.log.warn("Failed to make a psfCandidate from star %d: %s", star.getId(), err)
                continue

            vmax = afwMath.makeStatistics(im, afwMath.MAX).getValue()
            if not np.isfinite(vmax):
                continue
            psfCandidateList.append(psfCandidate)
            goodStarCat.append(star)

        return pipeBase.Struct(
            psfCandidates=psfCandidateList,
            goodStarCat=goodStarCat,
        )
