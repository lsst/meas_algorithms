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
__all__ = ["BaseStarSelectorConfig", "BaseStarSelectorTask", "starSelectorRegistry"]

import abc

import numpy as np

from lsst.afw.table import SourceCatalog, Schema
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from . import makePsfCandidate


class BaseStarSelectorConfig(pexConfig.Config):
    badFlags = pexConfig.ListField(
        doc="List of flags which cause a source to be rejected as bad",
        dtype=str,
        default=[
            "base_PixelFlags_flag_edge",
            "base_PixelFlags_flag_interpolatedCenter",
            "base_PixelFlags_flag_saturatedCenter",
            "base_PixelFlags_flag_crCenter",
            "base_PixelFlags_flag_bad",
        ],
    )


class BaseStarSelectorTask(pipeBase.Task, metaclass=abc.ABCMeta):
    """!Base class for star selectors

    Register all star selectors with the starSelectorRegistry using:
        starSelectorRegistry.register(name, class)
    """

    usesMatches = False  # Does the star selector use the "matches" argument in the "run method? Few do.
    ConfigClass = BaseStarSelectorConfig
    _DefaultName = "starSelector"

    def __init__(self, schema, **kwds):
        # catch code that passed config positionally before schema argument was added
        assert isinstance(schema, Schema)
        pipeBase.Task.__init__(self, **kwds)

    def run(self, exposure, sourceCat, matches=None, isStarField=None):
        """Select stars and set a flag field True for stars in the input catalog.

        Parameters
        ----------
        exposure :
            the exposure containing the sources
        sourceCat :
            catalog of sources that may be stars (an lsst.afw.table.SourceCatalog)
        matches :
            astrometric matches; ignored by this star selector
                (an lsst.afw.table.ReferenceMatchVector), or None. Some star selectors
                will ignore this argument, others may require it. See the usesMatches class variable.
        isStarField :
            name of flag field to set True for stars, or None to not set a field;
            the field is left unchanged for non-stars

        Returns
        -------
        struct : `lsst.pipe.base.Struct`
        Result struct containing:

            - starCat  catalog of stars that were selected as stars and successfuly made into PSF candidates
                        (a subset of sourceCat whose records are shallow copies)
            - psfCandidates  list of PSF candidates (lsst.meas.algorithms.PsfCandidate)
        """
        result = self.selectStars(exposure=exposure, sourceCat=sourceCat, matches=matches)

        if isStarField is not None:
            isStarKey = sourceCat.schema[isStarField].asKey()
            for star in result.starCat:
                star.set(isStarKey, True)

        return pipeBase.Struct(starCat=result.starCat)

    @abc.abstractmethod
    def selectStars(self, exposure, sourceCat, matches=None):
        """!Return a catalog of stars: a subset of sourceCat whose records are shallow copies

        @param[in] exposure  the exposure containing the sources
        @param[in] sourceCat  catalog of sources that may be stars (an lsst.afw.table.SourceCatalog)
        @param[in] matches  astrometric matches; ignored by this star selector
                (an lsst.afw.table.ReferenceMatchVector), or None. Some star selectors
                will ignore this argument, others may require it. See the usesMatches class variable.

        @warning The returned catalog must have records that are shallow copies
        (fortunately this is the default behavior when you add a record from one catalog to another);
        otherwise the run method cannot set the isStarField flag in the original source catalog.

        @return a pipeBase.Struct containing:
        - starCat  a catalog of stars (a subset of sourceCat whose records are shallow copies)
        """
        raise NotImplementedError("BaseStarSelectorTask is abstract, subclasses must override this method")


starSelectorRegistry = pexConfig.makeRegistry(
    doc="A registry of star selectors (subclasses of BaseStarSelectorTask)",
)
