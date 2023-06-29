# This file is part of meas_algorithms.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

__all__ = ["BasePsfDeterminerConfig", "BasePsfDeterminerTask", "psfDeterminerRegistry"]

import abc

import lsst.pipe.base as pipeBase
import lsst.pex.config as pexConfig


class BasePsfDeterminerConfig(pexConfig.Config):
    """Configuration that is likely to be shared by all PSF determiners

    This is fairly sparse; more fields can be moved here once it is clear they are universal.
    """
    stampSize = pexConfig.Field[int](
        doc="Size of the postage stamp (in native pixels) to render the PSF model. Should be odd.",
        default=None,
        optional=True,
        check=lambda x: (x > 0) & (x % 2 == 1),
    )
    kernelSize = pexConfig.Field[int](
        doc="Size of the kernel to create. "
            "If less than 15, then the median of the square root of the "
            "stellar quadrupole moments is multiplied by kernelSize and "
            "used as the kernel size.",
        default=None,
        optional=True,
        deprecated="This field is deprecated and will be removed after v25. "
                   "Use stampSize instead.",
    )
    kernelSizeMin = pexConfig.Field[int](
        doc="Minimum radius of the kernel. Relevant only if kernelSize < 15.",
        default=25,
        deprecated="This field is deprecated and will be removed after v25.",
    )
    kernelSizeMax = pexConfig.Field[int](
        doc="Maximum radius of the kernel. Relevant only if kernelSize < 15.",
        default=45,
        deprecated="This field is deprecated and will be removed after v25.",
    )

    def validate(self):
        # TODO: DM-36311: Lines involving kernelSize* (or the entire method)
        # may be removed after v25.
        super().validate()
        # Ensure that stampSize and kernelSize are not contradictory.
        if self.stampSize is not None and self.kernelSize is not None:
            raise pexConfig.FieldValidationError(self.__class__.kernelSize, self,
                                                 "Only one of stampSize (preferable) and kernelSize "
                                                 "(deprecated) must be set."
                                                 )
        if self.kernelSizeMin > self.kernelSizeMax:
            raise pexConfig.FieldValidationError(self.__class__.kernelSizeMax, self,
                                                 "kernelSizeMin must be <= kernelSizeMax"
                                                 )


class BasePsfDeterminerTask(pipeBase.Task, metaclass=abc.ABCMeta):
    """Base class for PSF determiners

    Register all PSF determiners with the psfDeterminerRegistry using:
        psfDeterminerRegistry.register(name, class)

    Parameters
    ----------
    config : `lsst.pexConfig.Config`
        Input for configuring the algorithm
    schema : `lsst.afw.table.Schema`
        Schema used for sources; passing a schema allows the
        determiner to reserve a flag field to mark stars used in
        PSF measurement, but some PSF determiners ignore this argument.
    """

    usesMatches = False  # Does the PSF determiner use the "matches" argument in the "run method? Few do.
    ConfigClass = BasePsfDeterminerConfig
    _DefaultName = "psfDeterminer"

    def __init__(self, config, schema=None, **kwds):
        pipeBase.Task.__init__(self, config=config, **kwds)

    @abc.abstractmethod
    def determinePsf(self, exposure, psfCandidateList, metadata=None, flagKey=None):
        """Determine a PSF model.

        Parameters
        ----------
        exposure : `lsst.afw.Exposure`
            Exposure containing the psf candidates.
        psdCandidateList : `list` [`lsst.meas.algorithms.PsfCandidate`]
            A sequence of PSF candidates; typically obtained by
            detecting sources and then running them through a star
            selector.
        metadata : `str`, optional
            A place to save interesting items.
        flagKey: `lsst.afw.table.Key`, optional
            Schema key used to mark sources actually used in PSF determination.

        Returns
        -------
        psf : `lsst.afw.detection.Psf`
            The fit PSF.
        cellSet : `lsst.afw.math.SpatialCellSet`
            The spatial cell set used to determine the PSF
        """
        raise NotImplementedError("BasePsfDeterminerTask is abstract, subclasses must override this method")


psfDeterminerRegistry = pexConfig.makeRegistry(
    doc="A registry of PSF determiners (subclasses of BasePsfDeterminerTask)",
)
