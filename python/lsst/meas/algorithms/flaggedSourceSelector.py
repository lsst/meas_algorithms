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

"""Select sources that have an existing flag field set."""
import lsst.pex.config
import lsst.afw.table
import lsst.pipe.base as pipeBase

from .sourceSelector import BaseSourceSelectorTask, sourceSelectorRegistry

__all__ = ["FlaggedSourceSelectorConfig", "FlaggedSourceSelectorTask"]


class FlaggedSourceSelectorConfig(BaseSourceSelectorTask.ConfigClass):
    field = lsst.pex.config.Field(
        dtype=str, default="calib_psf_used",
        doc="Name of a flag field that is True for Sources that should be used.",
    )


@lsst.pex.config.registerConfigurable("flagged", sourceSelectorRegistry)
class FlaggedSourceSelectorTask(BaseSourceSelectorTask):
    """
    A trivial SourceSelector that simply uses an existing flag field to filter
    a SourceCatalog.

    This is most frequently used in steps that occur after the a PSF model has
    been built, to allow other procedures that need Sources to use the set of
    Sources used to determine the PSF.

    Attributes
    ----------
    usesMatches : `bool`
        A boolean variable specify if the inherited source selector uses
        matches.
    """

    ConfigClass = FlaggedSourceSelectorConfig
    _DefaultName = "flagged"

    def selectSources(self, sourceCat, matches=None, exposure=None):
        """Return a bool array representing which sources to select from
        sourceCat.

        The input catalog must be contiguous in memory.

        Parameters
        ----------
        sourceCat : `lsst.afw.table.SourceCatalog`
            Catalog of sources to select from.
        matches : `list` of `lsst.afw.table.ReferenceMatch` or None
            Ignored in this SourceSelector.
        exposure : `lsst.afw.image.Exposure` or None
            The exposure the catalog was built from; used for debug display.

        Returns
        -------
        struct : `lsst.pipe.base.Struct`
            The struct contains the following data:

            ``selected``
                Boolean array of sources that were selected, same length as
                ``sourceCat``. (`numpy.ndarray` of `bool`)
        """
        key = sourceCat.schema.find(self.config.field).key
        return pipeBase.Struct(selected=sourceCat[key])
