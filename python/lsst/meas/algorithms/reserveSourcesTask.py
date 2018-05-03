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

__all__ = ["ReserveSourcesConfig", "ReserveSourcesTask"]

import numpy as np

from lsst.pex.config import Config, Field
from lsst.pipe.base import Task, Struct


class ReserveSourcesConfig(Config):
    """Configuration for reserving sources"""
    fraction = Field(dtype=float, default=0.0,
                     doc="Fraction of candidates to reserve from fitting; none if <= 0")
    seed = Field(dtype=int, default=1,
                 doc=("This number will be added to the exposure ID to set the random seed for "
                      "reserving candidates"))


class ReserveSourcesTask(Task):
    """Reserve sources from analysis

    We randomly select a fraction of sources that will be reserved
    from analysis. This allows evaluation of the quality of model fits
    using sources that were not involved in the fitting process.

    Constructor parameters
    ----------------------
    columnName : `str`, required
        Name of flag column to add; we will suffix this with "_reserved".
    schema : `lsst.afw.table.Schema`, required
        Catalog schema.
    doc : `str`
        Documentation for column to add.
    config : `ReserveSourcesConfig`
        Configuration.
    """
    ConfigClass = ReserveSourcesConfig
    _DefaultName = "reserveSources"

    def __init__(self, columnName=None, schema=None, doc=None, **kwargs):
        Task.__init__(self, **kwargs)
        assert columnName is not None, "columnName not provided"
        assert schema is not None, "schema not provided"
        self.columnName = columnName
        self.key = schema.addField(self.columnName + "_reserved", type="Flag", doc=doc)

    def run(self, sources, prior=None, expId=0):
        """Select sources to be reserved

        Reserved sources will be flagged in the catalog, and we will return
        boolean arrays that identify the sources to be reserved from and
        used in the analysis. Typically you'll want to use the sources
        from the `use` array in your fitting, and use the sources from the
        `reserved` array as an independent test of your fitting.

        Parameters
        ----------
        sources : `lsst.afw.table.Catalog` or `list` of `lsst.afw.table.Record`
            Sources from which to select some to be reserved.
        prior : `numpy.ndarray` of type `bool`, optional
            Prior selection of sources. Should have the same length as
            `sources`. If set, we will only consider for reservation sources
            that are flagged `True` in this array.
        expId : `int`
            Exposure identifier; used for seeding the random number generator.

        Return struct contents
        ----------------------
        reserved : `numpy.ndarray` of type `bool`
            Sources to be reserved are flagged `True` in this array.
        use : `numpy.ndarray` of type `bool`
            Sources the user should use in analysis are flagged `True`.
        """
        if prior is not None:
            assert len(prior) == len(sources), "Length mismatch: %s vs %s" % (len(prior), len(sources))
            numSources = prior.sum()
        else:
            numSources = len(sources)
        selection = self.select(numSources, expId)
        if prior is not None:
            selection = self.applySelectionPrior(prior, selection)
        self.markSources(sources, selection)
        self.log.info("Reserved %d/%d sources", selection.sum(), len(selection))
        return Struct(reserved=selection,
                      use=prior & ~selection if prior is not None else np.logical_not(selection))

    def select(self, numSources, expId=0):
        """Randomly select some sources

        We return a boolean array with a random selection. The fraction
        of sources selected is specified by the config parameter `fraction`.

        Parameters
        ----------
        numSources : `int`
            Number of sources in catalog from which to select.
        expId : `int`
            Exposure identifier; used for seeding the random number generator.

        Returns
        -------
        selection : `numpy.ndarray` of type `bool`
            Selected sources are flagged `True` in this array.
        """
        selection = np.zeros(numSources, dtype=bool)
        if self.config.fraction <= 0:
            return selection
        reserve = int(np.round(numSources*self.config.fraction))
        selection[:reserve] = True
        rng = np.random.RandomState(self.config.seed + expId)
        rng.shuffle(selection)
        return selection

    def applySelectionPrior(self, prior, selection):
        """Apply selection to full catalog

        The `select` method makes a random selection of sources. If those
        sources don't represent the full population (because a sub-selection
        has already been made), then we need to generate a selection covering
        the entire population.

        Parameters
        ----------
        prior : `numpy.ndarray` of type `bool`
            Prior selection of sources, identifying the subset from which the
            random selection has been made.
        selection : `numpy.ndarray` of type `bool`
            Selection of sources in subset identified by `prior`.

        Returns
        -------
        full : `numpy.ndarray` of type `bool`
            Selection applied to full population.
        """
        full = np.zeros(len(prior), dtype=bool)
        full[prior] = selection
        return full

    def markSources(self, sources, selection):
        """Mark sources in a list or catalog

        This requires iterating through the list and setting the flag in
        each source individually. Even if the `sources` is a `Catalog`
        with contiguous records, it's not currently possible to set a boolean
        column (DM-6981) so we need to iterate.

        Parameters
        ----------
        catalog : `lsst.afw.table.Catalog` or `list` of `lsst.afw.table.Record`
            Catalog in which to flag selected sources.
        selection : `numpy.ndarray` of type `bool`
            Selection of sources to mark.
        """
        for src, select in zip(sources, selection):
            if select:
                src.set(self.key, True)
