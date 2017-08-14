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
import lsst.pex.config as pexConfig
from .subtractBackground import SubtractBackgroundTask


class FindCosmicRaysConfig(pexConfig.Config):
    """Config for the findCosmicRays function
    """
    nCrPixelMax = pexConfig.Field(
        dtype=int,
        doc="maximum number of contaminated pixels",
        default=10000,
    )
    minSigma = pexConfig.Field(
        dtype=float,
        doc="CRs must be > this many sky-sig above sky",
        default=6.0,
    )
    min_DN = pexConfig.Field(
        dtype=float,
        doc="CRs must have > this many DN (== electrons/gain) in initial detection",
        default=150.0,
    )
    cond3_fac = pexConfig.Field(
        dtype=float,
        doc="used in condition 3 for CR; see CR.cc code",
        default=2.5,
    )
    cond3_fac2 = pexConfig.Field(
        dtype=float,
        doc="used in condition 3 for CR; see CR.cc code",
        default=0.6,
    )
    niteration = pexConfig.Field(
        dtype=int,
        doc="number of times to look for contaminated pixels near known CR pixels",
        default=3,
    )
    keepCRs = pexConfig.Field(
        dtype=bool,
        doc="Don't interpolate over CR pixels",
        default=False,
    )
    background = pexConfig.ConfigurableField(
        target=SubtractBackgroundTask,
        doc="Background estimation configuration"
    )

    def setDefaults(self):
        self.background.useApprox = False
        self.background.binSize = 100000
        self.background.statisticsProperty = "MEDIAN"
        self.background.undersampleStyle = "REDUCE_INTERP_ORDER"
        self.background.algorithm = "AKIMA_SPLINE"
