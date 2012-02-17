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

import math
import lsst.pex.logging       as pexLog
import lsst.meas.algorithms   as measAlg
import lsst.afw.detection     as afwDetection
import lsst.afw.geom          as afwGeom
import lsst.afw.table         as afwTable
import lsst.afw.coord         as afwCoord
import lsst.afw.display.ds9   as ds9


def sourceMeasurement(
    exposure,                 # exposure to analyse
    psf,                      # psf
    footprintSets,            # sequence of (FootprintSet, isNegative) pairs
    sourceConfig,             # instance of MeasureSourcesConfig
    ):
    """ Source Measurement """

    try:
        import lsstDebug

        display = lsstDebug.Info(__name__).display
    except ImportError, e:
        try:
            display
        except NameError:
            display = False

    if display:
        frame = 0
        ds9.mtv(exposure, title="input", frame=frame)

    exposure.setPsf(psf)
    measureSources = sourceConfig.makeMeasureSources()

    sourceVector = measureSources.makeSourceVector()
    sourceConfig.slots.setupTable(sourceVector.table)
    for footprintSet, isNegative in footprintSets:
        footprintSet.makeSources(sourceVector)
    
    for i, source in enumerate(sourceVector):

        # actually try to "apply" this object; MeasureSources should swallow and log all exceptions
        # (any that escape should be considered bugs in MeasureSources)
        measureSources.apply(source, exposure)

        if display:
            ds9.dot("+", source.getX(), source.getY(), size=3, ctype=ds9.RED)
            if display > 1:
                cov = source.getCentroidErr()
                ds9.dot("@:%.1f,%.1f,%1f" % (cov[0,0], 0, cov[1,1]),  # FIXME: should '0' be cov[0,1]?
                        source.getX(), source.getY(), size=3, ctype=ds9.RED)

    return sourceVector
