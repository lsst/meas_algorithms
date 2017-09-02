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
__all__ = ["makeCoaddApCorrMap", ]

from lsst.afw.image import ApCorrMap
from .coaddBoundedField import CoaddBoundedField, CoaddBoundedFieldElement


def makeCoaddApCorrMap(catalog, coaddBox, coaddWcs, weightFieldName="weight"):
    """Construct an ApCorrMap for a coadd

    @param catalog: Table of coadd inputs (lsst.afw.table.ExposureCatalog)
    @param coaddBox: Bounding box for coadd (lsst.afw.geom.Box2I)
    @param coaddWcs: Wcs for coadd
    @param weightFieldName: name of weight field in catalog
    @return aperture corrections
    """

    # Assemble the BoundedFields for each type
    everything = {}  # name --> list of CoaddBoundedFieldElement
    weightKey = catalog.schema[weightFieldName].asKey()
    for row in catalog:
        apCorrMap = row.getApCorrMap()
        if not apCorrMap:
            continue
        weight = row.get(weightKey)
        wcs = row.getWcs()
        validPolygon = row.getValidPolygon()
        for name, bf in apCorrMap.items():
            if name not in everything:
                everything[name] = []
            everything[name].append(CoaddBoundedFieldElement(bf, wcs, validPolygon, weight))

    # Construct a CoaddBoundedField for each type
    apCorrMap = ApCorrMap()
    for name, elements in everything.items():
        apCorrMap.set(name, CoaddBoundedField(coaddBox, coaddWcs, elements))

    return apCorrMap
