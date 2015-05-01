from __future__ import absolute_import, division, print_function
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
__all__ = ["setMatchDistance"]

def setMatchDistance(matches):
    """Set the distance field of the matches in a match list to the distance in radians on the sky

    @warning the coord field of the source in each match must be correct

    @param[in,out] matches  a list of matches, an instance of lsst.afw.table.ReferenceMatch
        reads the coord field of the source and reference object of each match
        writes the distance field of each match
    """
    if len(matches) < 1:
        return

    sourceCoordKey = matches[0].first.schema["coord"].asKey()
    refObjCoordKey = matches[0].second.schema["coord"].asKey()
    for match in matches:
        sourceCoord = match.first.get(sourceCoordKey)
        refObjCoord = match.second.get(refObjCoordKey)
        match.distance = refObjCoord.angularSeparation(sourceCoord).asRadians()
