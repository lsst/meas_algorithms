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

__all__ = ["getRefFluxField", "getRefFluxKeys", "LoadReferenceObjectsTask", "LoadReferenceObjectsConfig"]

import abc

import numpy

import lsst.geom
import lsst.afw.table as afwTable
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.daf.base import PropertyList


def getRefFluxField(schema, filterName=None):
    """!Get name of flux field in schema

    if filterName is specified:
        return *filterName*_camFlux if present
        else return *filterName*_flux if present (camera filter name matches reference filter name)
        else throw RuntimeError
    else:
        return camFlux, if present,
        else throw RuntimeError

    @param[in] schema  reference catalog schema
    @param[in] filterName  name of camera filter
    @return flux field name
    @throw RuntimeError if appropriate field is not found
    """
    if not isinstance(schema, afwTable.Schema):
        raise RuntimeError("schema=%s is not a schema" % (schema,))
    if filterName:
        fluxFieldList = [filterName + "_camFlux", filterName + "_flux"]
    else:
        fluxFieldList = ["camFlux"]
    for fluxField in fluxFieldList:
        if fluxField in schema:
            return fluxField

    raise RuntimeError("Could not find flux field(s) %s" % (", ".join(fluxFieldList)))


def getRefFluxKeys(schema, filterName=None):
    """!Return flux and flux error keys

    @param[in] schema  reference catalog schema
    @param[in] filterName  name of camera filter
    @return a pair of keys:
        flux key
        flux error key, if present, else None
    @throw RuntimeError if flux field not found
    """
    fluxField = getRefFluxField(schema, filterName)
    fluxErrField = fluxField + "Err"
    fluxKey = schema[fluxField].asKey()
    try:
        fluxErrKey = schema[fluxErrField].asKey()
    except Exception:
        fluxErrKey = None
    return (fluxKey, fluxErrKey)


class LoadReferenceObjectsConfig(pexConfig.Config):
    pixelMargin = pexConfig.RangeField(
        doc="Padding to add to 4 all edges of the bounding box (pixels)",
        dtype=int,
        default=300,
        min=0,
    )
    defaultFilter = pexConfig.Field(
        doc="Default reference catalog filter to use if filter not specified in exposure; " +
        "if blank then filter must be specified in exposure",
        dtype=str,
        default="",
    )
    filterMap = pexConfig.DictField(
        doc="Mapping of camera filter name: reference catalog filter name; " +
        "each reference filter must exist",
        keytype=str,
        itemtype=str,
        default={},
    )

# The following comment block adds a link to this task from the Task Documentation page.
## @addtogroup LSST_task_documentation
## @{
## @page LoadReferenceObjectsTask
## @ref LoadReferenceObjectsTask_ "LoadReferenceObjectsTask"
## @copybrief LoadReferenceObjectsTask
## @}


class LoadReferenceObjectsTask(pipeBase.Task, metaclass=abc.ABCMeta):
    """!Abstract base class to load objects from reference catalogs

    @anchor LoadReferenceObjectsTask_

    @section meas_algorithms_loadReferenceObjects_Contents Contents

     - @ref meas_algorithms_loadReferenceObjects_Purpose
     - @ref meas_algorithms_loadReferenceObjects_Initialize
     - @ref meas_algorithms_loadReferenceObjects_IO
     - @ref meas_algorithms_loadReferenceObjects_Schema
     - @ref meas_algorithms_loadReferenceObjects_Config

    @section meas_algorithms_loadReferenceObjects_Purpose  Description

    Abstract base class for tasks that load objects from a reference catalog
    in a particular region of the sky.

    Implementations must subclass this class, override the loadSkyCircle method,
    and will typically override the value of ConfigClass with a task-specific config class.

    @section meas_algorithms_loadReferenceObjects_Initialize   Task initialisation

    @copydoc \_\_init\_\_

    @section meas_algorithms_loadReferenceObjects_IO       Invoking the Task

    @copydoc loadPixelBox

    @section meas_algorithms_loadReferenceObjects_Schema       Schema of the reference object catalog

    Reference object catalogs are instances of lsst.afw.table.SimpleCatalog with the following schema
    (other fields may also be present):
    - coord: ICRS position of star on sky (an lsst.geom.SpherePoint)
    - centroid: position of star on an exposure, if relevant (an lsst.afw.Point2D)
    - hasCentroid: is centroid usable?
    - *referenceFilterName*_flux: brightness in the specified reference catalog filter (Jy)
        Note: the function lsst.afw.image.abMagFromFlux will convert flux in Jy to AB Magnitude.
    - *referenceFilterName*_fluxErr (optional): brightness standard deviation (Jy);
        omitted if no data is available; possibly nan if data is available for some objects but not others
    - camFlux: brightness in default camera filter (Jy); omitted if defaultFilter not specified
    - camFluxErr: brightness standard deviation for default camera filter;
        omitted if defaultFilter not specified or standard deviation not available that filter
    - *cameraFilterName*_camFlux: brightness in specified camera filter (Jy)
    - *cameraFilterName*_camFluxErr (optional): brightness standard deviation
        in specified camera filter (Jy); omitted if no data is available;
        possibly nan if data is available for some objects but not others
    - photometric (optional): is the object usable for photometric calibration?
    - resolved (optional): is the object spatially resolved?
    - variable (optional): does the object have variable brightness?

    @section meas_algorithms_loadReferenceObjects_Config       Configuration parameters

    See @ref LoadReferenceObjectsConfig for a base set of configuration parameters.
    Most subclasses will add configuration variables.
    """
    ConfigClass = LoadReferenceObjectsConfig
    _DefaultName = "LoadReferenceObjects"

    def __init__(self, butler=None, *args, **kwargs):
        """!Construct a LoadReferenceObjectsTask

        @param[in] butler  A daf.persistence.Butler object.  This allows subclasses to use the butler to
        access reference catalog files using the stack I/O abstraction scheme.
        """
        pipeBase.Task.__init__(self, *args, **kwargs)
        self.butler = butler

    @pipeBase.timeMethod
    def loadPixelBox(self, bbox, wcs, filterName=None, calib=None):
        """!Load reference objects that overlap a pixel-based rectangular region

        The search algorithm works by searching in a region in sky coordinates whose center is the center
        of the bbox and radius is large enough to just include all 4 corners of the bbox.
        Stars that lie outside the bbox are then trimmed from the list.

        @param[in] bbox  bounding box for pixels (an lsst.geom.Box2I or Box2D)
        @param[in] wcs  WCS (an lsst.afw.geom.SkyWcs)
        @param[in] filterName  name of camera filter, or None or blank for the default filter
        @param[in] calib  calibration, or None if unknown

        @return an lsst.pipe.base.Struct containing:
        - refCat a catalog of reference objects with the
            @link meas_algorithms_loadReferenceObjects_Schema standard schema @endlink
            as documented in LoadReferenceObjects, including photometric, resolved and variable;
            hasCentroid is False for all objects.
        - fluxField = name of flux field for specified filterName
        """
        circle = self._calculateCircle(bbox, wcs)

        # find objects in circle
        self.log.info("Loading reference objects using center %s and radius %s deg" %
                      (circle.coord, circle.radius.asDegrees()))
        loadRes = self.loadSkyCircle(circle.coord, circle.radius, filterName)
        refCat = loadRes.refCat
        numFound = len(refCat)

        # trim objects outside bbox
        refCat = self._trimToBBox(refCat=refCat, bbox=circle.bbox, wcs=wcs)
        numTrimmed = numFound - len(refCat)
        self.log.debug("trimmed %d out-of-bbox objects, leaving %d", numTrimmed, len(refCat))
        self.log.info("Loaded %d reference objects", len(refCat))

        # make sure catalog is contiguous
        if not refCat.isContiguous():
            loadRes.refCat = refCat.copy(deep=True)

        return loadRes

    @abc.abstractmethod
    def loadSkyCircle(self, ctrCoord, radius, filterName=None):
        """!Load reference objects that overlap a circular sky region

        @param[in] ctrCoord  ICRS center of search region (an lsst.geom.SpherePoint)
        @param[in] radius  radius of search region (an lsst.geom.Angle)
        @param[in] filterName  name of filter, or None for the default filter;
            used for flux values in case we have flux limits (which are not yet implemented)

        @return an lsst.pipe.base.Struct containing:
        - refCat a catalog of reference objects with the
            @link meas_algorithms_loadReferenceObjects_Schema standard schema @endlink
            as documented in LoadReferenceObjects, including photometric, resolved and variable;
            hasCentroid is False for all objects.
        - fluxField = name of flux field for specified filterName
        """
        return

    @staticmethod
    def _trimToBBox(refCat, bbox, wcs):
        """!Remove objects outside a given pixel-based bbox and set centroid and hasCentroid fields

        @param[in,out] refCat  a catalog of objects (an lsst.afw.table.SimpleCatalog,
            or other table type that has fields "coord", "centroid" and "hasCentroid").
            The "coord" field is read.
            The "centroid" and "hasCentroid" fields are set.
        @param[in] bbox  pixel region (an afwImage.Box2D)
        @param[in] wcs  WCS used to convert sky position to pixel position (an lsst.afw.math.WCS)

        @return a catalog of reference objects in bbox, with centroid and hasCentroid fields set
        """
        afwTable.updateRefCentroids(wcs, refCat)
        centroidKey = afwTable.Point2DKey(refCat.schema["centroid"])
        retStarCat = type(refCat)(refCat.table)
        for star in refCat:
            point = star.get(centroidKey)
            if bbox.contains(point):
                retStarCat.append(star)
        return retStarCat

    def _addFluxAliases(self, schema):
        """Add aliases for camera filter fluxes to the schema

        If self.config.defaultFilter then adds these aliases:
            camFlux:      <defaultFilter>_flux
            camFluxErr: <defaultFilter>_fluxErr, if the latter exists

        For each camFilter: refFilter in self.config.filterMap adds these aliases:
            <camFilter>_camFlux:      <refFilter>_flux
            <camFilter>_camFluxErr: <refFilter>_fluxErr, if the latter exists

        @throw RuntimeError if any reference flux field is missing from the schema
        """
        aliasMap = schema.getAliasMap()

        def addAliasesForOneFilter(filterName, refFilterName):
            """Add aliases for a single filter

            @param[in] filterName  camera filter name, or ""
                the name is <filterName>_camFlux or camFlux if filterName is None
            @param[in] refFilterName  reference filter name; <refFilterName>_flux must exist
            """
            camFluxName = filterName + "_camFlux" if filterName is not None else "camFlux"
            refFluxName = refFilterName + "_flux"
            if refFluxName not in schema:
                raise RuntimeError("Unknown reference filter %s" % (refFluxName,))
            aliasMap.set(camFluxName, refFluxName)
            refFluxErrName = refFluxName + "Err"
            if refFluxErrName in schema:
                camFluxErrName = camFluxName + "Err"
                aliasMap.set(camFluxErrName, refFluxErrName)

        if self.config.defaultFilter:
            addAliasesForOneFilter(None, self.config.defaultFilter)

        for filterName, refFilterName in self.config.filterMap.items():
            addAliasesForOneFilter(filterName, refFilterName)

    @staticmethod
    def makeMinimalSchema(filterNameList, addFluxErr=False,
                          addIsPhotometric=False, addIsResolved=False, addIsVariable=False):
        """!Make the standard schema for reference object catalogs

        @param[in] filterNameList  list of filter names; used to create *filterName*_flux fields
        @param[in] addFluxErr  if True then include flux sigma fields
        @param[in] addIsPhotometric  if True add field "photometric"
        @param[in] addIsResolved  if True add field "resolved"
        @param[in] addIsVariable  if True add field "variable"
        """
        schema = afwTable.SimpleTable.makeMinimalSchema()
        afwTable.Point2DKey.addFields(
            schema,
            "centroid",
            "centroid on an exposure, if relevant",
            "pixel",
        )
        schema.addField(
            field="hasCentroid",
            type="Flag",
            doc="is position known?",
        )
        for filterName in filterNameList:
            schema.addField(
                field="%s_flux" % (filterName,),
                type=numpy.float64,
                doc="flux in filter %s" % (filterName,),
                units="Jy",
            )
        if addFluxErr:
            for filterName in filterNameList:
                schema.addField(
                    field="%s_fluxErr" % (filterName,),
                    type=numpy.float64,
                    doc="flux uncertainty in filter %s" % (filterName,),
                    units="Jy",
                )
        if addIsPhotometric:
            schema.addField(
                field="photometric",
                type="Flag",
                doc="set if the object can be used for photometric calibration",
            )
        if addIsResolved:
            schema.addField(
                field="resolved",
                type="Flag",
                doc="set if the object is spatially resolved",
            )
        if addIsVariable:
            schema.addField(
                field="variable",
                type="Flag",
                doc="set if the object has variable brightness",
            )
        return schema

    def _calculateCircle(self, bbox, wcs):
        """!Compute on-sky center and radius of search region

        @param[in] bbox  bounding box for pixels (an lsst.geom.Box2I or Box2D)
        @param[in] wcs  WCS (an lsst.afw.geom.SkyWcs)
        @return an lsst.pipe.base.Struct containing:
        - coord: ICRS center of the search region (lsst.geom.SpherePoint)
        - radius: the radius of the search region (lsst.geom.Angle)
        - bbox: the bounding box used to compute the circle (lsst.geom.Box2D)
        """
        bbox = lsst.geom.Box2D(bbox)  # make sure bbox is double and that we have a copy
        bbox.grow(self.config.pixelMargin)
        coord = wcs.pixelToSky(bbox.getCenter())
        radius = max(coord.separation(wcs.pixelToSky(pp)) for pp in bbox.getCorners())
        return pipeBase.Struct(coord=coord, radius=radius, bbox=bbox)

    def getMetadataBox(self, bbox, wcs, filterName=None, calib=None):
        """!Return metadata about the load

        This metadata is used for reloading the catalog (e.g., for
        reconstituting a normalised match list.

        @param[in] bbox  bounding box for pixels (an lsst.geom.Box2I or Box2D)
        @param[in] wcs  WCS (an lsst.afw.geom.SkyWcs)
        @param[in] filterName  name of camera filter, or None or blank for the default filter
        @param[in] calib  calibration, or None if unknown
        @return metadata (lsst.daf.base.PropertyList)
        """
        circle = self._calculateCircle(bbox, wcs)
        return self.getMetadataCircle(circle.coord, circle.radius, filterName, calib)

    def getMetadataCircle(self, coord, radius, filterName, calib=None):
        """!Return metadata about the load

        This metadata is used for reloading the catalog (e.g., for
        reconstituting a normalised match list.

        @param[in] coord  ICRS centr of circle (lsst.geom.SpherePoint)
        @param[in] radius  radius of circle (lsst.geom.Angle)
        @param[in] filterName  name of camera filter, or None or blank for the default filter
        @param[in] calib  calibration, or None if unknown
        @return metadata (lsst.daf.base.PropertyList)
        """
        md = PropertyList()
        md.add('RA', coord.getRa().asDegrees(), 'field center in degrees')
        md.add('DEC', coord.getDec().asDegrees(), 'field center in degrees')
        md.add('RADIUS', radius.asDegrees(), 'field radius in degrees, minimum')
        md.add('SMATCHV', 1, 'SourceMatchVector version number')
        filterName = "UNKNOWN" if filterName is None else str(filterName)
        md.add('FILTER', filterName, 'filter name for photometric data')
        return md

    def joinMatchListWithCatalog(self, matchCat, sourceCat):
        """!Relink an unpersisted match list to sources and reference objects

        A match list is persisted and unpersisted as a catalog of IDs produced by
        afw.table.packMatches(), with match metadata (as returned by the astrometry tasks)
        in the catalog's metadata attribute.  This method converts such a match catalog
        into a match list (an lsst.afw.table.ReferenceMatchVector) with links to source
        records and reference object records.

        @param[in]     matchCat   Unperisted packed match list (an lsst.afw.table.BaseCatalog).
                                  matchCat.table.getMetadata() must contain match metadata,
                                  as returned by the astrometry tasks.
        @param[in,out] sourceCat  Source catalog (an lsst.afw.table.SourceCatalog).
                                  As a side effect, the catalog will be sorted by ID.

        @return the match list (an lsst.afw.table.ReferenceMatchVector)
        """
        matchmeta = matchCat.table.getMetadata()
        version = matchmeta.getInt('SMATCHV')
        if version != 1:
            raise ValueError('SourceMatchVector version number is %i, not 1.' % version)
        filterName = matchmeta.getString('FILTER').strip()
        ctrCoord = lsst.geom.SpherePoint(matchmeta.getDouble('RA'),
                                         matchmeta.getDouble('DEC'), lsst.geom.degrees)
        rad = matchmeta.getDouble('RADIUS') * lsst.geom.degrees
        refCat = self.loadSkyCircle(ctrCoord, rad, filterName).refCat
        refCat.sort()
        sourceCat.sort()
        return afwTable.unpackMatches(matchCat, refCat, sourceCat)
