from __future__ import absolute_import, division, print_function

import abc

import numpy

import lsst.afw.table as afwTable
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

__all__ = ["getRefFluxField", "getRefFluxKeys", "LoadReferenceObjectsTask", "LoadReferenceObjectsConfig"]

def getRefFluxField(schema, filterName=None):
    """!Get name of flux field in schema

    if filterName is specified:
        return <filterName>_camFlux if present
        else return <filterName>_flux if present (camera filter name matches reference filter name)
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
    fluxErrField = fluxField + "Sigma"
    fluxKey = schema[fluxField].asKey()
    try:
        fluxErrKey = schema[fluxErrField].asKey()
    except Exception:
        fluxErrKey = None
    return (fluxKey, fluxErrKey)

class LoadReferenceObjectsConfig(pexConfig.Config):
    pixelMargin = pexConfig.RangeField(
        doc = "Padding to add to 4 all edges of the bounding box (pixels)",
        dtype = int,
        default = 50,
        min = 0,
    )
    defaultFilter = pexConfig.Field(
        doc = "Default reference catalog filter to use if filter not specified in exposure; " + \
            "if blank then filter must be specified in exposure",
        dtype = str,
        default = "",
    )
    filterMap = pexConfig.DictField(
        doc = "Mapping of camera filter name: reference catalog filter name; " + \
            "each reference filter must exist",
        keytype = str,
        itemtype = str,
        default = {},
    )

class LoadReferenceObjectsTask(pipeBase.Task):
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

    Implementations must subclass this class, override the loadObjectsInBBox method,
    and will typically override the value of ConfigClass with a task-specific config class.

    @section meas_algorithms_loadReferenceObjects_Initialize   Task initialisation

    @copydoc \_\_init\_\_

    @section meas_algorithms_loadReferenceObjects_IO       Invoking the Task

    @copydoc loadObjectsInBBox

    @section meas_algorithms_loadReferenceObjects_Schema       Schema of the reference object catalog

    Reference object catalogs are instances of lsst.afw.table.SimpleCatalog with the following schema
    (other fields may also be present):
    - coord: position of star on sky (an lsst.afw.coord.IcrsCoord)
    - centroid: position of star on an exposure, if relevant (an lsst.afw.Point2D)
    - hasCentroid: is centroid usable?
    - *referenceFilterName*_flux: brightness in the specified reference catalog filter: 10^(-0.4*mag)
    - *referenceFilterName*_fluxSigma (optional): brightness standard deviation;
        omitted if no data is available; possibly nan if data is available for some objects but not others
    - *cameraFilterName*_camera_flux: brightness in specified camera filter (magnitude)
    - *cameraFilterName*_camera_fluxSigma (optional): brightness standard deviation
        in specified camera filter (magnitude); omitted if no data is available;
        possibly nan if data is available for some objects but not others
    - default_flux (optional): brightness to use if no camera filter is available (magnitude);
        omitted unless defaultFilter is specified in the config
    - default_fluxSigma (optional): brightness standard deviation to use if no camera filter is available
        (magnitude); omitted unless defaultFilter is specified in the config and the corresponding
        fluxSigma field exists
    - photometric (optional): is the object usable for photometric calibration?
    - resolved (optional): is the object spatially resolved?
    - variable (optional): does the object have variable brightness?

    @section meas_algorithms_loadReferenceObjects_Config       Configuration parameters

    See @ref LoadReferenceObjectsConfig for a base set of configuration parameters.
    Most subclasses will add configuration variables.
    """
    __metaclass__ = abc.ABCMeta
    ConfigClass = LoadReferenceObjectsConfig
    _DefaultName = "LoadReferenceObjects"

    @abc.abstractmethod
    def loadObjectsInBBox(self, bbox, wcs, filterName=None, calib=None):
        """!Load reference objects that overlap a pixel-based rectangular region

        The search algorith works by searching in a region in sky coordinates whose center is the center
        of the bbox and radius is large enough to just include all 4 corners of the bbox.
        Stars that lie outside the bbox are then trimmed from the list.

        @param[in] bbox  bounding box for pixels (an lsst.afw.geom.Box2I or Box2D)
        @param[in] wcs  WCS (an lsst.afw.image.Wcs)
        @param[in] filterName  name of camera filter, or None or blank for the default filter
        @param[in] calib  calibration, or None if unknown

        @return a catalog of reference objects using the standard schema (see the class doc string)
        """
        return

    def _addFluxAliases(self, schema):
        """Add aliases for camera filter fluxes to the schema

        If self.config.defaultFilter then adds these aliases:
            camFlux:      <defaultFilter>_flux
            camFluxSigma: <defaultFilter>_fluxSigma, if the latter exists

        For each camFilter: refFilter in self.config.filterMap adds these aliases:
            <camFilter>_camFlux:      <refFilter>_flux
            <camFilter>_camFluxSigma: <refFilter>_fluxSigma, if the latter exists

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
            refFluxErrName = refFluxName + "Sigma"
            if refFluxErrName in schema:
                camFluxErrName = camFluxName + "Sigma"
                aliasMap.set(camFluxErrName, refFluxErrName)

        if self.config.defaultFilter:
            addAliasesForOneFilter(None, self.config.defaultFilter)

        for filterName, refFilterName in self.config.filterMap.iteritems():
            addAliasesForOneFilter(filterName, refFilterName)

    @staticmethod
    def makeMinimalSchema(filterNameList, addFluxSigma=False,
            addIsPhotometric=False, addIsResolved=False, addIsVariable=False):
        """!Make the standard schema for reference object catalogs

        @param[in] filterNameList  list of filter names; used to create <filterName>_flux fields
        @param[in] addFluxSigma  if True then include flux sigma fields
        @param[in] addIsPhotometric  if True add field "photometric"
        @param[in] addIsResolved  if True add field "resolved"
        """
        schema = afwTable.SimpleTable.makeMinimalSchema()
        schema.addField(
            field = "centroid",
            type = "PointD",
            doc = "centroid on an exposure, if relevant",
            units = "pixels",
        )
        schema.addField(
            field = "hasCentroid",
            type = "Flag",
            doc = "is position known?",
        )
        for filterName in filterNameList:
            schema.addField(
                field = "%s_flux" % (filterName,),
                type = numpy.float64,
                doc = "flux in filter %s" % (filterName,),
                units = "?",
            )
        if addFluxSigma:
            for filterName in filterNameList:
                schema.addField(
                    field = "%s_fluxSigma" % (filterName,),
                    type = numpy.float64,
                    doc = "flux uncertainty in filter %s" % (filterName,),
                    units = "?",
                )
        if addIsPhotometric:
            schema.addField(
                field = "photometric",
                type = "Flag",
                doc = "set if the object can be used for photometric calibration",
            )
        if addIsResolved:
            schema.addField(
                field = "resolved",
                type = "Flag",
                doc = "set if the object is spatially resolved",
            )
        if addIsVariable:
            schema.addField(
                field = "variable",
                type = "Flag",
                doc = "set if the object has variable brightness",
            )
        return schema
