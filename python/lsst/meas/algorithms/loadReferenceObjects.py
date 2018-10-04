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

import astropy.time
import astropy.units
import numpy

import lsst.geom
import lsst.afw.table as afwTable
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pex.exceptions as pexExcept
from lsst.daf.base import PropertyList


def getRefFluxField(schema, filterName=None):
    """Get the name of a flux field from a schema.

    if filterName is specified:
        return *filterName*_camFlux if present
        else return *filterName*_flux if present (camera filter name matches reference filter name)
        else throw RuntimeError
    else:
        return camFlux, if present,
        else throw RuntimeError

    Parameters
    ----------
    schema : `lsst.afw.table.Schema`
        Reference catalog schema.
    filterName : `str`
        Name of camera filter.

    Returns
    -------
    fluxFieldName : `str`
        Name of flux field.

    Raises
    ------
    RuntimeError
        If an appropriate field is not found.
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
    """Return keys for flux and flux error.

    Parameters
    ----------
    schema : `lsst.afw.table.Schema`
        Reference catalog schema.
    filterName : `str`
        Name of camera filter.

    Returns
    -------
    keys : `tuple` of (`lsst.afw.table.Key`, `lsst.afw.table.Key`)
        Two keys:

        - flux key
        - flux error key, if present, else None

    Raises
    ------
    RuntimeError
        If flux field not found.
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
    requireProperMotion = pexConfig.Field(
        doc="Require that the fields needed to correct proper motion "
            "(epoch, pm_ra and pm_dec) are present?",
        dtype=bool,
        default=False,
    )

# The following comment block adds a link to this task from the Task Documentation page.
## @addtogroup LSST_task_documentation
## @{
## @page LoadReferenceObjectsTask
## @ref LoadReferenceObjectsTask_ "LoadReferenceObjectsTask"
## @copybrief LoadReferenceObjectsTask
## @}


class LoadReferenceObjectsTask(pipeBase.Task, metaclass=abc.ABCMeta):
    # Parameters
    # ----------
    # butler : `lsst.daf.persistence.Butler`
    #     Data butler, for access reference catalogs.
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
    (other fields may also be present).
    The units use astropy quantity conventions, so a 2 suffix means squared.
    See also makeMinimalSchema.
    - coord: ICRS position of star on sky (an lsst.geom.SpherePoint)
    - centroid: position of star on an exposure, if relevant (an lsst.afw.Point2D)
    - hasCentroid: is centroid usable? (a Flag)
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
    - photometric (optional): is the object usable for photometric calibration? (a Flag)
    - resolved (optional): is the object spatially resolved? (a Flag)
    - variable (optional): does the object have variable brightness? (a Flag)
    - coord_raErr: uncertainty in `coord` along the direction of right ascension (radian, an Angle)
                    = uncertainty in ra * cos(dec); nan if unknown
    - coord_decErr: uncertainty in `coord` along the direction of declination (radian, an Angle);
        nan if unknown
    The following are optional; fields should only be present if the
    information is available for at least some objects.
    Numeric values are `nan` if unknown:
    - epoch: date of observation as TAI MJD (day)
    - pm_ra: proper motion along the direction of right ascension (rad/year, an Angle) = dra/dt * cos(dec)
    - pm_dec: proper motion along the direction of declination (rad/year, and Angle)
    - pm_raErr: uncertainty in `pm_ra` (rad/year)
    - pm_decErr: uncertainty in `pm_dec` (rad/year)
    - pm_ra_dec_Cov: covariance between pm_ra and pm_dec (rad2/year2)
    - pm_flag: set if proper motion, error or covariance is bad
    - parallax: parallax (rad, an Angle)
    - parallaxErr: uncertainty in `parallax` (rad)
    - parallax_flag: set if parallax value or parallaxErr is bad
    - coord_ra_pm_ra_Cov: covariance between coord_ra and pm_ra (rad2/year)
    - coord_ra_pm_dec_Cov: covariance between coord_ra and pm_dec (rad2/year)
    - coord_ra_parallax_Cov: covariance between coord_ra and parallax (rad2/year)
    - coord_dec_pm_ra_Cov: covariance between coord_dec and pm_ra (rad2/year)
    - coord_dec_pm_dec_Cov: covariance between coord_dec and pm_dec (rad2/year)
    - coord_dec_parallax_Cov: covariance between coord_dec and parallax (rad2/year)
    - pm_ra_parallax_Cov: covariance between pm_ra and parallax (rad2/year)
    - pm_dec_parallax_Cov: covariance between pm_dec and parallax (rad2/year)
    @section meas_algorithms_loadReferenceObjects_Config       Configuration parameters
    See @ref LoadReferenceObjectsConfig for a base set of configuration parameters.
    Most subclasses will add configuration variables.
    """
    ConfigClass = LoadReferenceObjectsConfig
    _DefaultName = "LoadReferenceObjects"

    def __init__(self, butler=None, *args, **kwargs):
        pipeBase.Task.__init__(self, *args, **kwargs)
        self.butler = butler

    @pipeBase.timeMethod
    def loadPixelBox(self, bbox, wcs, filterName=None, calib=None, epoch=None):
        """Load reference objects that overlap a rectangular pixel region.

        Parameters
        ----------
        bbox : `lsst.geom.Box2I` or `lsst.geom.Box2D`
            Bounding box for pixels.
        wcs : `lsst.afw.geom.SkyWcs`
            WCS; used to convert pixel positions to sky coordinates
            and vice-versa.
        filterName : `str`
            Name of filter, or `None` or `""` for the default filter.
            This is used for flux values in case we have flux limits
            (which are not yet implemented).
        calib : `lsst.afw.image.Calib` (optional)
            Calibration, or `None` if unknown.
        epoch : `astropy.time.Time` (optional)
            Epoch to which to correct proper motion and parallax,
            or None to not apply such corrections.

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            A Struct containing the following fields:
            refCat : `lsst.afw.catalog.SimpleCatalog`
                A catalog of reference objects with the standard
                schema, as documented in the main doc string for
                `LoadReferenceObjects`.
                The catalog is guaranteed to be contiguous.
            fluxField : `str`
                Name of flux field for specified `filterName`.

        Notes
        -----
        The search algorithm works by searching in a region in sky
        coordinates whose center is the center of the bbox and radius
        is large enough to just include all 4 corners of the bbox.
        Stars that lie outside the bbox are then trimmed from the list.
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
    def loadSkyCircle(self, ctrCoord, radius, filterName=None, epoch=None):
        """Load reference objects that overlap a circular sky region.

        Parameters
        ----------
        ctrCoord : `lsst.geom.SpherePoint`
            ICRS center of search region.
        radius : `lsst.geom.Angle`
            Radius of search region.
        filterName : `str` (optional)
            Name of filter, or `None` or `""` for the default filter.
            This is used for flux values in case we have flux limits
            (which are not yet implemented).
        epoch : `astropy.time.Time` (optional)
            Epoch to which to correct proper motion and parallax,
            or None to not apply such corrections.

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            A Struct containing the following fields:
            refCat : `lsst.afw.catalog.SimpleCatalog`
                A catalog of reference objects with the standard
                schema, as documented in the main doc string for
                `LoadReferenceObjects`.
                The catalog is guaranteed to be contiguous.
            fluxField : `str`
                Name of flux field for specified `filterName`.

        Notes
        -----
        Note that subclasses are responsible for performing the proper motion
        correction, since this is the lowest-level interface for retrieving
        the catalog.
        """
        return

    @staticmethod
    def _trimToBBox(refCat, bbox, wcs):
        """Remove objects outside a given pixel bounding box and set
        centroid and hasCentroid fields.

        Parameters
        ----------
        refCat : `lsst.afw.table.SimpleCatalog`
            A catalog of objects. The schema must include fields
            "coord", "centroid" and "hasCentroid".
            The "coord" field is read.
            The "centroid" and "hasCentroid" fields are set.
        bbox : `lsst.geom.Box2D`
            Pixel region
        wcs : `lsst.afw.geom.SkyWcs`
            WCS; used to convert sky coordinates to pixel positions.

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
        """Add aliases for camera filter fluxes to the schema.

        If self.config.defaultFilter then adds these aliases:
            camFlux:      <defaultFilter>_flux
            camFluxErr: <defaultFilter>_fluxErr, if the latter exists

        For each camFilter: refFilter in self.config.filterMap adds these aliases:
            <camFilter>_camFlux:      <refFilter>_flux
            <camFilter>_camFluxErr: <refFilter>_fluxErr, if the latter exists

        Parameters
        ----------
        schema : `lsst.afw.table.Schema`
            Schema for reference catalog.

        Throws
        ------
        RuntimeError
            If any reference flux field is missing from the schema.
        """
        aliasMap = schema.getAliasMap()

        def addAliasesForOneFilter(filterName, refFilterName):
            """Add aliases for a single filter

            Parameters
            ----------
            filterName : `str` (optional)
                Camera filter name. The resulting alias name is
                <filterName>_camFlux, or simply "camFlux" if `filterName`
                is `None` or `""`.
            refFilterName : `str`
                Reference catalog filter name; the field
                <refFilterName>_flux must exist.
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
    def makeMinimalSchema(filterNameList, *, addFluxErr=False, addCentroid=True,
                          addIsPhotometric=False, addIsResolved=False,
                          addIsVariable=False, coordErrDim=2,
                          addProperMotion=False, properMotionErrDim=2,
                          addParallax=False, addParallaxErr=True):
        """Make a standard schema for reference object catalogs.

        Parameters
        ----------
        filterNameList : `list` of `str`
            List of filter names. Used to create <filterName>_flux fields.
        addFluxErr : `bool`
            If True then include flux sigma fields.
        addIsPhotometric : `bool`
            If True then add field "photometric".
        addIsResolved : `bool`
            If True then add field "resolved".
        addIsVariable : `bool`
            If True then add field "variable".
        coordErrDim : `int`
            Number of coord error fields; must be one of 0, 2, 3:

            - If 2 or 3: add fields "coord_raErr" and "coord_decErr".
            - If 3: also add field "coord_radecErr".
        addProperMotion : `bool`
            If True add fields "epoch", "pm_ra", "pm_dec" and "pm_flag".
        properMotionErrDim : `int`
            Number of proper motion error fields; must be one of 0, 2, 3;
            ignored if addProperMotion false:
            - If 2 or 3: add fields "pm_raErr" and "pm_decErr".
            - If 3: also add field "pm_radecErr".
        addParallax : `bool`
            If True add fields "epoch", "parallax", "parallaxErr"
            and "parallax_flag".
        addParallaxErr : `bool`
            If True add field "parallaxErr"; ignored if addParallax false.

        Returns
        -------
        schema : `lsst.afw.table.Schema`
            Schema for reference catalog, an
            `lsst.afw.table.SimpleCatalog`.

        Notes
        -----
        Reference catalogs support additional covariances, such as
        covariance between RA and proper motion in declination,
        that are not supported by this method, but can be added after
        calling this method.
        """
        schema = afwTable.SimpleTable.makeMinimalSchema()
        if addCentroid:
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
        if coordErrDim not in (0, 2, 3):
            raise ValueError("coordErrDim={}; must be (0, 2, 3)".format(coordErrDim))
        if coordErrDim > 0:
            afwTable.CovarianceMatrix2fKey.addFields(
                schema=schema,
                prefix="coord",
                names=["ra", "dec"],
                units=["rad", "rad"],
                diagonalOnly=(coordErrDim == 2),
            )

        if addProperMotion or addParallax:
            schema.addField(
                field="epoch",
                type=numpy.float64,
                doc="date of observation (TAI, MJD)",
                units="day",
            )

        if addProperMotion:
            schema.addField(
                field="pm_ra",
                type="Angle",
                doc="proper motion in the right ascension direction = dra/dt * cos(dec)",
                units="rad/year",
            )
            schema.addField(
                field="pm_dec",
                type="Angle",
                doc="proper motion in the declination direction",
                units="rad/year",
            )
            if properMotionErrDim not in (0, 2, 3):
                raise ValueError("properMotionErrDim={}; must be (0, 2, 3)".format(properMotionErrDim))
            if properMotionErrDim > 0:
                afwTable.CovarianceMatrix2fKey.addFields(
                    schema=schema,
                    prefix="pm",
                    names=["ra", "dec"],
                    units=["rad/year", "rad/year"],
                    diagonalOnly=(properMotionErrDim == 2),
                )
            schema.addField(
                field="pm_flag",
                type="Flag",
                doc="Set if proper motion or proper motion error is bad",
            )

        if addParallax:
            schema.addField(
                field="parallax",
                type=numpy.float64,
                doc="parallax",
                units="rad",
            )
            if addParallaxErr:
                schema.addField(
                    field="parallaxErr",
                    type=numpy.float64,
                    doc="uncertainty in parallax",
                    units="rad",
                )
            schema.addField(
                field="parallax_flag",
                type="Flag",
                doc="Set if parallax or parallax error is bad",
            )
        return schema

    def _calculateCircle(self, bbox, wcs):
        """Compute on-sky center and radius of search region.

        Parameters
        ----------
        bbox : `lsst.geom.Box2I` or `lsst.geom.Box2D`
            Pixel bounding box.
        wcs : `lsst.afw.geom.SkyWcs`
            WCS; used to convert pixel positions to sky coordinates.

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            A Struct containing:

            - coord : `lsst.geom.SpherePoint`
                ICRS center of the search region.
            - radius : `lsst.geom.Angle`
                Radius of the search region.
            - bbox : `lsst.afw.geom.Box2D`
                Bounding box used to compute the circle.
        """
        bbox = lsst.geom.Box2D(bbox)  # make sure bbox is double and that we have a copy
        bbox.grow(self.config.pixelMargin)
        coord = wcs.pixelToSky(bbox.getCenter())
        radius = max(coord.separation(wcs.pixelToSky(pp)) for pp in bbox.getCorners())
        return pipeBase.Struct(coord=coord, radius=radius, bbox=bbox)

    def getMetadataBox(self, bbox, wcs, filterName=None, calib=None, epoch=None):
        """Return metadata about the load.

        This metadata is used for reloading the catalog (e.g., for
        reconstituting a normalised match list.

        Parameters
        ----------
        bbox : `lsst.geom.Box2I` or `lsst.geom.Box2D`
            Pixel bounding box.
        wcs : `lsst.afw.geom.SkyWcs`
            WCS; used to convert pixel positions to sky coordinates.
        filterName : `str`
            Name of camera filter, or `None` or `""` for the default
            filter.
        calib : `lsst.afw.image.Calib` (optional)
            Calibration, or `None` if unknown.
        epoch : `astropy.time.Time` (optional)
            Epoch to which to correct proper motion and parallax,
            or None to not apply such corrections.

        Returns
        -------
        metadata : lsst.daf.base.PropertyList
            Metadata about the load.
        """
        circle = self._calculateCircle(bbox, wcs)
        return self.getMetadataCircle(circle.coord, circle.radius, filterName, calib)

    def getMetadataCircle(self, coord, radius, filterName, calib=None, epoch=None):
        """Return metadata about the load.

        This metadata is used for reloading the catalog (e.g., for
        reconstituting a normalised match list.

        Parameters
        ----------
        coord : `lsst.geom.SpherePoint`
            ICRS center of the search region.
        radius : `lsst.geom.Angle`
            Radius of the search region.
        filterName : `str`
            Name of camera filter, or `None` or `""` for the default
            filter.
        calib : `lsst.afw.image.Calib` (optional)
            Calibration, or `None` if unknown.
        epoch : `astropy.time.Time` (optional)
            Epoch to which to correct proper motion and parallax,
            or None to not apply such corrections.

        Returns
        -------
        metadata : lsst.daf.base.PropertyList
            Metadata about the load
        """
        md = PropertyList()
        md.add('RA', coord.getRa().asDegrees(), 'field center in degrees')
        md.add('DEC', coord.getDec().asDegrees(), 'field center in degrees')
        md.add('RADIUS', radius.asDegrees(), 'field radius in degrees, minimum')
        md.add('SMATCHV', 1, 'SourceMatchVector version number')
        filterName = "UNKNOWN" if filterName is None else str(filterName)
        md.add('FILTER', filterName, 'filter name for photometric data')
        md.add('EPOCH', "NONE" if epoch is None else epoch, 'Epoch (TAI MJD) for catalog')
        return md

    def joinMatchListWithCatalog(self, matchCat, sourceCat):
        """Relink an unpersisted match list to sources and reference
        objects.

        A match list is persisted and unpersisted as a catalog of IDs
        produced by afw.table.packMatches(), with match metadata
        (as returned by the astrometry tasks) in the catalog's metadata
        attribute. This method converts such a match catalog into a match
        list, with links to source records and reference object records.

        Parameters
        ----------
        matchCat : `lsst.afw.table.BaseCatalog`
            Unperisted packed match list.
            ``matchCat.table.getMetadata()`` must contain match metadata,
            as returned by the astrometry tasks.
        sourceCat : `lsst.afw.table.SourceCatalog`
            Source catalog. As a side effect, the catalog will be sorted
            by ID.

        Returns
        -------
        matchList : `lsst.afw.table.ReferenceMatchVector`
            Match list.
        """
        matchmeta = matchCat.table.getMetadata()
        version = matchmeta.getInt('SMATCHV')
        if version != 1:
            raise ValueError('SourceMatchVector version number is %i, not 1.' % version)
        filterName = matchmeta.getString('FILTER').strip()
        ctrCoord = lsst.geom.SpherePoint(matchmeta.getDouble('RA'),
                                         matchmeta.getDouble('DEC'), lsst.geom.degrees)
        rad = matchmeta.getDouble('RADIUS') * lsst.geom.degrees
        try:
            epoch = matchmeta.getDouble('EPOCH')
        except (pexExcept.NotFoundError, pexExcept.TypeError):
            epoch = None  # Not present, or not correct type means it's not set
        refCat = self.loadSkyCircle(ctrCoord, rad, filterName, epoch=epoch).refCat
        refCat.sort()
        sourceCat.sort()
        return afwTable.unpackMatches(matchCat, refCat, sourceCat)

    def applyProperMotions(self, catalog, epoch):
        """Apply proper motion correction to a reference catalog.

        Adjust position and position error in the ``catalog``
        for proper motion to the specified ``epoch``,
        modifying the catalong in place.

        Parameters
        ----------
        catalog : `lsst.afw.table.SimpleCatalog`
            Catalog of positions, containing:

            - Coordinates, retrieved by the table's coordinate key.
            - ``coord_raErr`` : Error in Right Ascension (rad).
            - ``coord_decErr`` : Error in Declination (rad).
            - ``pm_ra`` : Proper motion in Right Ascension (rad/yr,
                East positive)
            - ``pm_raErr`` : Error in ``pm_ra`` (rad/yr), optional.
            - ``pm_dec`` : Proper motion in Declination (rad/yr,
                North positive)
            - ``pm_decErr`` : Error in ``pm_dec`` (rad/yr), optional.
            - ``epoch`` : Mean epoch of object (an astropy.time.Time)
        epoch : `astropy.time.Time` (optional)
            Epoch to which to correct proper motion and parallax,
            or None to not apply such corrections.
        """
        if ("epoch" not in catalog.schema or "pm_ra" not in catalog.schema or "pm_dec" not in catalog.schema):
            if self.config.requireProperMotion:
                raise RuntimeError("Proper motion correction required but not available from catalog")
            self.log.warn("Proper motion correction not available from catalog")
            return
        if not catalog.isContiguous():
            raise RuntimeError("Catalog must be contiguous")
        catEpoch = astropy.time.Time(catalog["epoch"], scale="tai", format="mjd")
        self.log.debug("Correcting reference catalog for proper motion to %r", epoch)
        # Use `epoch.tai` to make sure the time difference is in TAI
        timeDiffsYears = (epoch.tai - catEpoch).to(astropy.units.yr).value
        coordKey = catalog.table.getCoordKey()
        # Compute the offset of each object due to proper motion
        # as components of the arc of a great circle along RA and Dec
        pmRaRad = catalog["pm_ra"]
        pmDecRad = catalog["pm_dec"]
        offsetsRaRad = pmRaRad*timeDiffsYears
        offsetsDecRad = pmDecRad*timeDiffsYears
        # Compute the corresponding bearing and arc length of each offset
        # due to proper motion, and apply the offset
        # The factor of 1e6 for computing bearing is intended as
        # a reasonable scale for typical values of proper motion
        # in order to avoid large errors for small values of proper motion;
        # using the offsets is another option, but it can give
        # needlessly large errors for short duration
        offsetBearingsRad = numpy.arctan2(pmDecRad*1e6, pmRaRad*1e6)
        offsetAmountsRad = numpy.hypot(offsetsRaRad, offsetsDecRad)
        for record, bearingRad, amountRad in zip(catalog, offsetBearingsRad, offsetAmountsRad):
            record.set(coordKey,
                       record.get(coordKey).offset(bearing=bearingRad*lsst.geom.radians,
                                                   amount=amountRad*lsst.geom.radians))
        # Increase error in RA and Dec based on error in proper motion
        if "coord_raErr" in catalog.schema:
            catalog["coord_raErr"] = numpy.hypot(catalog["coord_raErr"],
                                                 catalog["pm_raErr"]*timeDiffsYears)
        if "coord_decErr" in catalog.schema:
            catalog["coord_decErr"] = numpy.hypot(catalog["coord_decErr"],
                                                  catalog["pm_decErr"]*timeDiffsYears)
