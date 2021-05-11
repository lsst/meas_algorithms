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

__all__ = ["getRefFluxField", "getRefFluxKeys", "LoadReferenceObjectsTask", "LoadReferenceObjectsConfig",
           "ReferenceObjectLoader"]

import abc
import itertools

import astropy.time
import astropy.units
import numpy

import lsst.geom as geom
import lsst.afw.table as afwTable
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.log
from lsst import sphgeom
from lsst.daf.base import PropertyList


def isOldFluxField(name, units):
    """Return True if this name/units combination corresponds to an
    "old-style" reference catalog flux field.
    """
    unitsCheck = units != 'nJy'  # (units == 'Jy' or units == '' or units == '?')
    isFlux = name.endswith('_flux')
    isFluxSigma = name.endswith('_fluxSigma')
    isFluxErr = name.endswith('_fluxErr')
    return (isFlux or isFluxSigma or isFluxErr) and unitsCheck


def hasNanojanskyFluxUnits(schema):
    """Return True if the units of all flux and fluxErr are correct (nJy).
    """
    for field in schema:
        if isOldFluxField(field.field.getName(), field.field.getUnits()):
            return False
    return True


def getFormatVersionFromRefCat(refCat):
    """"Return the format version stored in a reference catalog header.

    Parameters
    ----------
    refCat : `lsst.afw.table.SimpleCatalog`
        Reference catalog to inspect.

    Returns
    -------
    version : `int` or `None`
        Format version integer, or `None` if the catalog has no metadata
        or the metadata does not include a "REFCAT_FORMAT_VERSION" key.
    """
    md = refCat.getMetadata()
    if md is None:
        return None
    try:
        return md.getScalar("REFCAT_FORMAT_VERSION")
    except KeyError:
        return None


def convertToNanojansky(catalog, log, doConvert=True):
    """Convert fluxes in a catalog from jansky to nanojansky.

    Parameters
    ----------
    catalog : `lsst.afw.table.SimpleCatalog`
        The catalog to convert.
    log : `lsst.log.Log`
        Log to send messages to.
    doConvert : `bool`, optional
        Return a converted catalog, or just identify the fields that need to be converted?
        This supports the "write=False" mode of `bin/convert_to_nJy.py`.

    Returns
    -------
    catalog : `lsst.afw.table.SimpleCatalog` or None
        The converted catalog, or None if ``doConvert`` is False.

    Notes
    -----
    Support for old units in reference catalogs will be removed after the
    release of late calendar year 2019.
    Use `meas_algorithms/bin/convert_to_nJy.py` to update your reference catalog.
    """
    # Do not share the AliasMap: for refcats, that gets created when the
    # catalog is read from disk and should not be propagated.
    mapper = lsst.afw.table.SchemaMapper(catalog.schema, shareAliasMap=False)
    mapper.addMinimalSchema(lsst.afw.table.SimpleTable.makeMinimalSchema())
    input_fields = []
    output_fields = []
    for field in catalog.schema:
        oldName = field.field.getName()
        oldUnits = field.field.getUnits()
        if isOldFluxField(oldName, oldUnits):
            units = 'nJy'
            # remap Sigma flux fields to Err, so we can drop the alias
            if oldName.endswith('_fluxSigma'):
                name = oldName.replace('_fluxSigma', '_fluxErr')
            else:
                name = oldName
            newField = lsst.afw.table.Field[field.dtype](name, field.field.getDoc(), units)
            mapper.addMapping(field.getKey(), newField)
            input_fields.append(field.field)
            output_fields.append(newField)
        else:
            mapper.addMapping(field.getKey())

    fluxFieldsStr = '; '.join("(%s, '%s')" % (field.getName(), field.getUnits()) for field in input_fields)

    if doConvert:
        newSchema = mapper.getOutputSchema()
        output = lsst.afw.table.SimpleCatalog(newSchema)
        output.extend(catalog, mapper=mapper)
        for field in output_fields:
            output[field.getName()] *= 1e9
        log.info(f"Converted refcat flux fields to nJy (name, units): {fluxFieldsStr}")
        return output
    else:
        log.info(f"Found old-style refcat flux fields (name, units): {fluxFieldsStr}")
        return None


class _FilterCatalog:
    """This is a private helper class which filters catalogs by
    row based on the row being inside the region used to initialize
    the class.

    Parameters
    ----------
    region : `lsst.sphgeom.Region`
        The spatial region which all objects should lie within
    """
    def __init__(self, region):
        self.region = region

    def __call__(self, refCat, catRegion):
        """This call method on an instance of this class takes in a reference
        catalog, and the region from which the catalog was generated.

        If the catalog region is entirely contained within the region used to
        initialize this class, then all the entries in the catalog must be
        within the region and so the whole catalog is returned.

        If the catalog region is not entirely contained, then the location for
        each record is tested against the region used to initialize the class.
        Records which fall inside this region are added to a new catalog, and
        this catalog is then returned.

        Parameters
        ---------
        refCat : `lsst.afw.table.SourceCatalog`
            SourceCatalog to be filtered.
        catRegion : `lsst.sphgeom.Region`
            Region in which the catalog was created
        """
        if catRegion.isWithin(self.region):
            # no filtering needed, region completely contains refcat
            return refCat

        filteredRefCat = type(refCat)(refCat.table)
        for record in refCat:
            if self.region.contains(record.getCoord().getVector()):
                filteredRefCat.append(record)
        return filteredRefCat


class ReferenceObjectLoaderBase:
    """Base class for reference object loaders, to facilitate gen2/gen3 code
    sharing.
    """
    def applyProperMotions(self, catalog, epoch):
        """Apply proper motion correction to a reference catalog.

        Adjust position and position error in the ``catalog``
        for proper motion to the specified ``epoch``,
        modifying the catalog in place.

        Parameters
        ----------
        catalog : `lsst.afw.table.SimpleCatalog`
            Catalog of positions, containing at least these fields:

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
        epoch : `astropy.time.Time`
            Epoch to which to correct proper motion.
            If None, do not apply PM corrections or raise if
            ``config.requireProperMotion`` is True.

        Raises
        ------
        RuntimeError
            Raised if ``config.requireProperMotion`` is set but we cannot
            apply the proper motion correction for some reason.
        """
        if epoch is None:
            if self.config.requireProperMotion:
                raise RuntimeError("requireProperMotion=True but epoch not provided to loader.")
            else:
                self.log.debug("No epoch provided: not applying proper motion corrections to refcat.")
                return

        # Warn/raise for a catalog in an incorrect format, if epoch was specified.
        if ("pm_ra" in catalog.schema
                and not isinstance(catalog.schema["pm_ra"].asKey(), lsst.afw.table.KeyAngle)):
            if self.config.requireProperMotion:
                raise RuntimeError("requireProperMotion=True but refcat pm_ra field is not an Angle.")
            else:
                self.log.warn("Reference catalog pm_ra field is not an Angle; cannot apply proper motion.")
                return

        if ("epoch" not in catalog.schema or "pm_ra" not in catalog.schema):
            if self.config.requireProperMotion:
                raise RuntimeError("requireProperMotion=True but PM data not available from catalog.")
            else:
                self.log.warn("Proper motion correction not available for this reference catalog.")
            return

        applyProperMotionsImpl(self.log, catalog, epoch)


class ReferenceObjectLoader(ReferenceObjectLoaderBase):
    """This class facilitates loading reference catalogs with gen 3 middleware

    The middleware preflight solver will create a list of datarefs that may
    possibly overlap a given region. These datarefs are then used to construct
    and instance of this class. The class instance should then be passed into
    a task which needs reference catalogs. These tasks should then determine
    the exact region of the sky reference catalogs will be loaded for, and
    call a corresponding method to load the reference objects.
    """
    def __init__(self, dataIds, refCats, config, log=None):
        """ Constructs an instance of ReferenceObjectLoader

        Parameters
        ----------
        dataIds : iterable of `lsst.daf.butler.DataIds`
            An iterable object of DataSetRefs which point to reference catalogs
            in a gen 3 repository.
        refCats : iterable of `lsst.daf.butler.DeferedDatasetHandle`
            Handles to load refCats on demand
        config : `lsst.pex.config.configurableField`
            Configuration for the loader.
        log : `lsst.log.Log` or `None`, optional
            Logger object used to write out messages. If `None` the default
            lsst logger will be used.
        """
        self.dataIds = dataIds
        self.refCats = refCats
        self.log = log or lsst.log.Log.getDefaultLogger()
        self.config = config

    @staticmethod
    def _makeBoxRegion(BBox, wcs, BBoxPadding):
        outerLocalBBox = geom.Box2D(BBox)
        innerLocalBBox = geom.Box2D(BBox)

        # Grow the bounding box to allow for effects not fully captured by the
        # wcs provided (which represents the current best-guess wcs solution
        # associated with the dataset for which the calibration is to be
        # computed using the loaded and trimmed reference catalog being defined
        # here).  These effects could include pointing errors and/or an
        # insufficient optical distorition model for the instrument.  The idea
        # is to ensure the spherical geometric region created contains the
        # entire region covered by the bbox.
        # Also create an inner region that is sure to be inside the bbox.
        outerLocalBBox.grow(BBoxPadding)
        innerLocalBBox.grow(-1*BBoxPadding)

        # Handle the case where the inner bounding box shrank to a zero sized
        # region (which will be the case if the shrunken size of either
        # dimension is less than or equal to zero).  In this case, the inner
        # bounding box is set to the original input bounding box.  This is
        # probably not the best way to handle an empty inner bounding box, but
        # it is what the calling code currently expects.
        if innerLocalBBox.getDimensions() == geom.Extent2D(0, 0):
            innerLocalBBox = geom.Box2D(BBox)

        # Convert the corners of the bounding boxes to sky coordinates.
        innerBoxCorners = innerLocalBBox.getCorners()
        innerSphCorners = [wcs.pixelToSky(corner).getVector() for corner in innerBoxCorners]
        innerSkyRegion = sphgeom.ConvexPolygon(innerSphCorners)

        outerBoxCorners = outerLocalBBox.getCorners()
        outerSphCorners = [wcs.pixelToSky(corner).getVector() for corner in outerBoxCorners]
        outerSkyRegion = sphgeom.ConvexPolygon(outerSphCorners)

        return innerSkyRegion, outerSkyRegion, innerSphCorners, outerSphCorners

    def loadPixelBox(self, bbox, wcs, filterName=None, epoch=None, photoCalib=None,
                     bboxToSpherePadding=100):
        """Load reference objects that are within a pixel-based rectangular
        region.

        This algorithm works by creating a spherical box whose corners
        correspond to the WCS converted corners of the input bounding box
        (possibly padded).  It then defines a filtering function which looks at
        the pixel position of the reference objects and accepts only those that
        lie within the specified bounding box.

        The spherical box region and filtering function are passed to the
        generic loadRegion method which loads and filters the reference objects
        from the datastore and returns a single catalog containing the filtered
        set of reference objects.

        Parameters
        ----------
        bbox : `lsst.geom.Box2I` or `lsst.geom.Box2D`
            Box which bounds a region in pixel space.
        wcs : `lsst.afw.geom.SkyWcs`
            Wcs object defining the pixel to sky (and inverse) transform for
            the supplied ``bbox``.
        filterName : `str` or `None`, optional
            Name of camera filter, or `None` or blank for the default filter.
        epoch : `astropy.time.Time` or `None`, optional
            Epoch to which to correct proper motion and parallax, or `None`
            to not apply such corrections.
        photoCalib : `None`
            Deprecated and ignored, only included for api compatibility.
        bboxToSpherePadding : `int`, optional
            Padding to account for translating a set of corners into a
            spherical (convex) boundary that is certain to encompase the
            enitre area covered by the box.

        Returns
        -------
        referenceCatalog : `lsst.afw.table.SimpleCatalog`
            Catalog containing reference objects inside the specified bounding
            box (padded by self.config.pixelMargin).

        Raises
        ------
        RuntimeError
            Raised if no reference catalogs could be found for the specified
            region.
        TypeError
            Raised if the loaded reference catalogs do not have matching
            schemas.
        """
        paddedBbox = geom.Box2D(bbox)
        paddedBbox.grow(self.config.pixelMargin)
        innerSkyRegion, outerSkyRegion, _, _ = self._makeBoxRegion(paddedBbox, wcs, bboxToSpherePadding)

        def _filterFunction(refCat, region):
            # Perform an initial "pre filter" step based on the refCat coords
            # and the outerSkyRegion created from the self.config.pixelMargin-
            # paddedBbox plus an "extra" padding of bboxToSpherePadding and the
            # raw wcs.  This should ensure a large enough projected area on the
            # sky that accounts for any projection/distortion issues, but small
            # enough to filter out loaded reference objects that lie well
            # beyond the projected detector of interest.  This step is required
            # due to the very local nature of the wcs available for the
            # sky <--> pixel conversions.
            preFiltFunc = _FilterCatalog(outerSkyRegion)
            refCat = preFiltFunc(refCat, region)

            # Add columns to the pre-filtered reference catalog relating their
            # coordinates to equivalent pixel positions for the wcs provided
            # and use to populate those columns.
            refCat = self.remapReferenceCatalogSchema(refCat, position=True)
            afwTable.updateRefCentroids(wcs, refCat)
            # No need to filter the catalog if it is entirely contained in the
            # region defined by the inner sky region.
            if innerSkyRegion.contains(region):
                return refCat
            # Create a new reference catalog, and populate it only with records
            # that fall inside the padded bbox.
            filteredRefCat = type(refCat)(refCat.table)
            centroidKey = afwTable.Point2DKey(refCat.schema['centroid'])
            for record in refCat:
                pixCoords = record[centroidKey]
                if paddedBbox.contains(geom.Point2D(pixCoords)):
                    filteredRefCat.append(record)
            return filteredRefCat
        return self.loadRegion(outerSkyRegion, filtFunc=_filterFunction, epoch=epoch, filterName=filterName)

    def loadRegion(self, region, filtFunc=None, filterName=None, epoch=None):
        """Load reference objects within a specified region.

        This function loads the DataIds used to construct an instance of this
        class which intersect or are contained within the specified region. The
        reference catalogs which intersect but are not fully contained within
        the input region are further filtered by the specified filter function.
        This function returns a single source catalog containing all reference
        objects inside the specified region.

        Parameters
        ----------
        region : `lsst.sphgeom.Region`
            This can be any type that is derived from `lsst.sphgeom.Region` and
            should define the spatial region for which reference objects are to
            be loaded.
        filtFunc : callable or `None`, optional
            This optional parameter should be a callable object that takes a
            reference catalog and its corresponding region as parameters,
            filters the catalog by some criteria and returns the filtered
            reference catalog. If `None`, an internal filter function is used
            which filters according to if a reference object falls within the
            input region.
        filterName : `str` or `None`, optional
            Name of camera filter, or `None` or blank for the default filter.
        epoch : `astropy.time.Time` or `None`, optional
            Epoch to which to correct proper motion and parallax, or `None` to
            not apply such corrections.

        Returns
        -------
        referenceCatalog : `lsst.afw.table.SourceCatalog`
            Catalog containing reference objects which intersect the input region,
            filtered by the specified filter function.

        Raises
        ------
        RuntimeError
            Raised if no reference catalogs could be found for the specified
            region.
        TypeError
            Raised if the loaded reference catalogs do not have matching
            schemas.
        """
        regionLat = region.getBoundingBox().getLat()
        regionLon = region.getBoundingBox().getLon()
        self.log.info("Loading reference objects from region bounded by "
                      "[{:.8f}, {:.8f}], [{:.8f}, {:.8f}] RA Dec".
                      format(regionLon.getA().asDegrees(), regionLon.getB().asDegrees(),
                             regionLat.getA().asDegrees(), regionLat.getB().asDegrees()))
        if filtFunc is None:
            filtFunc = _FilterCatalog(region)
        # filter out all the regions supplied by the constructor that do not overlap
        overlapList = []
        for dataId, refCat in zip(self.dataIds, self.refCats):
            # SphGeom supports some objects intersecting others, but is not symmetric,
            # try the intersect operation in both directions
            try:
                intersects = dataId.region.intersects(region)
            except TypeError:
                intersects = region.intersects(dataId.region)

            if intersects:
                overlapList.append((dataId, refCat))

        if len(overlapList) == 0:
            raise RuntimeError("No reference tables could be found for input region")

        firstCat = overlapList[0][1].get()
        refCat = filtFunc(firstCat, overlapList[0][0].region)
        trimmedAmount = len(firstCat) - len(refCat)

        # Load in the remaining catalogs
        for dataId, inputRefCat in overlapList[1:]:
            tmpCat = inputRefCat.get()

            if tmpCat.schema != firstCat.schema:
                raise TypeError("Reference catalogs have mismatching schemas")

            filteredCat = filtFunc(tmpCat, dataId.region)
            refCat.extend(filteredCat)
            trimmedAmount += len(tmpCat) - len(filteredCat)

        self.log.debug(f"Trimmed {trimmedAmount} refCat objects lying outside padded region, "
                       "leaving {len(refCat)}")
        self.log.info(f"Loaded {len(refCat)} reference objects")

        # Ensure that the loaded reference catalog is continuous in memory
        if not refCat.isContiguous():
            refCat = refCat.copy(deep=True)

        self.applyProperMotions(refCat, epoch)

        # Verify the schema is in the correct units and has the correct version; automatically convert
        # it with a warning if this is not the case.
        if not hasNanojanskyFluxUnits(refCat.schema) or not getFormatVersionFromRefCat(refCat) >= 1:
            self.log.warn("Found version 0 reference catalog with old style units in schema.")
            self.log.warn("run `meas_algorithms/bin/convert_refcat_to_nJy.py` to convert fluxes to nJy.")
            self.log.warn("See RFC-575 for more details.")
            refCat = convertToNanojansky(refCat, self.log)

        expandedCat = self.remapReferenceCatalogSchema(refCat, position=True)

        # Add flux aliases
        expandedCat = self.addFluxAliases(expandedCat, self.config.defaultFilter, self.config.filterMap)

        # Ensure that the returned reference catalog is continuous in memory
        if not expandedCat.isContiguous():
            expandedCat = expandedCat.copy(deep=True)

        fluxField = getRefFluxField(schema=expandedCat.schema, filterName=filterName)
        return pipeBase.Struct(refCat=expandedCat, fluxField=fluxField)

    def loadSkyCircle(self, ctrCoord, radius, filterName=None, epoch=None):
        """Load reference objects that lie within a circular region on the sky.

        This method constructs a circular region from an input center and
        angular radius, loads reference catalogs which are contained in or
        intersect the circle, and filters reference catalogs which intersect
        down to objects which lie within the defined circle.

        Parameters
        ----------
        ctrCoord : `lsst.geom.SpherePoint`
            Point defining the center of the circular region.
        radius : `lsst.geom.Angle`
            Defines the angular radius of the circular region.
        filterName : `str` or `None`, optional
            Name of camera filter, or `None` or blank for the default filter.
        epoch : `astropy.time.Time` or `None`, optional
            Epoch to which to correct proper motion and parallax, or `None` to
            not apply such corrections.

        Returns
        -------
        referenceCatalog : `lsst.afw.table.SourceCatalog`
            Catalog containing reference objects inside the specified search
            circle.
        """
        centerVector = ctrCoord.getVector()
        sphRadius = sphgeom.Angle(radius.asRadians())
        circularRegion = sphgeom.Circle(centerVector, sphRadius)
        return self.loadRegion(circularRegion, filterName=filterName, epoch=epoch)

    def joinMatchListWithCatalog(self, matchCat, sourceCat):
        """Relink an unpersisted match list to sources and reference objects.

        A match list is persisted and unpersisted as a catalog of IDs
        produced by afw.table.packMatches(), with match metadata
        (as returned by the astrometry tasks) in the catalog's metadata
        attribute. This method converts such a match catalog into a match
        list, with links to source records and reference object records.

        Parameters
        ----------
        matchCat : `lsst.afw.table.BaseCatalog`
            Unpersisted packed match list.
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
        return joinMatchListWithCatalogImpl(self, matchCat, sourceCat)

    def getMetadataBox(self, bbox, wcs, filterName=None, photoCalib=None, epoch=None,
                       bboxToSpherePadding=100):
        """Return metadata about the load

        This metadata is used for reloading the catalog (e.g., for
        reconstituting a normalised match list.)

        Parameters
        ----------
        bbox : `lsst.geom.Box2I` or `lsst.geom.Box2D`
            Bounding box for the pixels.
        wcs : `lsst.afw.geom.SkyWcs`
            The WCS object associated with ``bbox``.
        filterName : `str` or `None`, optional
            Name of the camera filter, or `None` or blank for the default
            filter.
        photoCalib : `None`
            Deprecated, only included for api compatibility.
        epoch : `astropy.time.Time` or `None`,  optional
            Epoch to which to correct proper motion and parallax, or `None` to
            not apply such corrections.
        bboxToSpherePadding : `int`, optional
            Padding to account for translating a set of corners into a
            spherical (convex) boundary that is certain to encompase the
            enitre area covered by the box.

        Returns
        -------
        md : `lsst.daf.base.PropertyList`
            The metadata detailing the search parameters used for this
            dataset.
        """
        paddedBbox = geom.Box2D(bbox)
        paddedBbox.grow(self.config.pixelMargin)
        _, _, innerCorners, outerCorners = self._makeBoxRegion(paddedBbox, wcs, bboxToSpherePadding)
        md = PropertyList()
        for box, corners in zip(("INNER", "OUTER"), (innerCorners, outerCorners)):
            for (name, corner) in zip(("UPPER_LEFT", "UPPER_RIGHT", "LOWER_LEFT", "LOWER_RIGHT"),
                                      corners):
                md.add(f"{box}_{name}_RA", geom.SpherePoint(corner).getRa().asDegrees(), f"{box}_corner")
                md.add(f"{box}_{name}_DEC", geom.SpherePoint(corner).getDec().asDegrees(), f"{box}_corner")
        md.add("SMATCHV", 1, 'SourceMatchVector version number')
        filterName = "UNKNOWN" if filterName is None else str(filterName)
        md.add('FILTER', filterName, 'filter name for photometric data')
        md.add('EPOCH', "NONE" if epoch is None else epoch.mjd, 'Epoch (TAI MJD) for catalog')
        return md

    @staticmethod
    def getMetadataCircle(coord, radius, filterName, photoCalib=None, epoch=None):
        """Return metadata about the load.

        This metadata is used for reloading the catalog (e.g. for
        reconstituting a normalized match list.)

        Parameters
        ----------
        coord : `lsst.geom.SpherePoint`
            ICRS center of the search region.
        radius : `lsst.geom.Angle`
            Radius of the search region.
        filterName : `str` or `None`
            Name of the camera filter, or `None` or blank for the default
            filter.
        photoCalib : `None`
            Deprecated, only included for api compatibility.
        epoch : `astropy.time.Time` or `None`, optional
            Epoch to which to correct proper motion and parallax, or `None` to
            not apply such corrections.

        Returns
        -------
        md : `lsst.daf.base.PropertyList`
        """
        md = PropertyList()
        md.add('RA', coord.getRa().asDegrees(), 'field center in degrees')
        md.add('DEC', coord.getDec().asDegrees(), 'field center in degrees')
        md.add('RADIUS', radius.asDegrees(), 'field radius in degrees, minimum')
        md.add('SMATCHV', 1, 'SourceMatchVector version number')
        filterName = "UNKNOWN" if filterName is None else str(filterName)
        md.add('FILTER', filterName, 'filter name for photometric data')
        md.add('EPOCH', "NONE" if epoch is None else epoch.mjd, 'Epoch (TAI MJD) for catalog')
        return md

    @staticmethod
    def addFluxAliases(refCat, defaultFilter, filterReferenceMap):
        """Add flux columns and aliases for camera to reference mapping.

        Creates a new catalog containing the information of the input refCat
        as well as added flux columns and aliases between camera and reference
        fluxes.

        Parameters
        ----------
        refCat : `lsst.afw.table.SimpleCatalog`
            Catalog of reference objects
        defaultFilter : `str`
            Name of the default reference filter
        filterReferenceMap : `dict` of `str`
            Dictionary with keys corresponding to a filter name and values
            which correspond to the name of the reference filter.

        Returns
        -------
        refCat : `lsst.afw.table.SimpleCatalog`
            Reference catalog with columns added to track reference filters.

        Raises
        ------
        `RuntimeError`
            If the specified reference filter name is not specifed as a
            key in the reference filter map.
        """
        refCat = ReferenceObjectLoader.remapReferenceCatalogSchema(refCat,
                                                                   filterNameList=filterReferenceMap.keys())
        aliasMap = refCat.schema.getAliasMap()
        if filterReferenceMap is None:
            filterReferenceMap = {}
        for filterName, refFilterName in itertools.chain([(None, defaultFilter)],
                                                         filterReferenceMap.items()):
            if refFilterName:
                camFluxName = filterName + "_camFlux" if filterName is not None else "camFlux"
                refFluxName = refFilterName + "_flux"
                if refFluxName not in refCat.schema:
                    raise RuntimeError("Unknown reference filter %s" % (refFluxName,))
                aliasMap.set(camFluxName, refFluxName)

                refFluxErrName = refFluxName + "Err"
                camFluxErrName = camFluxName + "Err"
                aliasMap.set(camFluxErrName, refFluxErrName)

        return refCat

    @staticmethod
    def remapReferenceCatalogSchema(refCat, *, filterNameList=None, position=False, photometric=False):
        """This function takes in a reference catalog and creates a new catalog with additional
        columns defined the remaining function arguments.

        Parameters
        ----------
        refCat : `lsst.afw.table.SimpleCatalog`
            Reference catalog to map to new catalog

        Returns
        -------
        expandedCat : `lsst.afw.table.SimpleCatalog`
            Deep copy of input reference catalog with additional columns added
        """
        mapper = afwTable.SchemaMapper(refCat.schema, True)
        mapper.addMinimalSchema(refCat.schema, True)
        mapper.editOutputSchema().disconnectAliases()
        if filterNameList:
            for filterName in filterNameList:
                mapper.editOutputSchema().addField(f"{filterName}_flux",
                                                   type=numpy.float64,
                                                   doc=f"flux in filter {filterName}",
                                                   units="Jy"
                                                   )
                mapper.editOutputSchema().addField(f"{filterName}_fluxErr",
                                                   type=numpy.float64,
                                                   doc=f"flux uncertanty in filter {filterName}",
                                                   units="Jy"
                                                   )

        if position:
            mapper.editOutputSchema().addField("centroid_x", type=float, doReplace=True)
            mapper.editOutputSchema().addField("centroid_y", type=float, doReplace=True)
            mapper.editOutputSchema().addField("hasCentroid", type="Flag", doReplace=True)
            mapper.editOutputSchema().getAliasMap().set("slot_Centroid", "centroid")

        if photometric:
            mapper.editOutputSchema().addField("photometric",
                                               type="Flag",
                                               doc="set if the object can be used for photometric"
                                                   "calibration",
                                               )
            mapper.editOutputSchema().addField("resolved",
                                               type="Flag",
                                               doc="set if the object is spatially resolved"
                                               )
            mapper.editOutputSchema().addField("variable",
                                               type="Flag",
                                               doc="set if the object has variable brightness"
                                               )

        expandedCat = afwTable.SimpleCatalog(mapper.getOutputSchema())
        expandedCat.setMetadata(refCat.getMetadata())
        expandedCat.extend(refCat, mapper=mapper)

        return expandedCat


def getRefFluxField(schema, filterName=None):
    """Get the name of a flux field from a schema.

    return the alias of "anyFilterMapsToThis", if present
    else if filterName is specified:
        return "*filterName*_camFlux" if present
        else return "*filterName*_flux" if present (camera filter name
            matches reference filter name)
        else throw RuntimeError
    else:
        return "camFlux", if present,
        else throw RuntimeError

    Parameters
    ----------
    schema : `lsst.afw.table.Schema`
        Reference catalog schema.
    filterName : `str`, optional
        Name of camera filter. If not specified, ``defaultFilter`` needs to be
        set in the refcat loader config.

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
    try:
        return schema.getAliasMap().get("anyFilterMapsToThis")
    except LookupError:
        pass  # try the filterMap next

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
        default=250,
        min=0,
    )
    defaultFilter = pexConfig.Field(
        doc=("Default reference catalog filter to use if filter not specified in exposure;"
             " if blank then filter must be specified in exposure."),
        dtype=str,
        default="",
        deprecated="defaultFilter is deprecated by RFC-716. Will be removed after v22."
    )
    anyFilterMapsToThis = pexConfig.Field(
        doc=("Always use this reference catalog filter, no matter whether or what filter name is "
             "supplied to the loader. Effectively a trivial filterMap: map all filter names to this filter."
             " This can be set for purely-astrometric catalogs (e.g. Gaia DR2) where there is only one "
             "reasonable choice for every camera filter->refcat mapping, but not for refcats used for "
             "photometry, which need a filterMap and/or colorterms/transmission corrections."),
        dtype=str,
        default=None,
        optional=True
    )
    filterMap = pexConfig.DictField(
        doc=("Mapping of camera filter name: reference catalog filter name; "
             "each reference filter must exist in the refcat."
             " Note that this does not perform any bandpass corrections: it is just a lookup."),
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

    def validate(self):
        super().validate()
        if self.filterMap != {} and self.anyFilterMapsToThis is not None:
            msg = "`filterMap` and `anyFilterMapsToThis` are mutually exclusive"
            raise pexConfig.FieldValidationError(LoadReferenceObjectsConfig.anyFilterMapsToThis,
                                                 self, msg)


class LoadReferenceObjectsTask(pipeBase.Task, ReferenceObjectLoaderBase, metaclass=abc.ABCMeta):
    """Abstract base class to load objects from reference catalogs.
    """
    ConfigClass = LoadReferenceObjectsConfig
    _DefaultName = "LoadReferenceObjects"

    def __init__(self, butler=None, *args, **kwargs):
        """Construct a LoadReferenceObjectsTask

        Parameters
        ----------
        butler : `lsst.daf.persistence.Butler`
            Data butler, for access reference catalogs.
        """
        pipeBase.Task.__init__(self, *args, **kwargs)
        self.butler = butler

    @pipeBase.timeMethod
    def loadPixelBox(self, bbox, wcs, filterName=None, photoCalib=None, epoch=None):
        """Load reference objects that overlap a rectangular pixel region.

        Parameters
        ----------
        bbox : `lsst.geom.Box2I` or `lsst.geom.Box2D`
            Bounding box for pixels.
        wcs : `lsst.afw.geom.SkyWcs`
            WCS; used to convert pixel positions to sky coordinates
            and vice-versa.
        filterName : `str` or `None`, optional
            Name of filter, or `None` or `""` for the default filter.
            This is used for flux values in case we have flux limits
            (which are not yet implemented).
        photoCalib : `None`
            Deprecated, only included for api compatibility.
        epoch : `astropy.time.Time` or `None`, optional
            Epoch to which to correct proper motion and parallax, or `None` to
            not apply such corrections.

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
        loadRes = self.loadSkyCircle(circle.coord, circle.radius, filterName=filterName, epoch=epoch,
                                     centroids=True)
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
    def loadSkyCircle(self, ctrCoord, radius, filterName=None, epoch=None, centroids=False):
        """Load reference objects that overlap a circular sky region.

        Parameters
        ----------
        ctrCoord : `lsst.geom.SpherePoint`
            ICRS center of search region.
        radius : `lsst.geom.Angle`
            Radius of search region.
        filterName : `str` or `None`, optional
            Name of filter, or `None` or `""` for the default filter.
            This is used for flux values in case we have flux limits
            (which are not yet implemented).
        epoch : `astropy.time.Time` or `None`, optional
            Epoch to which to correct proper motion and parallax, or `None` to
            not apply such corrections.
        centroids : `bool`, optional
            Add centroid fields to the loaded Schema. ``loadPixelBox`` expects
            these fields to exist.

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

        Returns
        -------
        catalog : `lsst.afw.table.SimpleCatalog`
            Reference objects in the bbox, with centroid and
            hasCentroid fields set.
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

        Raises
        ------
        RuntimeError
            If any reference flux field is missing from the schema.
        """
        aliasMap = schema.getAliasMap()

        if self.config.anyFilterMapsToThis is not None:
            refFluxName = self.config.anyFilterMapsToThis + "_flux"
            if refFluxName not in schema:
                msg = f"Unknown reference filter for anyFilterMapsToThis='{refFluxName}'"
                raise RuntimeError(msg)
            aliasMap.set("anyFilterMapsToThis", refFluxName)
            return  # this is mutually exclusive with filterMap

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
    def makeMinimalSchema(filterNameList, *, addCentroid=False,
                          addIsPhotometric=False, addIsResolved=False,
                          addIsVariable=False, coordErrDim=2,
                          addProperMotion=False, properMotionErrDim=2,
                          addParallax=False):
        """Make a standard schema for reference object catalogs.

        Parameters
        ----------
        filterNameList : `list` of `str`
            List of filter names. Used to create <filterName>_flux fields.
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
                units="nJy",
            )
        for filterName in filterNameList:
            schema.addField(
                field="%s_fluxErr" % (filterName,),
                type=numpy.float64,
                doc="flux uncertainty in filter %s" % (filterName,),
                units="nJy",
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
                type="Angle",
                doc="parallax",
                units="rad",
            )
            schema.addField(
                field="parallaxErr",
                type="Angle",
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
            - bbox : `lsst.geom.Box2D`
                Bounding box used to compute the circle.
        """
        bbox = geom.Box2D(bbox)  # make sure bbox is double and that we have a copy
        bbox.grow(self.config.pixelMargin)
        coord = wcs.pixelToSky(bbox.getCenter())
        radius = max(coord.separation(wcs.pixelToSky(pp)) for pp in bbox.getCorners())
        return pipeBase.Struct(coord=coord, radius=radius, bbox=bbox)

    def getMetadataBox(self, bbox, wcs, filterName=None, photoCalib=None, epoch=None):
        """Return metadata about the load.

        This metadata is used for reloading the catalog (e.g., for
        reconstituting a normalised match list.

        Parameters
        ----------
        bbox : `lsst.geom.Box2I` or `lsst.geom.Box2D`
            Pixel bounding box.
        wcs : `lsst.afw.geom.SkyWcs`
            WCS; used to convert pixel positions to sky coordinates.
        filterName : `str` or `None`, optional
            Name of camera filter, or `None` or `""` for the default
            filter.
        photoCalib : `None`
            Deprecated, only included for api compatibility.
        epoch : `astropy.time.Time` or `None`, optional
            Epoch to which to correct proper motion and parallax,
            or None to not apply such corrections.

        Returns
        -------
        metadata : `lsst.daf.base.PropertyList`
            Metadata about the load.
        """
        circle = self._calculateCircle(bbox, wcs)
        return self.getMetadataCircle(circle.coord, circle.radius, filterName, epoch=epoch)

    def getMetadataCircle(self, coord, radius, filterName, photoCalib=None, epoch=None):
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
        photoCalib : `None`
            Deprecated, only included for api compatibility.
        epoch : `astropy.time.Time` (optional)
            Epoch to which to correct proper motion and parallax, or `None` to
            not apply such corrections.

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
        md.add('EPOCH', "NONE" if epoch is None else epoch.mjd, 'Epoch (TAI MJD) for catalog')
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
        return joinMatchListWithCatalogImpl(self, matchCat, sourceCat)


def joinMatchListWithCatalogImpl(refObjLoader, matchCat, sourceCat):
    """Relink an unpersisted match list to sources and reference
    objects.

    A match list is persisted and unpersisted as a catalog of IDs
    produced by afw.table.packMatches(), with match metadata
    (as returned by the astrometry tasks) in the catalog's metadata
    attribute. This method converts such a match catalog into a match
    list, with links to source records and reference object records.

    Parameters
    ----------
    refObjLoader
        Reference object loader to use in getting reference objects
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
    try:
        epoch = matchmeta.getDouble('EPOCH')
    except (LookupError, TypeError):
        epoch = None  # Not present, or not correct type means it's not set
    if 'RADIUS' in matchmeta:
        # This is a circle style metadata, call loadSkyCircle
        ctrCoord = geom.SpherePoint(matchmeta.getDouble('RA'),
                                    matchmeta.getDouble('DEC'), geom.degrees)
        rad = matchmeta.getDouble('RADIUS')*geom.degrees
        refCat = refObjLoader.loadSkyCircle(ctrCoord, rad, filterName, epoch=epoch).refCat
    elif "INNER_UPPER_LEFT_RA" in matchmeta:
        # This is the sky box type (only triggers in the LoadReferenceObject class, not task)
        # Only the outer box is required to be loaded to get the maximum region, all filtering
        # will be done by the unpackMatches function, and no spatial filtering needs to be done
        # by the refObjLoader
        box = []
        for place in ("UPPER_LEFT", "UPPER_RIGHT", "LOWER_LEFT", "LOWER_RIGHT"):
            coord = geom.SpherePoint(matchmeta.getDouble(f"OUTER_{place}_RA"),
                                     matchmeta.getDouble(f"OUTER_{place}_DEC"),
                                     geom.degrees).getVector()
            box.append(coord)
        outerBox = sphgeom.ConvexPolygon(box)
        refCat = refObjLoader.loadRegion(outerBox, filterName=filterName, epoch=epoch).refCat

    refCat.sort()
    sourceCat.sort()
    return afwTable.unpackMatches(matchCat, refCat, sourceCat)


def applyProperMotionsImpl(log, catalog, epoch):
    """Apply proper motion correction to a reference catalog.

    Adjust position and position error in the ``catalog``
    for proper motion to the specified ``epoch``,
    modifying the catalog in place.

    Parameters
    ----------
    log : `lsst.log.Log`
        Log object to write to.
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
    epoch : `astropy.time.Time`
        Epoch to which to correct proper motion.
    """
    if "epoch" not in catalog.schema or "pm_ra" not in catalog.schema or "pm_dec" not in catalog.schema:
        log.warn("Proper motion correction not available from catalog")
        return
    if not catalog.isContiguous():
        raise RuntimeError("Catalog must be contiguous")
    catEpoch = astropy.time.Time(catalog["epoch"], scale="tai", format="mjd")
    log.info("Correcting reference catalog for proper motion to %r", epoch)
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
                   record.get(coordKey).offset(bearing=bearingRad*geom.radians,
                                               amount=amountRad*geom.radians))
    # Increase error in RA and Dec based on error in proper motion
    if "coord_raErr" in catalog.schema:
        catalog["coord_raErr"] = numpy.hypot(catalog["coord_raErr"],
                                             catalog["pm_raErr"]*timeDiffsYears)
    if "coord_decErr" in catalog.schema:
        catalog["coord_decErr"] = numpy.hypot(catalog["coord_decErr"],
                                              catalog["pm_decErr"]*timeDiffsYears)
