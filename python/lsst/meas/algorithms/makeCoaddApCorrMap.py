from lsst.afw.image import ApCorrMap
from . import CoaddBoundedField, CoaddBoundedFieldElement

__all__ = ["makeCoaddApCorrMap", ]


def makeCoaddApCorrMap(catalog, coaddBox, coaddWcs, weightFieldName="weight"):
    """Construct an ApCorrMap for a coadd

    @param catalog: Table of coadd inputs (lsst.afw.table.ExposureCatalog)
    @param coaddBox: Bounding box for coadd (lsst.afw.geom.Box2I)
    @param coaddWcs: Wcs for coadd
    @param weightFieldName: name of weight field in catalog
    @return aperture corrections
    """

    # Assemble the BoundedFields for each type
    everything = {} # name --> list of CoaddBoundedFieldElement
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
