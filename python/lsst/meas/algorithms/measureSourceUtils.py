"""Support utilities for Measuring sources"""

import re
import lsst.afw.detection as afwDet
import algorithmsLib as algorithmsLib
import lsst.afw.display.ds9 as ds9

def explainDetectionFlags(flags):
    """Return a string explaining Source's detectionFlags"""

    result = []
    for k in algorithmsLib.Flags.__dict__.keys():
        if not re.search(r"^[_A-Z0-9]+$", k): # flag names match this re
            continue

        if (flags & algorithmsLib.Flags.__dict__[k]):
            result += [k]

    result.sort()
    return " ".join(result)
    
def showSourceSet(sSet, xy0=(0, 0), frame=0, ctype=ds9.GREEN, symb="+", size=2):
    """Draw the (XAstrom, YAstrom) positions of a set of Sources.  Image has the given XY0"""
    for s in sSet:
        ds9.dot(symb, s.getXAstrom() - xy0[0], s.getYAstrom() - xy0[1], frame=frame, ctype=ctype, size=size)

def makeSource(measureSources, fp, sourceId, isNegative):
    source = afwDet.Source()
    source.setId(sourceId)
    detectionBits = source.getFlagForDetection() | algorithmsLib.Flags.BINNED1
    if isNegative:
        detectionBits |= algorithmsLib.Flags.DETECT_NEGATIVE
    source.setFlagForDetection(detectionBits)
    try:
        measureSources.apply(source, fp)
    except Exception, e:
        # don't worry about measurement exceptions
        pass            

    return source

def makeSourceSet(measureSources, fpPositive, fpNegative):
    """
    Apply the measureSources functor to each footprint 
    - measureSources: instance of an lsst.meas.algorithms.MeasureSources
    - fpPositive: list of positive lsst.afw.detection.Footprint
    - fpNegative: list of negative lsst.afw.detection Footprint
    """
    sourceSet = afwDet.SourceSet()
    sourceId = 0;
    if fpPositive != None:
        for fp in fpPositive:
            sourceSet.append(makeSource(measureSources, fp, sourceId, False))
            sourceId += 1
    if fpNegative != None:    
        for fp in fpNegative:
            sourceSet.append(makeSource(measureSources, fp, sourceId, True))
            sourceId += 1

    return sourceSet
