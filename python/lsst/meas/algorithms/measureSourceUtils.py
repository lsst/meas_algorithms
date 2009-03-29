"""Support utilities for Measuring sources"""

import re
import lsst.meas.algorithms as algorithms
import lsst.afw.display.ds9 as ds9

def explainDetectionFlags(flags):
    """Return a string explaining Source's detectionFlags"""

    result = []
    for k in algorithms.Flags.__dict__.keys():
        if not re.search(r"^[_A-Z0-9]+$", k): # flag names match this re
            continue

        if (flags & algorithms.Flags.__dict__[k]):
            result += [k]

    result.sort()
    return " ".join(result)
    
def showSourceSet(sSet, xy0=(0, 0), frame=0, ctype=ds9.GREEN):
    """Draw the (XAstrom, YAstrom) positions of a set of Sources.  Image has the given XY0"""
    for s in sSet:
        ds9.dot("+", s.getXAstrom() - xy0[0], s.getYAstrom() - xy0[1], frame=frame, ctype=ctype)
