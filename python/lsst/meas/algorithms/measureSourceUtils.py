"""Support utilities for Measuring sources"""

import re
import lsst.meas.algorithms as algorithms

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
    
