#!/usr/bin/env python
import numpy
import lsst.daf.base
from matplotlib import pyplot

def makeArrays(dataref):
    metadata = dataref.get("processCcd_metadata")
    timings = metadata.get("processCcd:measurement.algorithmTimings")
    arrays = {}
    for name in timings.names():
        arrays[name] = numpy.array(timings.get(name), dtype=float)
    return arrays

def plotSingle(arrays):
    pyplot.figure()
    for name, array in arrays.iteritems():
        pyplot.semilogy(array, alpha=0.7, label=name)
    pyplot.legend()

def plotComparison(arrays1, arrays2):
    pyplot.figure()
    shared = arrays1.viewkeys().intersection(arrays2.viewkeys())
    for name in shared:
        pyplot.loglog(arrays1[name], arrays2[name], ',', alpha=0.5, label=name)
    

if __name__ == "__main__":
    import lsst.daf.persistence
    import lsst.obs.lsstSim
    import sys
    b1 = lsst.daf.persistence.ButlerFactory(mapper=lsst.obs.lsstSim.LsstSimMapper(root=sys.argv[1])).create()
    arrays1 = makeArrays(b1)
    if len(sys.argv) == 3:
        b2 = lsst.daf.persistence.ButlerFactory(mapper=lsst.obs.lsstSim.LsstSimMapper(root=sys.argv[2])).create()
        arrays2 = makeArrays(b2)
        plotComparison(arrays1, arrays1)
    else:
        plotSingle(arrays1)
    pyplot.show()
