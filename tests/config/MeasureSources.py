# Config file for tests/measure.py; previously tests/MeasureSources.paf.

import lsst.meas.algorithms as measAlg

root.source.astrom = "NAIVE"
root.source.apFlux = "NAIVE"
root.source.modelFlux = "GAUSSIAN"
root.source.psfFlux = "PSF"
root.source.shape = "SDSS"
root.astrometry.names = ["GAUSSIAN", "NAIVE", "SDSS"]
root.shape.names = ["SDSS"]
root.photometry.names = ["NAIVE", "GAUSSIAN", "PSF", "SINC"]
root.photometry["NAIVE"].radius = 3.0
root.photometry["GAUSSIAN"].shiftmax = 10
root.photometry["SINC"].radius = 3.0
