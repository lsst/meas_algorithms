# Config file for tests/measure.py; previously tests/MeasureSources.paf.

# 'config' should be an instance of lsst.meas.algorithms.SourceMeasurementConfig (defined in measurement.py)

config.slots.centroid = "centroid.naive"
config.slots.apFlux = "flux.naive"
config.slots.modelFlux = "flux.gaussian"
config.slots.psfFlux = "flux.psf"
config.slots.shape = "shape.sdss"
config.algorithms["flux.naive"].radius = 3.0
config.algorithms["flux.gaussian"].shiftmax = 10.0
config.algorithms["flux.sinc"].radius = 3.0
config.algorithms.names = ["flags.pixel",
                         "centroid.gaussian", "centroid.naive",
                         "shape.sdss",
                         "flux.naive", "flux.gaussian", "flux.psf", "flux.sinc",
                         "classification.extendedness"]
config.centroider.name = "centroid.sdss"
