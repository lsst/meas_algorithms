# Config file for tests/measure.py; previously tests/MeasureSources.paf.

# 'root' should be an instance of lsst.meas.algorithms.MeasureSourcesConfig (defined in measurement.py)

root.slots.centroid = "centroid.naive"
root.slots.apFlux = "flux.naive"
root.slots.modelFlux = "flux.gaussian"
root.slots.psfFlux = "flux.psf"
root.algorithms["flux.naive"].radius = 3.0
root.algorithms["flux.gaussian"].shiftmax = 10
root.algorithms["flux.sinc"].radius = 3.0
root.algorithms.names = ["flags.pixel",
                         "centroid.gaussian", "centroid.naive",
                         "shape.sdss",
                         "flux.naive", "flux.gaussian", "flux.psf", "flux.sinc",
                         "classification.extendedness"]
root.centroider.name = "centroid.sdss"
