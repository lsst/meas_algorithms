# Config file for tests/measure.py; previously tests/MeasureSources.paf.

root.slots.centroid = "centroid.naive"
root.slots.apFlux = "flux.naive"
root.slots.modelFlux = "flux.gaussian"
root.slots.psfFlux = "flux.psf"
root.measurement.algorithms["flux.naive"].radius = 3.0
root.measurement.algorithms["flux.gaussian"].shiftmax = 10
root.measurement.algorithms["flux.sinc"].radius = 3.0
root.measurement.algorithms.names = ["flags.pixel",
                                     "centroid.gaussian", "centroid.naive",
                                     "shape.sdss",
                                     "flux.naive", "flux.gaussian", "flux.psf", "flux.sinc",
                                     "classification.extendedness"]
root.measurement.centroider.name = "centroid.sdss"
