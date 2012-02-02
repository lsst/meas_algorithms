# Config file for tests/measure.py; previously tests/MeasureSources.paf.

root.source.centroid = "centroid.naive"
root.source.apFlux = "flux.naive"
root.source.modelFlux = "flux.gaussian"
root.source.psfFlux = "flux.psf"
root.algorithms["flux.naive"].radius = 3.0
root.algorithms["flux.gaussian"].shiftmax = 10
root.algorithms["flux.sinc"].radius = 3.0
root.algorithms.names = ["centroid.gaussian", "centroid.naive", "centroid.sdss",
                         "shape.sdss",
                         "flux.naive", "flux.gaussian", "flux.psf", "flux.sinc"]
