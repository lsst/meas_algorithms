import lsst.pex.config as pexConf
from . import algorithmsLib

sourceMeasurementRegistry = pexConf.makeRegistry("Registry for all source measurement classes.")

class SourceSlotConfig(pexConf.Config):

    centroid = pexConf.Field(dtype=str, default="centroid.sdss", optional=True,
                             doc="the name of the centroiding algorithm used to set source (x,y)")
    shape = pexConf.Field(dtype=str, default="shape.sdss", optional=True,
                          doc="the name of the algorithm used to set source ellipticity and size parameters")
    apFlux = pexConf.Field(dtype=str, default="flux.sinc", optional=True,
                           doc="the name of the algorithm used to set the source aperture flux slot")
    modelFlux = pexConf.Field(dtype=str, default="flux.gaussian", optional=True,
                           doc="the name of the algorithm used to set the source model flux slot")
    psfFlux = pexConf.Field(dtype=str, default="flux.psf", optional=True,
                           doc="the name of the algorithm used to set the source psf flux slot")
    instFlux = pexConf.Field(dtype=str, default="flux.gaussian", optional=True,
                             doc="the name of the algorithm used to set the source inst flux slot")

    def apply(self, table):
        if self.centroid is not None: table.defineCentroid(self.centroid)
        if self.shape is not None: table.defineShape(self.shape)
        if self.apFlux is not None: table.defineApFlux(self.apFlux)
        if self.modelFlux is not None: table.defineModelFlux(self.modelFlux)
        if self.psfFlux is not None: table.definePsfFlux(self.psfFlux)
        if self.instFlux is not None: table.defineInstFlux(self.instFlux)

@pexConf.wrap(algorithmsLib.ClassificationControl)
class ClassificationConfig(pexConf.Config):
    pass

class MeasureSourcesConfig(pexConf.Config):

    source = pexConf.ConfigField("The mapping from algorithms to fields in Source", SourceSlotConfig)

    algorithms = sourceMeasurementRegistry.makeField(
        "Configurations for individual algorithms", multi=True,
        default=["flags.pixel",
                 "centroid.gaussian", "centroid.naive", "centroid.sdss",
                 "shape.sdss",
                 "flux.gaussian", "flux.naive", "flux.psf", "flux.sinc",
                 "extendedness",
                 ]
        )

    def validate(self):
        pexConf.Config.validate(self)
        if self.source.centroid is not None and self.source.centroid not in self.algorithms.names:
            raise ValueError("source centroid slot algorithm '%s' is not being run." % self.source.astrom)
        if self.source.shape is not None and self.source.shape not in self.algorithms.names:
            raise ValueError("source shape slot algorithm '%s' is not being run." % self.source.shape)
        for slot in (self.source.psfFlux, self.source.apFlux, self.source.modelFlux, self.source.instFlux):
            if slot is not None and slot not in self.algorithms.names:
                raise ValueError("source flux slot algorithm '%s' is not being run." % slot)

    def makeMeasureSources(self, schema=None):
        if schema is None:
            ms = algorithmsLib.MeasureSources()
        else:
            ms = algorithmsLib.MeasureSources(schema)
        ms.addAlgorithms(self.algorithms.apply())
        return ms

def declareMeasurement(control):
    """Declare a measurement algorithm from a C++ Control Class.

    @param control   A wrapped C++ control object that will be used to define the Config.

    This registers the control/config pair, using the config.makeControl method
    (automatically created) as the factory.  Calling applyFactory() on a RegistryField that uses
    this registry thus returns Control object instance, filled with the settings from the active
    Config object.

    The config class is added to the current scope, with its name set by replacing "Control"
    with "Config" in the __name__ of the control class.  In addition, both config and control
    classes will be given a 'name' attribute with the registered name, since each measurement
    algorithm is registered only once.
    """
    global sourceMeasurementRegistry
    ctrlInstance = control()
    name = ctrlInstance.name
    config = pexConf.makeConfigClass(control, base=pexConf.Config)
    sourceMeasurementRegistry.register(name, config.makeControl, config)
    config.name = name
    globals()[config.__name__] = config

declareMeasurement(algorithmsLib.ClassificationControl)
declareMeasurement(algorithmsLib.PixelFlagControl)

declareMeasurement(algorithmsLib.GaussianCentroidControl)
declareMeasurement(algorithmsLib.NaiveCentroidControl)
declareMeasurement(algorithmsLib.SdssCentroidControl)

declareMeasurement(algorithmsLib.SdssShapeControl)

declareMeasurement(algorithmsLib.ApertureFluxControl)
declareMeasurement(algorithmsLib.GaussianFluxControl)
declareMeasurement(algorithmsLib.NaiveFluxControl)
declareMeasurement(algorithmsLib.PsfFluxControl)

# Here's an example on how to declare a measurement config more manually, and add a property to the Config.
@pexConf.wrap(algorithmsLib.SincFluxControl)
class SincFluxConfig(pexConf.Config):
    name = "SINC"
    def _get_radius(self): return self.radius2
    def _set_radius(self, r): self.radius2 = r
    radius = property(_get_radius, _set_radius, doc="synonym for radius2")
sourceMeasurementRegistry.register("flux.sinc", SincFluxConfig.makeControl, SincFluxConfig)
