import lsst.pex.config as pexConf
from . import algorithmsLib

class registries: # class is really just a namespace; it will go away with new Source

    astrometry = pexConf.makeRegistry("Registry for all astrometry measurement classes.")

    photometry = pexConf.makeRegistry("Registry for all photometry measurement classes.")

    shape = pexConf.makeRegistry("Registry for all shape measurement classes.")

def declareMeasurement(name, registry, control):
    """Declare a measurement algorithm from a C++ Control Class.

    @param name      The name of the algorithm for use in RegistryFields (typically all-caps).
    @param registry  The registry in which to put the algorithm.
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
    config = pexConf.makeConfigClass(control, base=pexConf.Config)
    registry.register(name, config.makeControl, config)
    config.name = name
    control.name = name
    globals()[config.__name__] = config

declareMeasurement("GAUSSIAN", registries.astrometry, algorithmsLib.GaussianAstrometryControl)
declareMeasurement("NAIVE", registries.astrometry, algorithmsLib.NaiveAstrometryControl)
declareMeasurement("SDSS", registries.astrometry, algorithmsLib.SdssAstrometryControl)

declareMeasurement("SDSS", registries.shape, algorithmsLib.SdssShapeControl)

declareMeasurement("APERTURE", registries.photometry, algorithmsLib.AperturePhotometryControl)
declareMeasurement("GAUSSIAN", registries.photometry, algorithmsLib.GaussianPhotometryControl)
declareMeasurement("NAIVE", registries.photometry, algorithmsLib.NaivePhotometryControl)
declareMeasurement("PSF", registries.photometry, algorithmsLib.PsfPhotometryControl)

# Here's an example on how to declare a measurement config more manually, and add a property to the Config.
@pexConf.wrap(algorithmsLib.SincPhotometryControl)
class SincPhotometryConfig(pexConf.Config):
    name = "SINC"
    def _get_radius(self): return self.radius2
    def _set_radius(self, r): self.radius2 = r
    radius = property(_get_radius, _set_radius, doc="synonym for radius2")
registries.photometry.register("SINC", SincPhotometryConfig.makeControl, SincPhotometryConfig)


class SourceConfig(pexConf.Config):
    astrom = pexConf.Field("The name of the centroiding algorithm used to set Source.[XY]Astrom",
                           dtype=str, default="SDSS", optional=True)
    shape = pexConf.Field("The name of the centroiding algorithm used to set Source.Mxx etc.",
                          dtype=str, default="SDSS", optional=True)
    apFlux = pexConf.Field("The name of the algorithm used to set Source.apFlux(Err)",
                           dtype=str, default="SINC", optional=True)
    modelFlux = pexConf.Field("The name of the algorithm used to set Source.modelFlux(Err)",
                              dtype=str, default="GAUSSIAN", optional=True)
    psfFlux = pexConf.Field("The name of the algorithm used to set Source.psfFlux(Err)",
                            dtype=str, default="PSF", optional=True)
    instFlux = pexConf.Field("The name of the algorithm used to set Source.instFlux(Err)",
                             dtype=str, default="GAUSSIAN", optional=True)

class ClassificationConfig(pexConf.Config):
    sg_fac1 = pexConf.Field("First S/G parameter; critical ratio of inst to psf flux", dtype=float, 
                            default=0.925, optional=True)
    sg_fac2 = pexConf.Field("Second S/G parameter; correction for instFlux error", dtype=float,
                            default=0.0, optional=True)
    sg_fac3 = pexConf.Field("Third S/G parameter; correction for psfFlux error", dtype=float,
                            default=0.0, optional=True)

class MeasureSourcesConfig(pexConf.Config):

    source = pexConf.ConfigField("The mapping from algorithms to fields in Source", SourceConfig)

    astrometry = registries.astrometry.makeField("Configurations for individual astrometry algorithms.",
                                       multi=True)
    astrometry.defaults = ["GAUSSIAN", "NAIVE", "SDSS"]

    shape = registries.shape.makeField("Configurations for various shape-measurement algorithms.",
                                       multi=True)
    shape.defaults = ["SDSS"]
        
    photometry = registries.photometry.makeField("Configurations for individual photometry algorithms.",
                                                 multi=True)
    photometry.defaults = ["GAUSSIAN", "NAIVE", "PSF", "SINC"]

    classification = pexConf.ConfigField("Parameters to do with star/galaxy classification",
                                         ClassificationConfig)

    def __init__(self, *args, **kwds):
        pexConf.Config.__init__(self, *args, **kwds)
        self.astrometry.names = type(self).astrometry.defaults
        self.shape.names = type(self).shape.defaults
        self.photometry.names = type(self).photometry.defaults

    def validate(self):
        pexConf.Config.validate(self)
        if self.source.astrom is not None and self.source.astrom not in self.astrometry.names:
            raise ValueError("source astrometry slot algorithm '%s' is not being run." % self.source.astrom)
        if self.source.shape is not None and self.source.shape not in self.shape.names:
            raise ValueError("source shape slot algorithm '%s' is not being run." % self.source.shape)
        for slot in (self.source.psfFlux, self.source.apFlux, self.source.modelFlux, self.source.instFlux):
            if slot is not None and slot not in self.photometry.names:
                raise ValueError("source photometry slot algorithm '%s' is not being run." % slot)

    def makeMeasureSources(self, exposure):
        import lsst.pex.policy
        from lsst.pex.config.convert import makePolicy
        self.validate()
        policy = lsst.pex.policy.Policy()
        policy.set("source", makePolicy(self.source))
        policy.set("classification", makePolicy(self.classification))
        ms = algorithmsLib.makeMeasureSources(exposure, policy)
        ms.getMeasureAstrom().addAlgorithms(self.astrometry.apply())
        ms.getMeasureShape().addAlgorithms(self.shape.apply())
        ms.getMeasurePhotom().addAlgorithms(self.photometry.apply())
        return ms

