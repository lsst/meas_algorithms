import lsst.pex.config as pexConf
from . import algorithmsLib

class registries: # class is really just a namespace; it will go away with new Source

    astrometry = pexConf.makeConfigRegistry("Registry of all astrometry measurement Config classes.")

    photometry = pexConf.makeConfigRegistry("Registry of all photometry measurement Config classes.")

    shape = pexConf.makeConfigRegistry("Registry of all shape measurement Config classes.")

class SourceConfig(pexConf.Config):
    astrom = pexConf.Field("The name of the centroiding algorithm used to set Source.[XY]Astrom", str)
    shape = pexConf.Field("The name of the centroiding algorithm used to set Source.Mxx etc.", str)
    apFlux = pexConf.Field("The name of the algorithm used to set Source.apFlux(Err)", str)
    modelFlux = pexConf.Field("The name of the algorithm used to set Source.modelFlux(Err)", str)
    psfFlux = pexConf.Field("The name of the algorithm used to set Source.psfFlux(Err)", str)
    instFlux = pexConf.Field("The name of the algorithm used to set Source.instFlux(Err)", str)

class ClassificationConfig(pexConf.Config):
    sg_fac1 = pexConf.Field("First S/G parameter; critical ratio of inst to psf flux", dtype=float, 
                            default=0.925, optional=True)
    sg_fac2 = pexConf.Field("Second S/G parameter; correction for instFlux error", dtype=float,
                            default=0.0, optional=True)
    sg_fac3 = pexConf.Field("Third S/G parameter; correction for psfFlux error", dtype=float,
                            default=0, optional=True)

class MeasureSourcesConfig(pexConf.Config):

    source = pexConf.ConfigField("The mapping from algorithms to fields in Source", SourceConfig)

    astrometry = pexConf.RegistryField("Configurations for individual astrometry algorithms.",
                                       typemap=registries.astrometry, multi=True)
    astrometry.defaults = ["GAUSSIAN", "NAIVE", "SDSS"]

    shape = pexConf.RegistryField("Configurations for various shape-measurement algorithms.",
                                  typemap=registries.shape, multi=True)
    shape.defaults = ["SDSS"]
        
    photometry = pexConf.RegistryField("Configurations for individual photometry algorithms.",
                               typemap=registries.photometry, multi=True)
    photometry.defaults = ["GAUSSIAN", "NAIVE", "PSF", "SINC"]

    classification = pexConf.ConfigField("Parameters to do with star/galaxy classification",
                                         ClassificationConfig)

    def __init__(self, *args, **kwds):
        pexConf.Config.__init__(self, *args, **kwds)
        self.astrometry.names = type(self).astrometry.defaults
        self.shape.names = type(self).shape.defaults
        self.photometry.names = type(self).photometry.defaults

    def makeMeasureSources(self, exposure):
        import lsst.pex.policy
        from lsst.pex.config.convert import makePolicy
        policy = lsst.pex.policy.Policy()
        policy.set("source", makePolicy(self.source))
        policy.set("classification", makePolicy(self.classification))
        ms = algorithmsLib.makeMeasureSources(exposure, policy)
        for conf in self.astrometry.active:
            ms.getMeasureAstrom().addAlgorithm(conf.makeControl())
        for conf in self.shape.active:
            ms.getMeasureShape().addAlgorithm(conf.makeControl())
        for conf in self.photometry.active:
            ms.getMeasurePhotom().addAlgorithm(conf.makeControl())
        return ms

@pexConf.register("GAUSSIAN", registries.astrometry)
@pexConf.wrap(algorithmsLib.GaussianAstrometryControl)
class GaussianAstrometryConfig(pexConf.Config): pass

@pexConf.register("NAIVE", registries.astrometry)
@pexConf.wrap(algorithmsLib.NaiveAstrometryControl)
class NaiveAstrometryConfig(pexConf.Config): pass

@pexConf.register("SDSS", registries.astrometry)
@pexConf.wrap(algorithmsLib.SdssAstrometryControl)
class SdssAstrometryConfig(pexConf.Config): pass

@pexConf.register("SDSS", registries.shape)
@pexConf.wrap(algorithmsLib.SdssShapeControl)
class SdssShapeConfig(pexConf.Config): pass

@pexConf.register("APERTURE", registries.photometry)
@pexConf.wrap(algorithmsLib.AperturePhotometryControl)
class AperturePhotometryConfig(pexConf.Config): pass

@pexConf.register("GAUSSIAN", registries.photometry)
@pexConf.wrap(algorithmsLib.GaussianPhotometryControl)
class GaussianPhotometryConfig(pexConf.Config): pass

@pexConf.register("NAIVE", registries.photometry)
@pexConf.wrap(algorithmsLib.NaivePhotometryControl)
class NaivePhotometryConfig(pexConf.Config): pass

@pexConf.register("PSF", registries.photometry)
@pexConf.wrap(algorithmsLib.PsfPhotometryControl)
class PsfPhotometryConfig(pexConf.Config): pass

@pexConf.register("SINC", registries.photometry)
@pexConf.wrap(algorithmsLib.SincPhotometryControl)
class SincPhotometryConfig(pexConf.Config):
    def _get_radius(self): return self.radius2
    def _set_radius(self, r): self.radius2 = r
    radius = property(_get_radius, _set_radius, doc="synonym for radius2")
