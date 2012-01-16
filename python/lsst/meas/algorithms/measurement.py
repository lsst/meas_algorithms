from lsst.pex.config import Config, Field, ListField, ConfigField, RegistryField, register

class SourceConfig(Config):
    astrom = Field("The name of the centroiding algorithm used to set Source.[XY]Astrom", str)
    shape = Field("The name of the centroiding algorithm used to set Source.Mxx etc.", str)
    apFlux = Field("The name of the algorithm used to set Source.apFlux(Err)", str)
    modelFlux = Field("The name of the algorithm used to set Source.modelFlux(Err)", str)
    psfFlux = Field("The name of the algorithm used to set Source.psfFlux(Err)", str)
    instFlux = Field("The name of the algorithm used to set Source.instFlux(Err)", str)

class ClassificationConfig(Config):
    sg_fac1 = Field("First S/G parameter; critical ratio of inst to psf flux", dtype=float, 
                    default=0.925, optional=True)
    sg_fac2 = Field("Second S/G parameter; correction for instFlux error", dtype=float,
                    default=0.0, optional=True)
    sg_fac3 = Field("Third S/G parameter; correction for psfFlux error", dtype=float,
                    default=0, optional=True)

class MeasureSourcesConfig(Config):

    source = ConfigField("The mapping from algorithms to fields in Source", SourceConfig)

    astrometry = RegistryField("Configurations for individual astrometry algorithms.", multi=True)
    astrometry.defaults = ["GAUSSIAN", "NAIVE", "SDSS"]

    shape = RegistryField("Configurations for various shape-measurement algorithms.", multi=True)
    shape.defaults = ["SDSS"]
        
    photometry = RegistryField("Configurations for individual photometry algorithms.", multi=True)
    photometry.defaults = ["GAUSSIAN", "NAIVE", "PSF", "SINC"]

    classification = ConfigField("Parameters to do with star/galaxy classification", ClassificationConfig)

    def __init__(self, *args, **kwds):
        Config.__init__(self, *args, **kwds)
        self.astrometry.names = type(self).astrometry.defaults
        self.shape.names = type(self).shape.defaults
        self.photometry.names = type(self).photometry.defaults

@register("GAUSSIAN", MeasureSourcesConfig.astrometry)
class GaussianAstrometryConfig(Config):
    pass

@register("NAIVE", MeasureSourcesConfig.astrometry)
class NaiveAstrometryConfig(Config):
    background = Field("FIXME! NEVER DOCUMENTED!", dtype=float, optional=True)

@register("SDSS", MeasureSourcesConfig.astrometry)
class SdssAstrometryConfig(Config):
    binmax = Field("FIXME! NEVER DOCUMENTED!", dtype=int, optional=True)
    peakMin = Field("FIXME! NEVER DOCUMENTED!", dtype=float, optional=True)
    wfac = Field("FIXME! NEVER DOCUMENTED!", dtype=float, optional=True)

@register("SDSS", MeasureSourcesConfig.shape)
class SdssShapeConfig(Config):
    background = Field("FIXME! NEVER DOCUMENTED!", dtype=float, optional=True)

@register("APERTURE", MeasureSourcesConfig.photometry)
class AperturePhotometryConfig(Config):
    radius = ListField("Radii for apertures (in pixels)", itemCheck=(lambda x: x > 0.0), dtype=float)

@register("GAUSSIAN", MeasureSourcesConfig.photometry)
class GaussianPhotometryConfig(Config):
    fixed = Field("FIXME! NEVER DOCUMENTED!", dtype=bool, default=False, optional=True)
    background = Field("FIXME! NEVER DOCUMENTED!", dtype=float, default=0.0, optional=True)
    shiftmax = Field("FIXME! NEVER DOCUMENTED!", dtype=float, default=10.0, optional=True)

@register("NAIVE", MeasureSourcesConfig.photometry)
class NaivePhotometryConfig(Config):
    radius = Field("FIXME! NEVER DOCUMENTED!", dtype=float, default=9.0, optional=True)

@register("PSF", MeasureSourcesConfig.photometry)
class PsfPhotometryConfig(Config):
    pass

@register("SINC", MeasureSourcesConfig.photometry)
class SincPhotometryConfig(Config):
    radius2 = Field("FIXME! NEVER DOCUMENTED!", dtype=float)
    radius1 = Field("FIXME! NEVER DOCUMENTED!", dtype=float, optional=True)
    angle = Field("FIXME! NEVER DOCUMENTED!", dtype=float, optional=True)
    ellipticity = Field("FIXME! NEVER DOCUMENTED!", dtype=float, optional=True) 

    def _get_radius(self): return self.radius2
    def _set_radius(self, r): self.radius2 = r
    radius = property(_get_radius, _set_radius, doc="synonym for radius2")
