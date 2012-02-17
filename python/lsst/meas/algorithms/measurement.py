import collections

import lsst.pex.config as pexConf
from . import algorithmsLib

@pexConf.wrap(algorithmsLib.AlgorithmControl)
class AlgorithmConfig(pexConf.Config):
    order = 0     # used to sort algorithms; 0 means no dependence on other algorithms.

@pexConf.wrap(algorithmsLib.CentroidControl)
class CentroidConfig(AlgorithmConfig):
    pass

@pexConf.wrap(algorithmsLib.ShapeControl)
class ShapeConfig(AlgorithmConfig):
    order = 100   # shape algorithms usually want to be run after centroids; subclasses can override.

@pexConf.wrap(algorithmsLib.FluxControl)
class FluxConfig(AlgorithmConfig):
    order = 200   # flux algorithms usually want to be run after shapes; subclasses can override.

class AlgorithmRegistry(pexConf.Registry):
    """A customized registry for source measurement algorithms.

    Using a customizated registry allows us to avoid a lot of the
    boilerplate that would otherwise be necessary when implementing a
    new source measurement algorithm.

    First, this registry class has the ability to views and associated
    Fields that can only refer to algorithms with a particular
    intermediate base class, all while referencing the same underlying
    registry.  This allows some fields to be restricted to a
    particular subclass while others are not:

    class MyConfig(Config):
        any = AlgorithmRegistry.all.makeField(
            "field that accepts any algorithm"
        )
        centroids = AlgorithmRegistry.filter(CentroidConfig).makeField(
            "only allows centroiders"
        )

    The only instance of this registry is the class attribute 'all'
    (hence 'AlgorithmRegistry.all', above).  The 'filter' class method
    is used to create filtered views into it.  Note that the base
    class is a config class, not a control class.

    Second, because all algorithms must have swigged C++ control
    classes, and these are essentially duplicates of the Config
    classes users interact with in a RegistryField, the Config classes
    themselves can be created transparently from the control classes
    when the control class is registered.  This is only done if the
    control class does not already have a ConfigClass attribute and no
    explicit ConfigClass argument is passed to the register method.

    Third, the registry provides a customized instance dict for its
    associated registry fields, which ensures the list of active
    algorithms is always sorted according to their config's 'order'
    class member.

    The 'configurable' held in this registry is a callable that
    returns a control object instance, and takes no additional
    arguments.  This simply calls config.makeControl() and sets the
    control object's name data member (the config does not have a name
    field, since it would be confusing to allow it to differ from the
    registration name).

    All config classes registered must have a makeControl() method
    that returns a control instance (this is usually provided by using
    pex.config.makeConfigClass on the swigged control class).  """

    class SubclassRegistryView(collections.Mapping):
        """A read-only view into registry that filters out items whose
        config classes don't inherit from a particular base class."""

        def __init__(self, registry, base):
            self.registry = registry
            self.base = base

        def __getitem__(self, k):
            t = self.registry[k]
            if not issubclass(t.ConfigClass, self.base):
                raise TypeError("%r is not a subclass of %r" % (t.ConfigClass, self.base))
            return t

        def __len__(self):
            return len(tuple(iter(self)))

        def __iter__(self):
            return (k for k, v in self.registry.iteritems() if issubclass(v.ConfigClass, self.base))

        def __contains__(self, k):
            v = self.registry.get(k)
            return v is not None and issubclass(v.ConfigClass, self.base)

        def makeField(self, doc, default=None, optional=False, multi=False):
            return pexConf.RegistryField(doc, self, default, optional, multi,
                                         instanceDictClass=AlgorithmRegistry.InstanceDict)

    class Configurable(object):
        """Class used as the actual element in the registry; a
        callable that returns a swigged C++ control object with its
        name set from the registry when called.  """

        __slots__ = "ConfigClass", "name"

        def __init__(self, name, ConfigClass):
            self.name = name
            self.ConfigClass = ConfigClass

        def __call__(self, config):
            ctrl = config.makeControl()
            ctrl.name = self.name
            return ctrl

    def __init__(self):
        """Construct a registry of AlgorithmControl classes."""
        pexConf.Registry.__init__(self, AlgorithmConfig)

    class InstanceDict(pexConf.RegistryInstanceDict):
        """A custom instance dict for associated registry fields that
        keeps configs sorted.  """

        def _setSelection(self, value):
            pexConf.RegistryInstanceDict._setSelection(self, value)
            if self._multi and self._selection is not None:
                self._selection.sort(key=lambda name: self[name].order)

    @classmethod
    def register(cls, name, target, ConfigClass=None, order=0):
        """Register an AlgorithmControl subclass.
        
        This is a class method, so you can either use it on the
        registry instance or its class (this works because the
        registry is a singleton, so the class knows to use cls.all as
        the instance).

        If it does not have a ConfigClass attribute pointing to the
        corresponding Config class, a config class will be created
        using pex.config.makeConfigClass.  A new config class will
        also be created if the ConfigClass attribute was inherited
        from a base class,

        @param[in] name         Name the algorithm will be registered
                                with; also the name of the source fields
                                it will fill.
        @param[in] target       An AlgorithmControl subclass.
        @param[in] ConfigClass  A Config class to be paired with the
                                control class.
        @param[in] order        Sets the sort order for the algorithm if
                                the ConfigClass is created; ignored
                                otherwise.
        """
        self = cls.all
        if not issubclass(target, algorithmsLib.AlgorithmControl):
            raise TypeError("Registry targets must be subclasses of AlgorithmControl")
        if ConfigClass is None:
            if hasattr(target, "ConfigClass") and (not hasattr(target.__base__, "ConfigClass")
                                                   or target.ConfigClass != target.__base__.ConfigClass):
                ConfigClass = target.ConfigClass     # class attributes
            else:
                if not hasattr(target.__base__, "ConfigClass"):
                    raise ValueError("Cannot create a config class for %s unless its base class "
                                     "has a ConfigClass attribute." % target)
                ConfigClass = pexConf.makeConfigClass(target, base=target.__base__.ConfigClass, module=2)
        target = self.Configurable(name, ConfigClass)
        pexConf.Registry.register(self, name, target)

    def duplicate(self, oldName, newName):
        """Register an existing algorithm class with a new name."""
        old = self[oldName]
        target = self.Configurable(newName, old.ConfigClass)
        pexConf.Registry.register(self, newName, target)

    def makeField(self, doc, default=None, optional=False, multi=False):
        return pexConf.RegistryField(doc, self, default, optional, multi,
                                     instanceDictClass=AlgorithmRegistry.InstanceDict)

    @classmethod
    def filter(cls, base):
        """Return a lazy read-only view that only contains items with
        the given Config (not Control) base class.
        """
        return cls.SubclassRegistryView(cls.all, base)

AlgorithmRegistry.all = AlgorithmRegistry()

class SourceSlotConfig(pexConf.Config):

    centroid = pexConf.Field(dtype=str, default="centroid.sdss", optional=True,
                             doc="the name of the centroiding algorithm used to set source x,y")
    shape = pexConf.Field(dtype=str, default="shape.sdss", optional=True,
                          doc="the name of the algorithm used to set source moments parameters")
    apFlux = pexConf.Field(dtype=str, default="flux.sinc", optional=True,
                           doc="the name of the algorithm used to set the source aperture flux slot")
    modelFlux = pexConf.Field(dtype=str, default="flux.gaussian", optional=True,
                           doc="the name of the algorithm used to set the source model flux slot")
    psfFlux = pexConf.Field(dtype=str, default="flux.psf", optional=True,
                           doc="the name of the algorithm used to set the source psf flux slot")
    instFlux = pexConf.Field(dtype=str, default="flux.gaussian", optional=True,
                             doc="the name of the algorithm used to set the source inst flux slot")

    def setupTable(self, table):
        if self.centroid is not None: table.defineCentroid(self.centroid)
        if self.shape is not None: table.defineShape(self.shape)
        if self.apFlux is not None: table.defineApFlux(self.apFlux)
        if self.modelFlux is not None: table.defineModelFlux(self.modelFlux)
        if self.psfFlux is not None: table.definePsfFlux(self.psfFlux)
        if self.instFlux is not None: table.defineInstFlux(self.instFlux)

class MeasureSourcesConfig(pexConf.Config):
    """
    Configuration for the MeasureSources C++ class.

    A configured instance of MeasureSources can be created using the
    makeMeasureSources method.
    """

    slots = pexConf.ConfigField(
        dtype = SourceSlotConfig,
        doc="Mapping from algorithms to special aliases in Source.\n"\
            "WARNING: the MeasureSources object does not setup the\n"\
            "slots automatically (because the slot definition is never\n"\
            "passed to C++).  This should be done by calling \n"\
            "slots.setupTable() after the table is created."
        )

    algorithms = AlgorithmRegistry.all.makeField(
        multi=True,
        default=["flags.pixel",
                 "centroid.gaussian", "centroid.naive",
                 "shape.sdss",
                 "flux.gaussian", "flux.naive", "flux.psf", "flux.sinc",
                 "classification.extendedness",
                 ],
        doc="Configuration and selection of measurement algorithms."
        )
    
    centroider = AlgorithmRegistry.filter(CentroidConfig).makeField(
        multi=False, default="centroid.sdss", optional=True,
        doc="Configuration for the initial centroid algorithm used to\n"\
            "feed center points to other algorithms.\n\n"\
            "Note that this is in addition to the centroider listed in\n"\
            "the 'algorithms' field; the same name should not appear in\n"\
            "both.\n\n"\
            "This field DOES NOT set which field name will be used to define\n"\
            "the alias for source.getX(), source.getY(), etc.\n"
        )

    def validate(self):
        pexConf.Config.validate(self)
        if self.centroider.name in self.algorithms.names:
            raise ValueError("The algorithm in the 'centroider' field must not also appear in the "\
                                 "'algorithms' field.")
        if self.slots.centroid is not None and (self.slots.centroid not in self.algorithms.names
                                                and self.slots.centroid != self.centroider.name):
            raise ValueError("source centroid slot algorithm '%s' is not being run." % self.slots.astrom)
        if self.slots.shape is not None and self.slots.shape not in self.algorithms.names:
            raise ValueError("source shape slot algorithm '%s' is not being run." % self.slots.shape)
        for slot in (self.slots.psfFlux, self.slots.apFlux, self.slots.modelFlux, self.slots.instFlux):
            if slot is not None and slot not in self.algorithms.names:
                raise ValueError("source flux slot algorithm '%s' is not being run." % slot)

    def makeMeasureSources(self, schema=None):
        """ Convenience method to make a MeasureSources instance and
        fill it with the configured algorithms.
        """
        if schema is None:
            ms = algorithmsLib.MeasureSources()
        else:
            ms = algorithmsLib.MeasureSources(schema)
        if self.centroider is not None:
            ms.setCentroider(self.centroider.apply())
        ms.addAlgorithms(self.algorithms.apply())
        return ms

AlgorithmRegistry.register("classification.extendedness", algorithmsLib.ClassificationControl, order=300)
AlgorithmRegistry.register("flags.pixel", algorithmsLib.PixelFlagControl)
AlgorithmRegistry.register("centroid.gaussian", algorithmsLib.GaussianCentroidControl)
AlgorithmRegistry.register("centroid.naive", algorithmsLib.NaiveCentroidControl)
AlgorithmRegistry.register("centroid.sdss", algorithmsLib.SdssCentroidControl)
AlgorithmRegistry.register("shape.sdss", algorithmsLib.SdssShapeControl)
AlgorithmRegistry.register("flux.aperture", algorithmsLib.ApertureFluxControl, order=200)
AlgorithmRegistry.register("flux.gaussian", algorithmsLib.GaussianFluxControl)
AlgorithmRegistry.register("flux.naive", algorithmsLib.NaiveFluxControl)
AlgorithmRegistry.register("flux.psf", algorithmsLib.PsfFluxControl)

# Here's an example on how to declare a measurement config more manually, and add a property to the Config.
@pexConf.wrap(algorithmsLib.SincFluxControl)
class SincFluxConfig(FluxConfig):
    def _get_radius(self): return self.radius2
    def _set_radius(self, r): self.radius2 = r
    radius = property(_get_radius, _set_radius, doc="synonym for radius2")
AlgorithmRegistry.register("flux.sinc", target=algorithmsLib.SincFluxControl, ConfigClass=SincFluxConfig)
