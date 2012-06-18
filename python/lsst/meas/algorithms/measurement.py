# 
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2011 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import math
import numpy

import lsst.pex.config as pexConfig
import lsst.afw.table as afwTable
import lsst.pipe.base as pipeBase
import lsst.afw.display.ds9 as ds9

from . import algorithmsLib
from .algorithmRegistry import *
from .replaceWithNoise import *

__all__ = "SourceSlotConfig", "SourceMeasurementConfig", "SourceMeasurementTask"

class SourceSlotConfig(pexConfig.Config):

    centroid = pexConfig.Field(dtype=str, default="centroid.sdss", optional=True,
                             doc="the name of the centroiding algorithm used to set source x,y")
    shape = pexConfig.Field(dtype=str, default="shape.sdss", optional=True,
                          doc="the name of the algorithm used to set source moments parameters")
    apFlux = pexConfig.Field(dtype=str, default="flux.sinc", optional=True,
                           doc="the name of the algorithm used to set the source aperture flux slot")
    modelFlux = pexConfig.Field(dtype=str, default="flux.gaussian", optional=True,
                           doc="the name of the algorithm used to set the source model flux slot")
    psfFlux = pexConfig.Field(dtype=str, default="flux.psf", optional=True,
                            doc="the name of the algorithm used to set the source psf flux slot")
    instFlux = pexConfig.Field(dtype=str, default="flux.gaussian", optional=True,
                             doc="the name of the algorithm used to set the source inst flux slot")

    def setupTable(self, table, prefix=None):
        """Convenience method to setup a table's slots according to the config definition.

        This is defined in the Config class to support use in unit tests without needing
        to construct a Task object.
        """
        if prefix is None: prefix = ""
        if self.centroid is not None: table.defineCentroid(prefix + self.centroid)
        if self.shape is not None: table.defineShape(prefix + self.shape)
        if self.apFlux is not None: table.defineApFlux(prefix + self.apFlux)
        if self.modelFlux is not None: table.defineModelFlux(prefix + self.modelFlux)
        if self.psfFlux is not None: table.definePsfFlux(prefix + self.psfFlux)
        if self.instFlux is not None: table.defineInstFlux(prefix + self.instFlux)

class SourceMeasurementConfig(pexConfig.Config):
    """
    Configuration for SourceMeasurementTask.
    A configured instance of MeasureSources can be created using the
    makeMeasureSources method.
    """

    slots = pexConfig.ConfigField(
        dtype = SourceSlotConfig,
        doc="Mapping from algorithms to special aliases in Source.\n"
        )

    algorithms = AlgorithmRegistry.all.makeField(
        multi=True,
        default=["flags.pixel",
                 "centroid.gaussian", "centroid.naive",
                 "shape.sdss",
                 "flux.gaussian", "flux.naive", "flux.psf", "flux.sinc",
                 "classification.extendedness",
                 "skycoord",
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

    apCorrFluxes = pexConfig.ListField(
        dtype=str, optional=False, default=["flux.psf", "flux.gaussian"],
        doc="Fields to which we should apply the aperture correction.  Elements in this list"\
            "are silently ignored if they are not in the algorithms list, to make it unnecessary"\
            "to always keep them in sync."
        )
    doApplyApCorr = pexConfig.Field(dtype=bool, default=True, optional=False, doc="Apply aperture correction?")

    # We might want to make this default to True once we have battle-tested it
    # Formerly known as "doRemoveOtherSources"
    doReplaceWithNoise = pexConfig.Field(dtype=bool, default=False, optional=False,
                                       doc='When measuring, replace other detected footprints with noise?')

    replaceWithNoise = pexConfig.ConfigurableField(
        target = ReplaceWithNoiseTask,
        doc = ("Task for replacing other sources by noise when measuring sources; run when " +
               "'doReplaceWithNoise' is set."),
    )

    prefix = pexConfig.Field(dtype=str, optional=True, default=None, doc="prefix for all measurement fields")

    def setDefaults(self):
        self.slots.centroid = self.centroider.name
        self.slots.shape = "shape.sdss"
        self.slots.psfFlux = "flux.psf"
        self.slots.apFlux = "flux.naive"
        self.slots.modelFlux = "flux.gaussian"
        self.slots.instFlux = "flux.gaussian"

    def validate(self):
        pexConfig.Config.validate(self)
        if self.centroider.name in self.algorithms.names:
            raise ValueError("The algorithm in the 'centroider' field must not also appear in the "\
                                 "'algorithms' field.")
        if self.slots.centroid is not None and (self.slots.centroid not in self.algorithms.names
                                                and self.slots.centroid != self.centroider.name):
            raise ValueError("source centroid slot algorithm '%s' is not being run." % self.slots.astrom)
        if self.slots.shape is not None and self.slots.shape not in self.algorithms.names:
            raise ValueError("source shape slot algorithm '%s' is not being run." % self.slots.shape)
        for slot in (self.slots.psfFlux, self.slots.apFlux, self.slots.modelFlux, self.slots.instFlux):
            if slot is not None:
                for name in self.algorithms.names:
                    if len(name) <= len(slot) and name == slot[:len(name)]:
                        break
                else:
                    raise ValueError("source flux slot algorithm '%s' is not being run." % slot)
        if self.doApplyApCorr:
            for flux in self.apCorrFluxes:
                for name in self.algorithms.names:
                    if len(name) <= len(flux) and name == flux[:len(name)]:
                        break
                else:
                    raise ValueError(
                        "flux algorithm '%s' is to be aperture-corrected but is not being run." % flux)
                

    def makeMeasureSources(self, schema, metadata=None):
        """ Convenience method to make a MeasureSources instance and
        fill it with the configured algorithms.

        This is defined in the Config class to support use in unit tests without needing
        to construct a Task object.
        """
        builder = algorithmsLib.MeasureSourcesBuilder(self.prefix if self.prefix is not None else "")
        if self.centroider is not None:
            builder.setCentroider(self.centroider.apply())
        builder.addAlgorithms(self.algorithms.apply())
        return builder.build(schema, metadata)

class SourceMeasurementTask(pipeBase.Task):
    """Measure the properties of sources on a single exposure.

    This task has no return value; it only modifies the SourceCatalog in-place.
    """
    ConfigClass = SourceMeasurementConfig
    _DefaultName = "sourceMeasurement"

    def __init__(self, schema, algMetadata=None, **kwds):
        """Create the task, adding necessary fields to the given schema.

        @param[in,out] schema        Schema object for measurement fields; will be modified in-place.
        @param[in,out] algMetadata   Passed to MeasureSources object to be filled with initialization
                                     metadata by algorithms (e.g. radii for aperture photometry).
        @param         **kwds        Passed to Task.__init__.
        """
        pipeBase.Task.__init__(self, **kwds)
        self.measurer = self.config.makeMeasureSources(schema, algMetadata)
        if self.config.doApplyApCorr:
            self.fluxKeys = [(schema.find(f).key, schema.find(f + ".err").key)
                             for f in self.config.apCorrFluxes]
            self.corrKey = schema.addField("aperturecorrection", type=float,
                                           doc="aperture correction factor applied to fluxes")
            self.corrErrKey = schema.addField("aperturecorrection.err", type=float,
                                              doc="aperture correction uncertainty")
        else:
            self.corrKey = None
            self.corrErrKey = None
        if self.config.doReplaceWithNoise:
            self.makeSubtask('replaceWithNoise')


    @pipeBase.timeMethod
    def run(self, exposure, sources, apCorr=None, noiseImage=None,
            noiseMeanVar=None, references=None, refWcs=None):
        """Run measure() and applyApCorr().

        @param[in]     exposure Exposure to process
        @param[in,out] sources  SourceCatalog containing sources detected on this exposure.
        @param[in]     apCorr   ApertureCorrection object to apply.
        @param[in]     noiseImage   (passed to measure(); see there for documentation)
        @param[in]     noiseMeanVar (passed to measure(); see there for documentation)
        @param[in]     references SourceCatalog containing reference sources detected on reference exposure.
        @param[in]     refWcs     Wcs for the reference exposure.
        @return None

        The aperture correction is only applied if config.doApplyApCorr is True and the apCorr
        argument is not None.
        """
        self.measure(exposure, sources, noiseImage=noiseImage, noiseMeanVar=noiseMeanVar,
                     references=references, refWcs=refWcs)
        if self.config.doApplyApCorr and apCorr:
            self.applyApCorr(sources, apCorr)

    def preMeasureHook(self, exposure, sources):
        '''A hook, for debugging purposes, that is called at the start of the
        measure() method.'''

        # pipe_base's Task provides self._display.
        if self._display:
            frame = 0
            ds9.mtv(exposure, title="input", frame=frame)
            ds9.cmdBuffer.pushSize()

    def postMeasureHook(self, exposure, sources):
        '''A hook, for debugging purposes, that is called at the end of the
        measure() method.'''
        if self._display:
            ds9.cmdBuffer.popSize()

    def preSingleMeasureHook(self, exposure, sources, i):
        '''A hook, for debugging purposes, that is called immediately before
        the measurement algorithms for each source.

        Note that this will also be called with i=-1 just before entering the
        loop over measuring sources.'''
        pass

    def postSingleMeasureHook(self, exposure, sources, i):
        '''A hook, for debugging purposes, that is called immediately after
        the measurement algorithms.'''
        self.postSingleMeasurementDisplay(exposure, sources[i])

    def postSingleMeasurementDisplay(self, exposure, source):
        if self._display:
            if self._display > 1:
                ds9.dot(str(source.getId()), source.getX() + 2, source.getY(),
                        size=3, ctype=ds9.RED)
                cov = source.getCentroidErr()
                ds9.dot(("@:%.1f,%.1f,%1f" % (cov[0,0], cov[0,1], cov[0,0])),
                        source.getX(), source.getY(), size=3, ctype=ds9.RED)
                symb = "%d" % source.getId()
            else:
                symb = "+"
                ds9.dot(symb, source.getX(), source.getY(), size=3, ctype=ds9.RED)
            print source.getX(), source.getY(), source.getPsfFlux(), source.getModelFlux()
    
    @pipeBase.timeMethod
    def measure(self, exposure, sources, noiseImage=None, noiseMeanVar=None, references=None, refWcs=None):
        """Measure sources on an exposure, with no aperture correction.

        @param[in]     exposure Exposure to process
        @param[in,out] sources  SourceCatalog containing sources detected on this exposure.
        @param[in]     noiseImage If 'config.doReplaceWithNoise = True', you can pass in
                       an Image containing noise.  This overrides the "config.noiseSource" setting.
        @param[in]     noiseMeanVar: if 'config.doReplaceWithNoise = True', you can specify
                       the mean and variance of the Gaussian noise that will be added, by passing
                       a tuple of (mean, variance) floats.  This overrides the "config.noiseSource"
                       setting (but is overridden by noiseImage).
        @param[in]     references SourceCatalog containing reference sources detected on reference exposure.
        @param[in]     refWcs     Wcs for the reference exposure.
        @return None
        """
        if references is None:
            references = [None] * len(sources)
        if len(sources) != len(references):
            raise RuntimeError("Number of sources (%d) and references (%d) don't match" %
                               (len(sources), len(references)))

        if self.config.doReplaceWithNoise and not hasattr(self, 'replaceWithNoise'):
            self.makeSubtask('replaceWithNoise')

        self.log.info("Measuring %d sources" % len(sources))
        self.config.slots.setupTable(sources.table, prefix=self.config.prefix)

        self.preMeasureHook(exposure, sources)

        # "noiseout": we will replace all the pixels within detected
        # Footprints with noise, and then add sources in one at a
        # time, measure them, then replace with noise again.  The idea
        # is that measurement algorithms might look outside the
        # Footprint, and we don't want other sources to interfere with
        # the measurements.  The faint wings of sources are still
        # there, but that's life.
        noiseout = self.config.doReplaceWithNoise
        if noiseout:
            self.replaceWithNoise.begin(exposure, sources, noiseImage, noiseMeanVar)
            # At this point the whole image should just look like noise.

        # Call the hook, with source id = -1, before we measure anything.
        # (this is *after* the sources have been replaced by noise, if noiseout)
        self.preSingleMeasureHook(exposure, sources, -1)

        for i, (source, ref) in enumerate(zip(sources, references)):
            if noiseout:
                self.replaceWithNoise.insertSource(exposure, i)

            self.preSingleMeasureHook(exposure, sources, i)

            # Make the measurement
            if ref is None:
                self.measurer.apply(source, exposure)
            else:
                self.measurer.apply(source, exposure, ref, refWcs)

            self.postSingleMeasureHook(exposure, sources, i)

            if noiseout:
                # Replace this source's pixels by noise again.
                self.replaceWithNoise.removeSource(exposure, sources, source)

        if noiseout:
            # Put the exposure back the way it was
            self.replaceWithNoise.end(exposure, sources)

        self.postMeasureHook(exposure, sources)
            
    @pipeBase.timeMethod
    def applyApCorr(self, sources, apCorr):
        import numpy
        self.log.log(self.log.INFO, "Applying aperture correction to %d sources" % len(sources))
        for source in sources:
            corr, corrErr = apCorr.computeAt(source.getX(), source.getY())
            for fluxKey, fluxErrKey in self.fluxKeys:
                flux = source.get(fluxKey)
                fluxErr = source.get(fluxErrKey)
                source.set(fluxKey, flux * corr)
                source.set(fluxErrKey, (fluxErr**2 * corr**2 + flux**2 * corrErr**2)**0.5)
            source.set(self.corrKey, corr)
            source.set(self.corrErrKey, corrErr)

