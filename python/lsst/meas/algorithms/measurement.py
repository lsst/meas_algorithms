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
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDet

from . import algorithmsLib
from .algorithmRegistry import *

__all__ = "SourceSlotConfig", "SourceMeasurementConfig", "SourceMeasurementTask"

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

class SourceMeasurementConfig(pexConf.Config):
    """
    Configuration for SourceMeasurementTask.
    A configured instance of MeasureSources can be created using the
    makeMeasureSources method.
    """

    slots = pexConf.ConfigField(
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

    apCorrFluxes = pexConf.ListField(
        dtype=str, optional=False, default=["flux.psf", "flux.gaussian"],
        doc="Fields to which we should apply the aperture correction.  Elements in this list"\
            "are silently ignored if they are not in the algorithms list, to make it unnecessary"\
            "to always keep them in sync."
        )
    doApplyApCorr = pexConf.Field(dtype=bool, default=True, optional=False, doc="Apply aperture correction?")

    # We might want to make this default to True once we have battle-tested it
    doRemoveOtherSources = pexConf.Field(dtype=bool, default=False, optional=False,
                                         doc='When measuring, replace other detected footprints with noise?')

    prefix = pexConf.Field(dtype=str, optional=True, default=None, doc="prefix for all measurement fields")

    def setDefaults(self):
        self.slots.centroid = self.centroider.name
        self.slots.shape = "shape.sdss"
        self.slots.psfFlux = "flux.psf"
        self.slots.apFlux = "flux.naive"
        self.slots.modelFlux = "flux.gaussian"
        self.slots.instFlux = "flux.gaussian"

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
                             for f in self.config.apCorrFluxes if f in self.config.algorithms.names]
            self.corrKey = schema.addField("aperturecorrection", type=float,
                                           doc="aperture correction factor applied to fluxes")
            self.corrErrKey = schema.addField("aperturecorrection.err", type=float,
                                              doc="aperture correction uncertainty")
        else:
            self.corrKey = None
            self.corrErrKey = None


    @pipeBase.timeMethod
    def run(self, exposure, sources, apCorr=None, noiseImage=None):
        """Run measure() and applyApCorr().

        @param[in]     exposure Exposure to process
        @param[in,out] sources  SourceCatalog containing sources detected on this exposure.
        @param[in]     apCorr   ApertureCorrection object to apply.
        @param[in]     noiseImage  (passed to measure(); see there for documentation)

        @return None

        The aperture correction is only applied if config.doApplyApCorr is True and the apCorr
        argument is not None.
        """
        self.measure(exposure, sources, noiseImage=noiseImage)
        if self.config.doApplyApCorr and apCorr:
            self.applyApCorr(sources, apCorr)

    def preMeasureHook(self, exposure, sources):
        '''A hook, for debugging purposes, that is called at the start of the
        measure() method.'''
        try:
            import lsstDebug
            self.display = lsstDebug.Info(__name__).display
        except ImportError, e:
            try:
                self.display = display
            except NameError:
                self.display = False
        if self.display:
            frame = 0
            ds9.mtv(exposure, title="input", frame=frame)
            ds9.cmdBuffer.pushSize()

    def postMeasureHook(self, exposure, sources):
        '''A hook, for debugging purposes, that is called at the end of the
        measure() method.'''
        if self.display:
            ds9.cmdBuffer.popSize()

    def preSingleMeasureHook(self, exposure, sources, i):
        '''A hook, for debugging purposes, that is called immediately before
        the measurement algorithms for each source'''
        pass

    def postSingleMeasureHook(self, exposure, sources, i):
        '''A hook, for debugging purposes, that is called immediately after
        the measurement algorithms.'''
        self.postSingleMeasurementDisplay(exposure, sources[i])

    def postSingleMeasurementDisplay(self, exposure, source):
        if self.display:
            if self.display > 1:
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
    def measure(self, exposure, sources,
                noiseImage=None):
        """Measure sources on an exposure, with no aperture correction.

        @param[in]     exposure Exposure to process
        @param[in,out] sources  SourceCatalog containing sources detected on this exposure.
        @param[in]     noiseImage If 'config.doRemoveOtherSources = True', you can pass in
                       an Image containing noise; if None, noise will be generated
                       automatically.
        @return None
        """
        self.preMeasureHook(exposure, sources)
                                                                                    
        self.log.info("Measuring %d sources" % len(sources))
        self.config.slots.setupTable(sources.table, prefix=self.config.prefix)

        # "noiseout": we will replace all the pixels within detected
        # Footprints with noise, and then add sources in one at a
        # time, measure them, then replace with noise again.  The idea
        # is that measurement algorithms might look outside the
        # Footprint, and we don't want other sources to interfere with
        # the measurements.  The faint wings of sources are still
        # there, but that's life.
        noiseout = self.config.doRemoveOtherSources
        if noiseout:
            # We need the source table to be sorted by ID to do the parent lookups
            # (sources.find() below)
            if not sources.isSorted():
                sources.sort()
            mi = exposure.getMaskedImage()
            im = mi.getImage()

            # Start by creating HeavyFootprints for each source.
            #
            # The "getParent()" checks are here because top-level
            # sources (ie, those with no parents) are not supposed to
            # have HeavyFootprints, but child sources (ie, those that
            # have been deblended) should have HeavyFootprints
            # already.
            heavies = []
            for source in sources:
                fp = source.getFootprint()
                if source.getParent():
                    # this source has been deblended; "fp" should
                    # already be a HeavyFootprint.
                    ### FIXME -- This cast shouldn't be necessary!!
                    #heavies.append(fp)
                    heavies.append(afwDet.cast_HeavyFootprintF(fp))

                else:
                    # top-level source: copy pixels from the input
                    # image.
                    ### FIXME: the heavy footprint includes the mask
                    ### and variance planes, which we shouldn't need
                    ### (I don't think we ever want to modify them in
                    ### the input image).  Copying them around is
                    ### wasteful.
                    heavy = afwDet.makeHeavyFootprint(fp, mi)
                    heavies.append(heavy)

            # We now create a noise HeavyFootprint for each top-level Source.
            if noiseImage is None:
                rand = afwMath.Random()
                # We compute an image-wide noise standard deviation.
                # We could instead scale each pixel by its variance.
                # This could be a config switch (or the user could
                # pass in an appropriate noise image via the
                # "noiseImage" parameter)
                s = afwMath.makeStatistics(mi.getVariance(), afwMath.MEDIAN)
                skystd = math.sqrt(s.getValue(afwMath.MEDIAN))
                self.log.logdebug("Measured median sky standard deviation: %g" % skystd)
            # We'll put the noisy footprints in a map from id -> HeavyFootprint:
            heavyNoise = {}
            for source in sources:
                if source.getParent():
                    continue
                fp = source.getFootprint()
                bb = fp.getBBox()
                if noiseImage is None:
                    # Create an Image and fill it with Gaussian noise.
                    rim = afwImage.ImageF(bb.getWidth(), bb.getHeight())
                    rim.setXY0(bb.getMinX(), bb.getMinY())
                    afwMath.randomGaussianImage(rim, rand)
                    rim *= skystd
                else:
                    # Use the given noiseImage.
                    rim = noiseImage
                # Pull the HeavyFootprint out of the random image.
                ### FIXME: As above, notice that here we have to
                ### create a MaskedImage with bogus mask and variance
                ### planes.
                heavy = afwDet.makeHeavyFootprint(fp, afwImage.MaskedImageF(rim))
                heavyNoise[source.getId()] = heavy
                # Also insert the noisy footprint into the image now.
                # Notice that we're just inserting it into "im", ie,
                # the Image, not the MaskedImage.
                heavy.insert(im)
            # At this point the whole image should just look like noise.

        for i,source in enumerate(sources):

            if noiseout:
                # Copy this source's pixels into the image
                heavies[i].insert(im)

            self.preSingleMeasureHook(exposure, sources, i)
            self.measurer.apply(source, exposure)
            self.postSingleMeasureHook(exposure, sources, i)

            if noiseout:
                # Replace this source's pixels by noise again.
                # Do this by finding the source's top-level ancestor
                ancestor = source
                j = 0
                while ancestor.getParent():
                    ancestor = sources.find(ancestor.getParent())
                    j += 1
                    if not ancestor or j == 100:
                        raise RuntimeError('Source hierarchy too deep, or (more likely) your Source table is botched.')
                # Re-insert the noise pixels
                heavyNoise[ancestor.getId()].insert(im)

        if noiseout:
            # Put the exposure back the way it was (ie, replace all the top-level pixels)
            for source,heavy in zip(sources,heavies):
                if source.getParent():
                    continue
                heavy.insert(im)

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

