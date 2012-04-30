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

    noiseSource = pexConf.ChoiceField(doc='If "doRemoveOtherSources" is set, how do choose the mean and variance of the Gaussian noise we generate?',
                                      dtype=str, allowed={
                                          'measure': 'Measure clipped mean and variance from the whole image',
                                          'meta': 'Mean = 0, variance = the "BGMEAN" metadata entry',
                                          'variance': "Mean = 0, variance = the image's variance",
                                          },
                                      default='measure',
                                      optional=False)

    noiseOffset = pexConf.Field(dtype=float, optional=False, default=0.,
                                doc='If "doRemoveOtherSources" is set, add this value to the noise.')

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
    def run(self, exposure, sources, apCorr=None, noiseImage=None,
            noiseMeanVar=None):
        """Run measure() and applyApCorr().

        @param[in]     exposure Exposure to process
        @param[in,out] sources  SourceCatalog containing sources detected on this exposure.
        @param[in]     apCorr   ApertureCorrection object to apply.
        @param[in]     noiseImage   (passed to measure(); see there for documentation)
        @param[in]     noiseMeanVar (passed to measure(); see there for documentation)

        @return None

        The aperture correction is only applied if config.doApplyApCorr is True and the apCorr
        argument is not None.
        """
        self.measure(exposure, sources, noiseImage=noiseImage, noiseMeanVar=noiseMeanVar)
        if self.config.doApplyApCorr and apCorr:
            self.applyApCorr(sources, apCorr)

    def preMeasureHook(self, exposure, sources):
        '''A hook, for debugging purposes, that is called at the start of the
        measure() method.'''
        try:
            import lsstDebug
            self._smt_display = lsstDebug.Info(__name__).display
        except ImportError, e:
            try:
                self._smt_display = display
            except NameError:
                self._smt_display = False
        if self._smt_display:
            frame = 0
            ds9.mtv(exposure, title="input", frame=frame)
            ds9.cmdBuffer.pushSize()

    def postMeasureHook(self, exposure, sources):
        '''A hook, for debugging purposes, that is called at the end of the
        measure() method.'''
        if hasattr(self, '_smt_display') and self._smt_display:
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
        if hasattr(self, '_smt_display') and self._smt_display:
            if self._smt_display > 1:
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
                noiseImage=None, noiseMeanVar=None):
        """Measure sources on an exposure, with no aperture correction.

        @param[in]     exposure Exposure to process
        @param[in,out] sources  SourceCatalog containing sources detected on this exposure.
        @param[in]     noiseImage If 'config.doRemoveOtherSources = True', you can pass in
                       an Image containing noise.  This overrides the "config.noiseSource" setting.
        @param[in]     noiseMeanVar: if 'config.doRemoveOtherSources = True', you can specify
                       the mean and variance of the Gaussian noise that will be added, by passing
                       a tuple of (mean, variance) floats.  This overrides the "config.noiseSource"
                       setting (but is overridden by noiseImage).
                       
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

        print 'self.__module__:', self.__module__
        print '__name__:', __name__
        print exposure.getMetadata().toString()

        noiseout = self.config.doRemoveOtherSources
        if noiseout:
            # We need the source table to be sorted by ID to do the parent lookups
            # (sources.find() below)
            if not sources.isSorted():
                sources.sort()
            mi = exposure.getMaskedImage()
            im = mi.getImage()
            mask = mi.getMask()

            # Add Mask planes for THISDET and OTHERDET
            removeplanes = []
            bitmasks = []
            for maskname in ['THISDET', 'OTHERDET']:
                try:
                    # does it already exist?
                    plane = mask.getMaskPlane(maskname)
                    self.log.logdebug('Mask plane "%s" already existed' % maskname)
                except:
                    # if not, add it; we should delete it when done.
                    plane = mask.addMaskPlane(maskname)
                    removeplanes.append(maskname)
                mask.clearMaskPlane(plane)
                bitmask = mask.getPlaneBitMask(maskname)
                bitmasks.append(bitmask)
                self.log.logdebug('Mask plane "%s": plane %i, bitmask %i = 0x%x' %
                                  (maskname, plane, bitmask, bitmask))
            thisbitmask,otherbitmask = bitmasks
            del bitmasks

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
                    # Swig downcasts it to Footprint, so we have to re-cast.
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
            # We'll put the noisy footprints in a map from id -> HeavyFootprint:
            heavyNoise = {}
            noisegen = self.getNoiseGenerator(exposure, noiseImage, noiseMeanVar)
            self.log.logdebug('Using noise generator: %s' % (str(noisegen)))
            for source in sources:
                if source.getParent():
                    continue
                fp = source.getFootprint()
                heavy = noisegen.getHeavyFootprint(fp)
                heavyNoise[source.getId()] = heavy
                # Also insert the noisy footprint into the image now.
                # Notice that we're just inserting it into "im", ie,
                # the Image, not the MaskedImage.
                heavy.insert(im)
                # Also set the OTHERDET bit
                afwDet.setMaskFromFootprint(mask, fp, otherbitmask)
            # At this point the whole image should just look like noise.

        # Call the hook before we measure anything...
        self.preSingleMeasureHook(exposure, sources, -1)

        for i,source in enumerate(sources):

            if noiseout:
                # Copy this source's pixels into the image
                fp = heavies[i]
                fp.insert(im)
                afwDet.setMaskFromFootprint(mask, fp, thisbitmask)
                afwDet.clearMaskFromFootprint(mask, fp, otherbitmask)

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
                fp = heavyNoise[ancestor.getId()]
                fp.insert(im)
                # Clear the THISDET mask plane.
                afwDet.clearMaskFromFootprint(mask, fp, thisbitmask)
                afwDet.setMaskFromFootprint(mask, fp, otherbitmask)

        if noiseout:
            # Put the exposure back the way it was (ie, replace all the top-level pixels)
            for source,heavy in zip(sources,heavies):
                if source.getParent():
                    continue
                heavy.insert(im)

            for maskname in removeplanes:
                mask.removeAndClearMaskPlane(maskname, True)

        self.postMeasureHook(exposure, sources)

    def getNoiseGenerator(self, exposure, noiseImage, noiseMeanVar):
        if noiseImage is not None:
            return ImageNoiseGenerator(noiseImage)

        if noiseMeanVar is not None:
            try:
                # Assume noiseMeanVar is an iterable of floats
                noiseMean,noiseVar = noiseMeanVar
                noiseMean = float(noiseMean)
                noiseVar = float(noiseVar)
                noiseStd = math.sqrt(noiseVar)
                self.log.logdebug('Using passed-in noise mean = %g, variance = %g -> stdev %g' %
                                  (noiseMean, noiseVar, noiseStd))
                return FixedGaussianNoiseGenerator(noiseMean, noiseStd)
            except:
                self.log.logdebug('Failed to cast passed-in noiseMeanVar to floats: %s' %
                                  (str(noiseMeanVar)))

        offset = self.config.noiseOffset
        noiseSource = self.config.noiseSource

        if noiseSource == 'meta':
            # check the exposure metadata
            meta = exposure.getMetadata()
            # this key name correspond to estimateBackground() in detection.py
            try:
                bgMean = meta.getAsDouble('BGMEAN')
                # Do we need to correct for gain?  Probably not if our ip_isr has run, since it
                # renormalizes the CCD to have gain = 1.
                noiseVar = bgMean
                noiseStd = math.sqrt(noiseVar)
                self.log.logdebug('Using noise variance = (BGMEAN = %g) from exposure metadata' %
                                  (bgMean))
                return FixedGaussianNoiseGenerator(offset, noiseStd)
            except:
                self.log.logdebug('Failed to get BGMEAN from exposure metadata')

        if noiseSource == 'variance':
            self.log.logdebug('Will draw noise according to the variance plane.')
            var = exposure.getMaskedImage().getVariance()
            return VariancePlaneNoiseGenerator(var, mean=offset)

        # Compute an image-wide clipped variance.
        im = exposure.getMaskedImage().getImage()
        s = afwMath.makeStatistics(im, afwMath.MEANCLIP | afwMath.STDEVCLIP)
        noiseMean = s.getValue(afwMath.MEANCLIP)
        noiseStd = s.getValue(afwMath.STDEVCLIP)
        self.log.logdebug("Measured from image: clipped mean = %g, stdev = %g" %
                          (noiseMean,noiseStd))
        return FixedGaussianNoiseGenerator(noiseMean + offset, noiseStd)

            
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


class NoiseGenerator(object):
    '''
    Base class for noise generators used by the "doRemoveOtherSources" routine:
    these produce HeavyFootprints filled with noise generated in various ways.

    This is an abstract base class.
    '''
    def getHeavyFootprint(self, fp):
        bb = fp.getBBox()
        mim = self.getMaskedImage(bb)
        return afwDet.makeHeavyFootprint(fp, mim)
    def getMaskedImage(self, bb):
        im = self.getImage(bb)
        return afwImage.MaskedImageF(im)
    def getImage(self, bb):
        return None

class ImageNoiseGenerator(NoiseGenerator):
    '''
    "Generates" noise by cutting out a subimage from a user-supplied noise Image.
    '''
    def __init__(self, img):
        '''
        img: an afwImage.ImageF
        '''
        self.mim = afwImage.MaskedImageF(img)
    def getMaskedImage(self, bb):
        return self.mim

class GaussianNoiseGenerator(NoiseGenerator):
    '''
    Generates noise using the afwMath.Random() and afwMath.randomGaussianImage() routines.

    This is an abstract base class.
    '''
    def __init__(self, rand=None):
        if rand is None:
            rand = afwMath.Random()
        self.rand = rand
    def getRandomImage(self, bb):
        # Create an Image and fill it with Gaussian noise.
        rim = afwImage.ImageF(bb.getWidth(), bb.getHeight())
        rim.setXY0(bb.getMinX(), bb.getMinY())
        afwMath.randomGaussianImage(rim, self.rand)
        return rim

class FixedGaussianNoiseGenerator(GaussianNoiseGenerator):
    '''
    Generates Gaussian noise with a fixed mean and standard deviation.
    '''
    def __init__(self, mean, std, rand=None):
        super(FixedGaussianNoiseGenerator, self).__init__(rand=rand)
        self.mean = mean
        self.std = std
    def __str__(self):
        return 'FixedGaussianNoiseGenerator: mean=%g, std=%g' % (self.mean, self.std)
    def getImage(self, bb):
        rim = self.getRandomImage(bb)
        rim *= self.std
        rim += self.mean
        return rim

class VariancePlaneNoiseGenerator(GaussianNoiseGenerator):
    '''
    Generates Gaussian noise whose variance matches that of the variance plane of the image.
    '''
    def __init__(self, var, mean=None, rand=None):
        '''
        var: an afwImage.ImageF; the variance plane.
        mean: floating-point or afwImage.Image
        '''
        super(VariancePlaneNoiseGenerator, self).__init__(rand=rand)
        self.var = var
        if mean is not None and mean == 0.:
            mean = None
        self.mean = mean
    def __str__(self):
        return 'VariancePlaneNoiseGenerator: mean=' + str(self.mean)
    def getImage(self, bb):
        rim = self.getRandomImage(bb)
        # Use the image's variance plane to scale the noise.
        stdev = afwImage.ImageF(self.var, bb, afwImage.LOCAL, True)
        stdev.sqrt()
        rim *= stdev
        if self.mean is not None:
            rim += self.mean
        return rim

