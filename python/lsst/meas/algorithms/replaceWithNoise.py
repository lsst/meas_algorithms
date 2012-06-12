import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.detection as afwDet
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage

class ReplaceWithNoiseConfig(pexConfig.Config):
    """
    Configuration for ReplaceWithNoiseTask.
    """
    noiseSource = pexConfig.ChoiceField(doc='How do we choose the mean and variance of the Gaussian noise we generate?',
                                      dtype=str, allowed={
                                          'measure': 'Measure clipped mean and variance from the whole image',
                                          'meta': 'Mean = 0, variance = the "BGMEAN" metadata entry',
                                          'variance': "Mean = 0, variance = the image's variance",
                                          },
                                      default='measure',
                                      optional=False)

    noiseOffset = pexConfig.Field(dtype=float, optional=False, default=0.,
                                  doc='Add ann offset to the generated noise.')

    noiseSeed = pexConfig.Field(dtype=int, default=0, doc='The seed value to use for random number generation.')

class ReplaceWithNoiseTask(pipeBase.Task):
    ConfigClass = ReplaceWithNoiseConfig
    _DefaultName = "replaceWithNoise"

    def begin(self, exposure, sources, noiseImage=None, noiseMeanVar=None):
        # creates heavies, replaces all footprints with noise
        # We need the source table to be sorted by ID to do the parent lookups
        # (sources.find() below)
        if not sources.isSorted():
            sources.sort()
        mi = exposure.getMaskedImage()
        im = mi.getImage()
        mask = mi.getMask()

        # Add temporary Mask planes for THISDET and OTHERDET
        self.removeplanes = []
        bitmasks = []
        for maskname in ['THISDET', 'OTHERDET']:
            try:
                # does it already exist?
                plane = mask.getMaskPlane(maskname)
                self.log.logdebug('Mask plane "%s" already existed' % maskname)
            except:
                # if not, add it; we should delete it when done.
                plane = mask.addMaskPlane(maskname)
                self.removeplanes.append(maskname)
            mask.clearMaskPlane(plane)
            bitmask = mask.getPlaneBitMask(maskname)
            bitmasks.append(bitmask)
            self.log.logdebug('Mask plane "%s": plane %i, bitmask %i = 0x%x' %
                              (maskname, plane, bitmask, bitmask))
        self.thisbitmask,self.otherbitmask = bitmasks
        del bitmasks

        # Start by creating HeavyFootprints for each source.
        #
        # The "getParent()" checks are here because top-level
        # sources (ie, those with no parents) are not supposed to
        # have HeavyFootprints, but child sources (ie, those that
        # have been deblended) should have HeavyFootprints
        # already.
        self.heavies = []
        for source in sources:
            fp = source.getFootprint()
            if source.getParent():
                # this source has been deblended; "fp" should
                # already be a HeavyFootprint.
                # Swig downcasts it to Footprint, so we have to re-cast.
                self.heavies.append(afwDet.cast_HeavyFootprintF(fp))
            else:
                # top-level source: copy pixels from the input
                # image.
                ### FIXME: the heavy footprint includes the mask
                ### and variance planes, which we shouldn't need
                ### (I don't think we ever want to modify them in
                ### the input image).  Copying them around is
                ### wasteful.
                heavy = afwDet.makeHeavyFootprint(fp, mi)
                self.heavies.append(heavy)

        # We now create a noise HeavyFootprint for each top-level Source.
        # We'll put the noisy footprints in a map from id -> HeavyFootprint:
        self.heavyNoise = {}
        noisegen = self.getNoiseGenerator(exposure, noiseImage, noiseMeanVar)
        self.log.logdebug('Using noise generator: %s' % (str(noisegen)))
        for source in sources:
            if source.getParent():
                continue
            fp = source.getFootprint()
            heavy = noisegen.getHeavyFootprint(fp)
            self.heavyNoise[source.getId()] = heavy
            # Also insert the noisy footprint into the image now.
            # Notice that we're just inserting it into "im", ie,
            # the Image, not the MaskedImage.
            heavy.insert(im)
            # Also set the OTHERDET bit
            afwDet.setMaskFromFootprint(mask, fp, self.otherbitmask)

    def insertSource(self, exposure, sourcei):
        # Copy this source's pixels into the image
        mi = exposure.getMaskedImage()
        im = mi.getImage()
        mask = mi.getMask()
        fp = self.heavies[sourcei]
        fp.insert(im)
        afwDet.setMaskFromFootprint(mask, fp, self.thisbitmask)
        afwDet.clearMaskFromFootprint(mask, fp, self.otherbitmask)

    def removeSource(self, exposure, sources, source):
        # remove a single source
        # (Replace this source's pixels by noise again.)
        # Do this by finding the source's top-level ancestor
        mi = exposure.getMaskedImage()
        im = mi.getImage()
        mask = mi.getMask()
        ancestor = source
        j = 0
        while ancestor.getParent():
            ancestor = sources.find(ancestor.getParent())
            j += 1
            if not ancestor or j == 100:
                raise RuntimeError('Source hierarchy too deep, or (more likely) your Source table is botched.')
        # Re-insert the noise pixels
        fp = self.heavyNoise[ancestor.getId()]
        fp.insert(im)
        # Clear the THISDET mask plane.
        afwDet.clearMaskFromFootprint(mask, fp, self.thisbitmask)
        afwDet.setMaskFromFootprint(mask, fp, self.otherbitmask)

    def end(self, exposure, sources):
        # restores original image, cleans up temporaries
        # (ie, replace all the top-level pixels)
        mi = exposure.getMaskedImage()
        im = mi.getImage()
        mask = mi.getMask()
        for source,heavy in zip(sources,self.heavies):
            if source.getParent():
                continue
            heavy.insert(im)
        for maskname in self.removeplanes:
            mask.removeAndClearMaskPlane(maskname, True)

        del self.removeplanes
        del self.thisbitmask
        del self.otherbitmask
        del self.heavies
        del self.heavyNoise

    def getNoiseGenerator(self, exposure, noiseImage, noiseMeanVar):
        if noiseImage is not None:
            return ImageNoiseGenerator(noiseImage)
        rand = None
        if self.config.noiseSeed:
            # default algorithm, our seed
            rand = afwMath.Random(afwMath.Random.MT19937, self.config.noiseSeed)
        if noiseMeanVar is not None:
            try:
                # Assume noiseMeanVar is an iterable of floats
                noiseMean,noiseVar = noiseMeanVar
                noiseMean = float(noiseMean)
                noiseVar = float(noiseVar)
                noiseStd = math.sqrt(noiseVar)
                self.log.logdebug('Using passed-in noise mean = %g, variance = %g -> stdev %g' %
                                  (noiseMean, noiseVar, noiseStd))
                return FixedGaussianNoiseGenerator(noiseMean, noiseStd, rand=rand)
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
                # We would have to adjust for GAIN if ip_isr didn't make it 1.0
                noiseStd = math.sqrt(bgMean)
                self.log.logdebug('Using noise variance = (BGMEAN = %g) from exposure metadata' %
                                  (bgMean))
                return FixedGaussianNoiseGenerator(offset, noiseStd, rand=rand)
            except:
                self.log.logdebug('Failed to get BGMEAN from exposure metadata')

        if noiseSource == 'variance':
            self.log.logdebug('Will draw noise according to the variance plane.')
            var = exposure.getMaskedImage().getVariance()
            return VariancePlaneNoiseGenerator(var, mean=offset, rand=rand)

        # Compute an image-wide clipped variance.
        im = exposure.getMaskedImage().getImage()
        s = afwMath.makeStatistics(im, afwMath.MEANCLIP | afwMath.STDEVCLIP)
        noiseMean = s.getValue(afwMath.MEANCLIP)
        noiseStd = s.getValue(afwMath.STDEVCLIP)
        self.log.logdebug("Measured from image: clipped mean = %g, stdev = %g" %
                          (noiseMean,noiseStd))
        return FixedGaussianNoiseGenerator(noiseMean + offset, noiseStd, rand=rand)

class NoiseGenerator(object):
    '''
    Base class for noise generators used by the "doReplaceWithNoise" routine:
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

