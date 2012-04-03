import math
import lsst.pipe.base as pipeBase
import lsst.meas.algorithms as measAlg
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDet

class SingleSourceMeasurementTask(measAlg.SourceMeasurementTask):
    '''
    A SourceMeasurementTask subclass that replaces all detected Footprints
    with noise before running the measurement algorithms.
    '''
    _DefaultName = "singleSourceMeasurement"

    @pipeBase.timeMethod
    def measure(self, exposure, sources):
        """Measure sources on an exposure, with no aperture correction.

        @param[in]     exposure Exposure to process
        @param[in,out] sources  SourceCatalog containing sources detected on this exposure.
        @return None
        """

        print 'Sources:', sources

        if not sources.isSorted():
            sources.sort()

        self.log.log(self.log.INFO, "Measuring %d sources" % len(sources))
        self.config.slots.setupTable(sources.table, prefix=self.config.prefix)
        mi = exposure.getMaskedImage()
        # Create HeavyFootprints
        heavies = []
        for source in sources:
            print 'source', source
            fp = source.getFootprint()
            print '  footprint', fp
            myid = source.getId()
            print '  id', myid
            parent = source.getParent()
            print '  parent', parent
            # I was going to just noise out the parent footprints, but since we
            # need to make HeavyFootprints anyway, we need them all.
            #if parent != 0:
            #    continue
            # What about deblended Sources that are already carrying their own
            # HeavyFootprints?
            heavy = afwDet.makeHeavyFootprint(fp, mi)
            heavies.append(heavy)

        # Noise out the parent footprints
        rand = afwMath.Random()
        im = mi.getImage()
        # Image-wide standard deviation
        #skystd = 100.
        s = afwMath.makeStatistics(mi.getVariance(), afwMath.MEDIAN)
        skystd = math.sqrt(s.getValue(afwMath.MEDIAN))
        print 'Median standard devation:', skystd
        # mapping from (parent) ID to random image.
        heavyNoise = {}
        for source in sources:
            parent = source.getParent()
            if parent != 0:
                continue
            fp = source.getFootprint()
            #for sp in fp.getSpans():
            bb = fp.getBBox()
            print 'BBox w,h', bb.getWidth(), bb.getHeight()
            rim = afwImage.ImageF(bb.getWidth(), bb.getHeight())
            rim.setXY0(bb.getMinX(), bb.getMinY())
            afwMath.randomGaussianImage(rim, rand)
            rim *= skystd
            print 'Random image:', rim

            # TWEAK HeavyFootprintCtrl to insert only the PIXELS ?
            
            heavy = afwDet.makeHeavyFootprint(fp, rim)
            heavyNoise[source.getId()] = heavy

        for source,heavy in zip(sources, heavies):
            # Insert this source's pixels...
            heavy.insert(mi)
            self.measurer.apply(source, exposure)

            # Re-noise out this source's pixel (via its top-level ancestor)
            ancestor = source
            while ancestor.getParent() != 0:
                ancestor = sources.find(ancestor.getParent())

            heavyNoise[ancestor.getId()].insert(mi)
            
