import lsst.afw.math        as math
import lsst.afw.table       as afwTable
import lsst.afw.image       as afwImg
import lsst.meas.algorithms as measAlg
import lsst.afw.display.ds9 as ds9

class GaussianFluxTest():
    """This class reads sensor data from PhoSim and processes it."""

    def readFile(self):
        """Read the file and make the exposure etc"""

        self.exposure    = afwImg.ExposureF('../../meas_base/tests/calexp_one.fits')
        self.truth_flux = [60000.0,]
        # self.exposure    = afwImg.ExposureF('../..//meas_base/tests/calexp.fits')
        # self.truth_flux = [100000.0,75000.0,1E15]

        self.maskedImage = self.exposure.getMaskedImage()

    def doStats(self):
        """Do pixel based processing"""

        statFlags = math.NPOINT | math.MEAN | math.STDEV | math.MAX | math.MIN | math.ERRORS
        control = math.StatisticsControl()

        imageStatistics = math.makeStatistics(self.maskedImage, statFlags, control)
        numBins         = imageStatistics.getResult(math.NPOINT)[0]
        mean            = imageStatistics.getResult(math.MEAN)[0]

        print "The image has dimensions %i x %i pixels" %(self.maskedImage.getWidth(),
                                                          self.maskedImage.getHeight())
        print "Number of analyzed bins in image is %i"  %numBins
        print "Max    = %9d"                            %imageStatistics.getResult(math.MAX)[0]
        print "Min    = %9d"                            %imageStatistics.getResult(math.MIN)[0]
        print "Mean   = %9.8f +- %3.1f"                 %imageStatistics.getResult(math.MEAN)
        print "StdDev = %9.2f"                          %imageStatistics.getResult(math.STDEV)[0]

    def detectAndMeasureSources(self,fixed,display):
        """Find and measure the sources"""

        # Configure the detection and measurement algorithms
        schema                = afwTable.SourceTable.makeMinimalSchema()
        schema.setVersion(0)

        detectSourcesConfig   = measAlg.SourceDetectionConfig()
        detectSourcesConfig.thresholdValue = 5.0
        measureSourcesConfig  = measAlg.SourceMeasurementConfig()
        measureSourcesConfig.algorithms["flux.gaussian"].fixed = fixed

        # Setup the detection and measurement tasks
        detect  = measAlg.SourceDetectionTask  (config=detectSourcesConfig,  schema=schema)
        measure = measAlg.SourceMeasurementTask(config=measureSourcesConfig, schema=schema)

        # Detect the sources,then put them into a catalog (the table is where the catalog atually stores stuff)
        table   = afwTable.SourceTable.make(schema)
        catalog = detect.makeSourceCatalog(table, self.exposure, sigma=5)

        # Get the sources out of the catalog and apply Measurement routines
        sources = catalog.sources
        measure.run(self.exposure, sources)

        # Now let's look at the output from some of the measurment algorithms.
        fields = ['centroid.sdss', 'shape.sdss', 'flux.gaussian']
        keys   = [schema.find(f).key for f in fields]

        if True :
            i = -1
            for source in sources:
                i = i+1
                print "\tSource ", source.get('id'), ' ==>  fixed = ', fixed
                for f,k in zip(fields, keys):
                    print "\t",f, source.get(k)
                gaussianFlux = source.get('flux.gaussian')
                print "   % diff = ", 200.0*(gaussianFlux-self.truth_flux[i])/(gaussianFlux+self.truth_flux[i])
                print " rel diff = ", gaussianFlux/self.truth_flux[i]

        if display:
            frame = 1
            ds9.mtv(self.exposure, frame=frame)
            ds9.setMaskTransparency(50)
            with ds9.Buffering():
                for s in sources:
                    xy = s.getCentroid()
                    ds9.dot('+', *xy, ctype=ds9.CYAN if s.get("flags.negative") else ds9.GREEN, frame=frame)
                    ds9.dot(s.getShape(), *xy, ctype=ds9.RED, frame=frame)


    def run(self,fixed,display):
        self.readFile()
        # self.doStats()
        self.detectAndMeasureSources(fixed,display)

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description="Testing GaussianFlux")
    parser.add_argument('--fixed', action="store_true", help="Fix centroid and shape to slot values", default=False)
    parser.add_argument('--ds9', action="store_true", help="Display sources on ds9", default=False)
    args = parser.parse_args()

    GaussianFluxTest().run(fixed=args.fixed,display=args.ds9)
