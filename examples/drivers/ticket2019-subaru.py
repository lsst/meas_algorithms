debugPlots = True
if debugPlots:
    import matplotlib
    matplotlib.use('Agg')


import lsst.pipe.base as pipeBase
import lsst.pipe.tasks.processCcd as procCcd
import lsst.daf.persistence as dafPersist
import lsst.obs.suprimecam as obsSc
import lsst.pex.logging as pexLog

import lsst.daf.base as dafBase
import lsst.afw.table as afwTable
from lsst.ip.isr import IsrTask
from lsst.pipe.tasks.calibrate import CalibrateTask

if debugPlots:
    #
    # Plot the image that will be passed to the measurement algorithms,
    # for a few bright sources in a small region.
    #
    import pylab as plt
    import numpy as np
    import lsst.meas.algorithms as measAlg
    import lsst.afw.math as afwMath
    class DebugSourceMeasTask(measAlg.SourceMeasurementTask):
        def __str__(self):
            return 'DebugSourceMeasTask'
        def run(self, *args, **kwargs):
            print 'DebugSourceMeasTask running.'
            super(DebugSourceMeasTask,self).run(*args, **kwargs)
        def _plotimage(self, im):
            xlo,xhi,ylo,yhi = self.plotregion
            plt.clf()
            if not isinstance(im, np.ndarray):
                im = im.getArray()
            plt.imshow(im[ylo:yhi, xlo:xhi],
                       extent=[xlo,xhi,ylo,yhi], **self.plotargs)
            plt.gray()
            
        def preMeasureHook(self, exposure, sources):
            mi = exposure.getMaskedImage()
            im = mi.getImage()
            s = afwMath.makeStatistics(im, afwMath.STDEVCLIP + afwMath.MEANCLIP)
            mn = s.getValue(afwMath.MEANCLIP)
            std = s.getValue(afwMath.STDEVCLIP)
            lo = mn -  3*std
            hi = mn + 20*std
            self.plotargs = dict(interpolation='nearest', origin='lower',
                                 vmin=lo, vmax=hi)
            self.plotregion = (100, 500, 100, 500)
            self.nplots = 0
            self._plotimage(im)
            plt.savefig('pre.png')
        def postMeasureHook(self, exposure, sources):
            mi = exposure.getMaskedImage()
            im = mi.getImage()
            self._plotimage(im)
            plt.savefig('post.png')
        def preSingleMeasureHook(self, exposure, sources, i):
            if i != -1:
                return
            mi = exposure.getMaskedImage()
            mask = mi.getMask()
            print 'Mask planes:'
            mask.printMaskPlanes()
            bitmask = mask.getPlaneBitMask('THISDET')
            print 'THISDET bitmask:', bitmask
            ma = mask.getArray()
            
            oldargs = self.plotargs
            args = oldargs.copy()
            args.update(vmin=0, vmax=1)
            self.plotargs = args

            mim = ((ma & bitmask) > 0)
            self._plotimage(mim)
            plt.savefig('mask-thisdet.png')

            for i in range(mask.getNumPlanesUsed()):
                bitmask = (1 << i)
                mim = ((ma & bitmask) > 0)
                self._plotimage(mim)
                plt.savefig('mask-bit%02i.png' % i)

            self.plotargs = oldargs

        def postSingleMeasureHook(self, exposure, sources, i):
            xlo,xhi,ylo,yhi = self.plotregion
            x,y = sources[i].getX(), sources[i].getY()
            if x < xlo or x > xhi or y < ylo or y > yhi:
                return
            if (not np.isfinite(x)) or (not np.isfinite(y)):
                return
            if sources[i].getPsfFlux() < 20000.:
                return
            print 'Source at', x,y
            print 'PSF flux', sources[i].getPsfFlux()
            #print 'model flux', sources[i].getModelFlux()
            #print 'ap flux', sources[i].getApFlux()
            #print 'inst flux', sources[i].getInstFlux()
            if self.nplots >= 5:
                return
            mi = exposure.getMaskedImage()
            im = mi.getImage()
            self._plotimage(im)
            plt.savefig('meas%02i.png' % self.nplots)

            mask = mi.getMask()
            bitmask = mask.getPlaneBitMask('THISDET')
            print 'THISDET bitmask:', bitmask
            mim = ((mask.getArray() & bitmask) > 0)
            oldargs = self.plotargs
            args = oldargs.copy()
            args.update(vmin=0, vmax=1)
            self.plotargs = args
            self._plotimage(mim)
            self.plotargs = oldargs
            plt.savefig('meas%02i-mask.png' % self.nplots)

            self.nplots += 1


if __name__ == '__main__':
    mapperArgs = dict(root='/lsst/home/dstn/lsst/ACT-data/rerun/dstn',
                      calibRoot='/lsst/home/dstn/lsst/ACT-data/CALIB')
    mapper = obsSc.SuprimecamMapper(**mapperArgs)
    butlerFactory = dafPersist.ButlerFactory(mapper = mapper)
    butler = butlerFactory.create()
    print 'Butler', butler
    dataRef = butler.subset('raw', dataId = dict(visit=126969, ccd=5))
    print 'dataRef:', dataRef
    #dataRef.butlerSubset = dataRef
    #print 'dataRef:', dataRef
    print 'len(dataRef):', len(dataRef)
    for dr in dataRef:
        print '  ', dr

            
    conf = procCcd.ProcessCcdConfig()
    conf.measurement.doRemoveOtherSources = True

    conf.doDetection = True
    conf.doMeasurement = True
    conf.doWriteSources = True
    conf.doWriteCalibrate = True

    conf.calibrate.doComputeApCorr = False
    conf.calibrate.doAstrometry = False
    conf.calibrate.doPhotoCal = False
    # the default Cosmic Ray parameters don't work for Subaru images
    cr = conf.calibrate.repair.cosmicray
    cr.minSigma = 10.
    cr.min_DN = 500.
    cr.niteration = 3
    cr.nCrPixelMax = 1000000

    if debugPlots:
        conf.doMeasurement = False

    proc = procCcd.ProcessCcdTask(config=conf, name='ProcessCcd')
    proc.log.setThreshold(pexLog.Log.DEBUG)
    if debugPlots:
        proc.measurement = DebugSourceMeasTask(proc.schema,
                                               algMetadata=proc.algMetadata,
                                               config=conf.measurement)
        conf.doMeasurement = True

    proc.measurement.log.setThreshold(pexLog.Log.DEBUG)

    conf.calibrate.measurement.doApplyApCorr = False
    conf.measurement.doApplyApCorr = False
    conf.validate()

    for dr in dataRef:
        print 'dr', dr
        # Only do ISR and Calibration if necessary...
        doIsr   = False
        doCalib = False
        print 'calexp', mapper.map('calexp', dr.dataId)
        print 'psf', mapper.map('psf', dr.dataId)
        try:
            psf = dr.get('psf')
            print 'PSF:', psf
        except:
            print 'No PSF'
            doCalib = True
        try:
            expo = dr.get('calexp')
            print 'Calexp:', expo
        except:
            print 'No calexp'
            doCalib = True

        if doCalib:
            # calib needs 'postISRCCD'
            try:
                pisr = dr.get('postISRCCD')
                print 'postISRCCD:', pisr
            except:
                print 'No postISRCCD'
                doIsr = True
                conf.doWriteIsr = True

        conf.doIsr = doIsr
        conf.doCalibrate = doCalib

        print 'proc.run'
        proc.run(dr)

        if False:
            conf.doMeasurement = False
            res = proc.run(dr)
            proc.measurement.run(res.exposure, res.sources, noiseMeanVar='variance')
            proc.measurement.run(res.exposure, res.sources, noiseMeanVar='measure')
            conf.doMeasurement = True

        
