import lsst.pipe.base as pipeBase
import lsst.pipe.tasks.processCcd as procCcd
#from lsst.meas.algorithms.SingleSourceMeasurementTask import *
import lsst.daf.persistence as dafPersist
import lsst.obs.suprimecam as obsSc

import lsst.daf.base as dafBase
import lsst.afw.table as afwTable
from lsst.ip.isr import IsrTask
from lsst.pipe.tasks.calibrate import CalibrateTask

class ProcessCcdTask(procCcd.ProcessCcdTask):
    _DefaultName = 'hackProcessCcdTask'
    def __init__(self, **kwargs):
        #### ARGH, I can't make this cooperate with ProcessCcdTask.
        # so that super.__init__ doesn't init the "measurement" subtask:
        # self.config.doMeasurement = False
        # super(ProcessCcdTask, self).__init__(**kwargs)
        # if self.config.doMeasurement:
        #     self.makeSubtask("measurement", SingleSourceMeasurementTask,
        #                      schema=self.schema, algMetadata=self.algMetadata)

        pipeBase.Task.__init__(self, **kwargs)
        self.makeSubtask("isr", IsrTask)
        self.makeSubtask("calibrate", CalibrateTask)
        self.schema = afwTable.SourceTable.makeMinimalSchema()
        self.algMetadata = dafBase.PropertyList()
        if self.config.doDetection:
            self.makeSubtask("detection", measAlg.SourceDetectionTask, schema=self.schema)
        if self.config.doMeasurement:
            self.makeSubtask("measurement", SingleSourceMeasurementTask,
                             schema=self.schema, algMetadata=self.algMetadata)

        #print 'self.config:'
        #printConfig(self.config)
        #print 'self.config.doIsr:', self.config.doIsr

    # def run(self, dataRef):
    #     print 'hack ProcessCcd.run():', dataRef
    #     print 'doIsr', self.config.doIsr
    # 
    #     super(procCcd.ProcessCcdTask.run,self).run(dataRef)


def printConfig(conf):
    for k,v in conf.items():
        print '  ', k, v

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

    #conf.measurement.doRemoveOtherSources = True
    conf.measurement.doRemoveOtherSources = False

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

    #proc = ProcessCcdTask(config=conf)
    proc = procCcd.ProcessCcdTask(config=conf, name='ProcessCcd')

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
        conf.doIsr = doIsr
        conf.doCalibrate = doCalib

        proc.run(dr)
    
