from math import sqrt, log, exp, pi
import os, sys, unittest
import eups
import lsst.utils.tests as utilsTests
import lsst.afw.detection as afwDet
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.meas.algorithms as measAlg

import lsst.meas.extensions.shapeHSM

class RemoveOtherSourcesTestCase(unittest.TestCase):

    def test1(self):
        FWHM = 5
        ksize = 25
        psf = afwDet.createPsf("DoubleGaussian", ksize, ksize,
                               FWHM/(2*sqrt(2*log(2))), 1, 0.1)

        im = afwImage.ImageF(100, 100)

        sky = 100
        rand = afwMath.Random()
        afwMath.randomGaussianImage(im, rand)
        im *= sky

        xc,yc = 40,50
        sx,sy = 10,3
        W,H = im.getWidth(), im.getHeight()
        for y in xrange(H):
            for x in xrange(W):
                g = exp(-0.5 * ((x - xc)**2/(float(sx)**2) + (y - yc)**2/(float(sy)**2)))
                g /= (2.*pi*sx*sy)
                im.set(x,y, im.get(x,y) + g * 1e6)

        # shape.sdss (ixx=99.2656584815, iyy=25.0290814687, ixy=-0.3852295282)

        for x,y,flux in [(80,50,1000)]: #, (50,50,1000), (55,50,1000)]:
            psfim = psf.computeImage(afwGeom.Point2D(x,y)).convertF()
            psfim *= float(flux)
            bbox = psfim.getBBox(afwImage.PARENT)
            imview = afwImage.ImageF(im, bbox)
            imview += psfim

        im.writeFits('im.fits')

        mi = afwImage.MaskedImageF(im)
        exposure = afwImage.makeExposure(mi)
        exposure.setPsf(psf)

        detconf = measAlg.SourceDetectionConfig()
        detconf.reEstimateBackground = False

        measconf = measAlg.SourceMeasurementConfig()
        measconf.doApplyApCorr = False

        newalgs = [ 'shape.hsm.ksb', 'shape.hsm.bj', 'shape.hsm.linear' ]
        measconf.algorithms = list(measconf.algorithms.names) + newalgs
        #measconf.algorithms += newalgs

        schema = afwTable.SourceTable.makeMinimalSchema()
        detect = measAlg.SourceDetectionTask(config=detconf, schema=schema)
        measure = measAlg.SourceMeasurementTask(config=measconf, schema=schema)

        print 'Running detection...'
        table = afwTable.SourceTable.make(schema)
        detected = detect.makeSourceCatalog(table, exposure)
        sources = detected.sources

        print 'Running measurement...'
        measure.run(exposure, sources)

        #fields = schema.getNames()
        #print 'Fields:', fields
        fields = ['centroid.sdss', 'shape.sdss',
                  'shape.hsm.bj.moments',
                  'shape.hsm.ksb.moments',
                  'shape.hsm.linear.moments',
                  ]
        keys = [schema.find(f).key for f in fields]
        #print 'Fields:', fields
        #print 'Keys:', keys
        for source in sources:
            print '  ', source
            for f,k in zip(fields, keys):
                val = source.get(k)
                print '    ', f, val
        

#-=-=-=-=-=-=-=-=-=-=-=-=silly boilerplate line-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(RemoveOtherSourcesTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)
def run(exit = False):
    """Run the utilsTests"""
    utilsTests.run(suite(), exit)
if __name__ == "__main__":
    run(True)


            
