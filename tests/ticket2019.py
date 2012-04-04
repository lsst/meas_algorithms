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

plots = True
if plots:
    try:
        import matplotlib
        matplotlib.use('Agg')
        import pylab as plt
        from matplotlib.patches import Ellipse
    except:
        plots = False

# Add an axis-aligned Gaussian to an image, in a horribly inefficient way.
def addGaussian(im, xc, yc, sx, sy, flux):
    W,H = im.getWidth(), im.getHeight()
    for y in xrange(H):
        for x in xrange(W):
            g = exp(-0.5 * ((x - xc)**2/(float(sx)**2) + (y - yc)**2/(float(sy)**2)))
            g /= (2.*pi*sx*sy)
            im.set(x,y, im.get(x,y) + g * flux)

def addPsf(im, psf, xc, yc, flux):
    psfim = psf.computeImage(afwGeom.Point2D(xc,yc)).convertF()
    psfim *= float(flux)
    bbox = psfim.getBBox(afwImage.PARENT)
    # clip in case the PSF goes outside the image
    bbox.clip(im.getBBox())
    psfim = afwImage.ImageF(psfim, bbox, afwImage.PARENT)
    imview = afwImage.ImageF(im, bbox)
    imview += psfim


def plotSources(im, sources, schema):
    plt.clf()
    plt.imshow(im.getArray(), origin='lower', interpolation='nearest',
               vmin=-100, vmax=500)
    plt.gray()
    shapekey = schema.find('shape.sdss').key
    xykey = schema.find('centroid.sdss').key
    for source in sources:
        quad = source.get(shapekey)
        ixx,iyy = quad.getIxx(), quad.getIyy()
        x,y = source.get(xykey)
        sx,sy = sqrt(ixx),sqrt(iyy)
        plt.gca().add_artist(Ellipse([x,y], 2.*sx, 2.*sy, angle=0.,
                                     ec='r', fc='none', lw=2, alpha=0.5))
    

class RemoveOtherSourcesTestCase(unittest.TestCase):

    def test1(self):
        FWHM = 5
        ksize = 25
        psf = afwDet.createPsf("DoubleGaussian", ksize, ksize,
                               FWHM/(2*sqrt(2*log(2))), 1, 0.1)

        im = afwImage.ImageF(120, 200)

        sky = 100
        rand = afwMath.Random()
        afwMath.randomGaussianImage(im, rand)
        im *= sky

        x = 60
        #for dx,y in [(30, 20), (35, 40), (40, 60), (45, 80)]:
        for i in range(5):
            dx = 28 + i*4
            #dx = 24 + i*4
            y = 20 + i*40
            #                x   y  sx sy flux
            #addGaussian(im, x, y, 10, 3, 1.5e5)
            addGaussian(im, x, y, 10, 3, 1e5)
            addPsf(im, psf, x+dx, y, 1000)
            addPsf(im, psf, x-dx, y, 1000)

        im.writeFits('im.fits')

        mi = afwImage.MaskedImageF(im)
        var = mi.getVariance()
        print 'Variance:', var
        var = sky
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
        

        if plots:
            plotSources(im, sources, schema)
            plt.savefig('1a.png')

        measconf.doRemoveOtherSources = True
        print 'Running measurement...'
        measure.run(exposure, sources)

        fields = schema.getNames()
        print 'Fields:', fields
        fields = ['centroid.sdss', 'shape.sdss',
                  'shape.hsm.bj.moments',
                  'shape.hsm.ksb.moments',
                  'shape.hsm.linear.moments',
                  #'flags.pixel',
                  ]
        keys = [schema.find(f).key for f in fields]
        #print 'Fields:', fields
        #print 'Keys:', keys
        for source in sources:
            print '  ', source
            for f,k in zip(fields, keys):
                val = source.get(k)
                print '    ', f, val

        if plots:
            plotSources(im, sources, schema)
            plt.savefig('1b.png')



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


            
