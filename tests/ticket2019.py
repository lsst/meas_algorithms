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
import lsst.pex.config as pexConfig
import numpy as np

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

    flagkeys = [schema.find(x).key for x in [
        'shape.sdss.flags.maxiter', 'shape.sdss.flags.shift',
        'shape.sdss.flags.unweighted', 'shape.sdss.flags.unweightedbad']]
    flagabbr = ['M','S','U','B']

    for source in sources:
        quad = source.get(shapekey)
        ixx,iyy = quad.getIxx(), quad.getIyy()
        x,y = source.get(xykey)
        sx,sy = sqrt(ixx),sqrt(iyy)
        plt.gca().add_artist(Ellipse([x,y], 2.*sx, 2.*sy, angle=0.,
                                     ec='r', fc='none', lw=2, alpha=0.5))
        fs = ''
        for j,key in enumerate(flagkeys):
            val = source.get(key)
            if val:
                fs += flagbbr[j]
        if len(fs):
            plt.text(x+5, y, fs, ha='left', va='center')
        

# We define a SourceMeasurementTask subclass to use the debug hooks.
class MySourceMeasurementTask(measAlg.SourceMeasurementTask):
    def __init__(self, *args, **kwargs):
        self.plotpat = kwargs.pop('plotpat', 'postmeas-%(sourcenum)i.png')
        self.doplot = kwargs.pop('doplot', True)
        super(MySourceMeasurementTask,self).__init__(*args, **kwargs)

    def postSingleMeasureHook(self, exposure, sources, i):
        if self.doplot:
            im = exposure.getMaskedImage().getImage()
            plotSources(im, [sources[i]], sources.getSchema())
            fn = self.plotpat % dict(sourcenum=i, sourceid=sources[i].getId())
            plt.savefig(fn)
            print 'Wrote', fn

class RemoveOtherSourcesTestCase(unittest.TestCase):
    def test2(self):
        # Check that doRemoveOtherSources works with deblended source
        # hierarchies.
        seed = 42
        rand = afwMath.Random(afwMath.Random.MT19937, seed)

        psf = self.getpsf()
        im = afwImage.ImageF(200, 50)
        skystd = 100
        afwMath.randomGaussianImage(im, rand)
        im *= skystd
        imorig = afwImage.ImageF(im, True)

        mi = afwImage.MaskedImageF(im)
        mi.getVariance().set(skystd**2)
        exposure = afwImage.makeExposure(mi)
        exposure.setPsf(psf)

        detconf = measAlg.SourceDetectionConfig()
        detconf.reEstimateBackground = False
        measconf = measAlg.SourceMeasurementConfig()
        measconf.doApplyApCorr = False
            
        schema = afwTable.SourceTable.makeMinimalSchema()
        detect = measAlg.SourceDetectionTask(config=detconf, schema=schema)
        measure = MySourceMeasurementTask(config=measconf, schema=schema)
        table = afwTable.SourceTable.make(schema)
        table.preallocate(10)

        # We're going to fake up a perfect deblend hierarchy here, by
        # creating individual images containing single sources and
        # measuring them, and then creating a deblend hierarchy where
        # the children have the correct HeavyFootprints.  We want to
        # find that the measurements on the deblend hierarchy and the
        # blended image are equal to the individual images.
        #
        # Note that we don't expect the measurements to be *identical*
        # because the measurement routine adds its own noise, so the pixels
        # within the parent footprint but outside the child footprint will
        # be different.
        
        fullim = None
        sources = None
        xx,yy,vx,vy = [],[],[],[]

        for i in range(5):
            imcopy = afwImage.ImageF(imorig, True)
            y = 25
            # Put all four sources in the parent (i==0), and one
            # source in each child (i=[1 to 4])
            if i in [0,1]:
                addPsf(imcopy, psf, 20, y, 1000)
            if i in [0,2]:
                addGaussian(imcopy, 40, y, 10, 3, 2e5)
            if i in [0,3]:
                addGaussian(imcopy, 75, y, 10, 3, 2e5)
            if i in [0,4]:
                addPsf(imcopy, psf, 95, y, 1000)

            # copy the pixels into the exposure object
            im <<= imcopy

            print 'Running detection...'
            if i == 0:
                detected = detect.makeSourceCatalog(table, exposure)
                sources = detected.sources
                print len(sources), 'sources'
                self.assertEqual(len(sources), 1)
            else:
                fpSets = detect.detectFootprints(exposure)
                print fpSets.numPos, 'detected'
                fpSets.positive.makeSources(sources)
                self.assertEqual(fpSets.numPos, 1)
                print len(sources), 'sources total'

            s = sources[-1]

            # Remember the centroid provided by the detection alg.
            #pk = s.getFootprint().getPeaks()[0]
            #print 'Initial peak pos', pk.getFx(), pk.getFy()
            #cx.append(pk.getFx())
            #cy.append(pk.getFy())

            print 'Running measurement...'
            measure.plotpat = 'single-%i.png' % i
            measure.run(exposure, sources[-1:])

            fp = s.getFootprint()
            if i == 0:
                # This is the real image
                fullim = imcopy
            else:
                print 'Creating heavy footprint...'
                heavy = afwDet.makeHeavyFootprint(fp, mi)
                s.setFootprint(heavy)

            #pk = s.getFootprint().getPeaks()[0]
            #print 'Final peak pos', pk.getFx(), pk.getFy()

            # Record the single-source measurements.
            xx.append(s.getX())
            yy.append(s.getY())
            vx.append(s.getIxx())
            vy.append(s.getIyy())

        parent = sources[0]
        kids = sources[1:]
        #print 'parent id', parent.getId()
        #print 'kid ids', [k.getId() for k in kids]

        parentid = parent.getId()
        pfp = parent.getFootprint()
        for s in kids:
            s.setParent(parentid)
            # Ensure that the parent footprint contains all the child footprints
            for span in s.getFootprint().getSpans():
                pfp.addSpan(span)
            #print 'kid:'
            #print '  ', s
            #print '  id', s.getId()
            #print '  parent', s.getParent()
            #print '  fp', s.getFootprint()
            #print '  hfp', afwDet.cast_HeavyFootprintF(s.getFootprint())
        pfp.normalize()
        parent.setFootprint(pfp)

        # Reset the centroids and blot out the measurements we're going to check...
        #for s,cxi,cyi in zip(sources, cx, cy):
            #print 's is', s
            #s.setX(cxi)
            #s.setY(cyi)
        #key = schema.find(measconf.slots.shape).key
        key = sources.getTable().getShapeKey()
        for s in sources:
            print 'shape', s.get(key)
            sh = s.get(key)
            sh.setIxx(np.nan)
            sh.setIyy(np.nan)
            sh.setIxy(np.nan)
            s.set(key, sh)

        for s in sources:
            print 'before final meas:', s.getIxx(), s.getIyy(), s.getIxy()
            
        im <<= fullim
        measconf.doRemoveOtherSources = True
        measure.plotpat = 'joint-%(sourcenum)i.png'
        measure.run(exposure, sources)

        fields = ['centroid.sdss', 'shape.sdss', ]
        keys = [schema.find(f).key for f in fields]

        flagkeys = [schema.find(x).key for x in [
            'flags.badcentroid',
            'shape.sdss.flags.maxiter', 'shape.sdss.flags.shift',
            'shape.sdss.flags.unweighted', 'shape.sdss.flags.unweightedbad']]
        flagabbr = ['C', 'M','S','U','B']

        xx2,yy2,vx2,vy2 = [],[],[],[]
        for s in sources:
            print 'Measured: id', s.getId()
            print '  centroid', s.getX(), s.getY()
            print '  shape(ixx,iyy)', s.getIxx(), s.getIyy()
            print '  flags',
            for ab,key in zip(flagabbr, flagkeys):
                if s.get(key):
                    print ab, ' '
            print

            xx2.append(s.getX())
            yy2.append(s.getY())
            vx2.append(s.getIxx())
            vy2.append(s.getIyy())

        print 'xx', xx, xx2
        print 'yy', yy, yy2
        print 'vx', vx, vx2
        print 'vy', vy, vy2


    ### FIXME
    def tst1(self):
        seed = 42
        rand = afwMath.Random(afwMath.Random.MT19937, seed)
        
        #for k in range(5):
        self.runone(1, rand)

    def getpsf(self):
        FWHM = 5
        ksize = 25
        psf = afwDet.createPsf("DoubleGaussian", ksize, ksize,
                               FWHM/(2*sqrt(2*log(2))), 1, 0.1)
        return psf

    def runone(self, kk, rand):
        psf = self.getpsf()
        
        im = afwImage.ImageF(120, 200)
        skystd = 100
        afwMath.randomGaussianImage(im, rand)
        im *= skystd

        x = 60
        y0 = 16
        ystep = 33
        for i in range(6):
            dx = [28,29,30, 35,36,37][i]
            y = y0 + i*ystep
            #                x   y  sx sy flux
            addGaussian(im, x, y, 10, 3, 2e5)
            #addGaussian(im, x, y, 10, 3, 1e5)
            addPsf(im, psf, x+dx, y, 1000)
            addPsf(im, psf, x-dx, y, 1000)

        #im.writeFits('im.fits')

        # This doesn't seem to be used...
        mi = afwImage.MaskedImageF(im)
        var = mi.getVariance()
        #print 'Variance:', var
        var = skystd**2
        exposure = afwImage.makeExposure(mi)
        exposure.setPsf(psf)

        detconf = measAlg.SourceDetectionConfig()
        detconf.reEstimateBackground = False

        measconf = measAlg.SourceMeasurementConfig()
        measconf.doApplyApCorr = False

        #newalgs = [ 'shape.hsm.ksb', 'shape.hsm.bj', 'shape.hsm.linear' ]
        #measconf.algorithms = list(measconf.algorithms.names) + newalgs
        #measconf.algorithms += newalgs
        #measconf.algorithms = list(measconf.algorithms.names)

        schema = afwTable.SourceTable.makeMinimalSchema()
        detect = measAlg.SourceDetectionTask(config=detconf, schema=schema)
        measure = measAlg.SourceMeasurementTask(config=measconf, schema=schema)

        print 'Running detection...'
        table = afwTable.SourceTable.make(schema)
        detected = detect.makeSourceCatalog(table, exposure)
        sources = detected.sources

        self.assertEqual(len(sources), 18)

        for jj in range(2):

            print 'Running measurement...'
            measure.run(exposure, sources)

            #fields = schema.getNames()
            #print 'Fields:', fields
            fields = ['centroid.sdss', 'shape.sdss',
                      #'shape.hsm.bj.moments',
                      #'shape.hsm.ksb.moments',
                      #'shape.hsm.linear.moments',
                      #'shape.sdss.flags.maxiter', 'shape.sdss.flags.shift',
                      #'shape.sdss.flags.unweighted', 'shape.sdss.flags.unweightedbad'
                      ]
            keys = [schema.find(f).key for f in fields]
            xx,yy,vx,vy = [],[],[],[]
            for source in sources:
                #print '  ', source
                #for f,k in zip(fields, keys):
                #    val = source.get(k)
                #    print '    ', f, val
                xx.append(source.getX())
                yy.append(source.getY())
                vx.append(source.getIxx())
                vy.append(source.getIyy())
        
            if plots:
                plotSources(im, sources, schema)
                plt.savefig('%i%s.png' % (kk, chr(ord('a')+jj)))

            #print 'vx:', vx
            # Sort, first vertically then horizontally
            iy = [int(round((y - y0) / float(ystep))) for y in yy]
            iy = np.array(iy)
            xx = np.array(xx)
            vx = np.array(vx)
            vy = np.array(vy)
            I = np.argsort(iy * 1000 + xx)
            vx = vx[I]
            vy = vy[I]
            #print 'vx', vx
            vx = vx[slice(1, 18, 3)]
            #print 'vx', vx

            bad = vx[:3]
            good = vx[3:]

            # When SdssShape fails, we get variance ~ 11

            I = np.flatnonzero(bad > 50.)
            self.assertTrue(len(I) > 0)
            bad = bad[I]
            I = np.flatnonzero(good > 50.)
            self.assertTrue(len(I) > 0)
            good = good[I]

            print 'bad:', bad
            print 'good:', good

            # bad: [ 209.78476672  192.35271583  176.76274525]
            # good: [  99.40557099  110.5701382 ]

            oklo,okhi = 80,120
            self.assertTrue(all((good > oklo) * (good < okhi)))
            if jj == 0:
                # Without "doRemoveOtherSources", we expect to find the variances
                # overestimated.
                self.assertTrue(all(bad > okhi))
            else:
                # With "doRemoveOtherSources", no problem!
                self.assertTrue(all((bad > oklo) * (bad < okhi)))

            measconf.doRemoveOtherSources = True




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


            
