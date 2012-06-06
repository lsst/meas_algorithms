#!/usr/bin/env python
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

'''
Ticket 2019:

Before running the measurement algorithms, we should replace other
detection footprints with noise, then put footprints in one at a time
and measure.
'''

class Ticket2139TestCase(unittest.TestCase):
    def test1(self):
        task = measAlg.SourceMeasurementTask()


#plots = True
plots = False
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
        noiseim = imorig

        mi = afwImage.MaskedImageF(im)
        mi.getVariance().set(skystd**2)
        exposure = afwImage.makeExposure(mi)
        exposure.setPsf(psf)

        detconf = measAlg.SourceDetectionConfig()
        detconf.reEstimateBackground = False
        measconf = measAlg.SourceMeasurementConfig()
        measconf.doApplyApCorr = False
        measconf.doRemoveOtherSources = True
        measconf.noiseSeed = 42

        schema = afwTable.SourceTable.makeMinimalSchema()
        detect = measAlg.SourceDetectionTask(config=detconf, schema=schema)
        measure = MySourceMeasurementTask(config=measconf, schema=schema,
                                          doplot=plots)
        table = afwTable.SourceTable.make(schema)
        table.preallocate(10)

        # We're going to fake up a perfect deblend hierarchy here, by
        # creating individual images containing single sources and
        # measuring them, and then creating a deblend hierarchy where
        # the children have the correct HeavyFootprints.  We want to
        # find that the measurements on the deblend hierarchy and the
        # blended image are equal to the individual images.
        #
        # Note that in the normal setup we don't expect the
        # measurements to be *identical* because of the faint wings of
        # the objects; when measuring a deblended child, we pick up
        # the wings of the other objects.
        #
        # In order to get exactly equal measurements, we'll fake some
        # sources that have no wings -- we'll copy just the source
        # pixels within the footprint.  This means that all the
        # footprints are the same, and the pixels inside the footprint
        # are the same.
        
        fullim = None
        sources = None
        # "normal" measurements
        xx0,yy0,vx0,vy0 = [],[],[],[]
        # "no-wing" measurements
        xx1,yy1,vx1,vy1 = [],[],[],[]

        y = 25
        for i in range(5):
            # no-noise source image
            sim = afwImage.ImageF(imorig.getWidth(), imorig.getHeight())
            # Put all four sources in the parent (i==0), and one
            # source in each child (i=[1 to 4])
            if i in [0,1]:
                addPsf(sim, psf, 20, y, 1000)
            if i in [0,2]:
                addGaussian(sim, 40, y, 10, 3, 2e5)
            if i in [0,3]:
                addGaussian(sim, 75, y, 10, 3, 2e5)
            if i in [0,4]:
                addPsf(sim, psf, 95, y, 1000)
            imcopy = afwImage.ImageF(imorig, True)
            imcopy += sim
            # copy the pixels into the exposure object
            im <<= imcopy

            if i == 0:
                detected = detect.makeSourceCatalog(table, exposure)
                sources = detected.sources
                print 'detected', len(sources), 'sources'
                self.assertEqual(len(sources), 1)
            else:
                fpSets = detect.detectFootprints(exposure)
                print 'detected', fpSets.numPos, 'sources'
                fpSets.positive.makeSources(sources)
                self.assertEqual(fpSets.numPos, 1)
                print len(sources), 'sources total'

            measure.plotpat = 'single-%i.png' % i
            measure.run(exposure, sources[-1:])
            s = sources[-1]
            fp = s.getFootprint()
            if i == 0:
                # This is the blended image
                fullim = imcopy
            else:
                print 'Creating heavy footprint...'
                heavy = afwDet.makeHeavyFootprint(fp, mi)
                s.setFootprint(heavy)

            # Record the single-source measurements.
            xx0.append(s.getX())
            yy0.append(s.getY())
            vx0.append(s.getIxx())
            vy0.append(s.getIyy())

            # "no-wings": add just the source pixels within the footprint
            im <<= sim
            h = afwDet.makeHeavyFootprint(fp, mi)
            sim2 = afwImage.ImageF(imorig.getWidth(), imorig.getHeight())
            h.insert(sim2)
            imcopy = afwImage.ImageF(imorig, True)
            imcopy += sim2
            im <<= imcopy
            measure.plotpat = 'single2-%i.png' % i
            measure.run(exposure, sources[i:i+1], noiseImage=noiseim)
            s = sources[i]
            xx1.append(s.getX())
            yy1.append(s.getY())
            vx1.append(s.getIxx())
            vy1.append(s.getIyy())
            if i == 0:
                fullim2 = imcopy

        # Now we'll build the fake deblended hierarchy.
        parent = sources[0]
        kids = sources[1:]
        # Ensure that the parent footprint contains all the child footprints
        pfp = parent.getFootprint()
        for s in kids:
            for span in s.getFootprint().getSpans():
                pfp.addSpan(span)
        pfp.normalize()
        #parent.setFootprint(pfp)
        # The parent-child relationship is established through the IDs
        parentid = parent.getId()
        for s in kids:
            s.setParent(parentid)

        # Reset all the measurements
        shkey = sources.getTable().getShapeKey()
        ckey = sources.getTable().getCentroidKey()
        for s in sources:
            sh = s.get(shkey)
            sh.setIxx(np.nan)
            sh.setIyy(np.nan)
            sh.setIxy(np.nan)
            s.set(shkey, sh)
            c = s.get(ckey)
            c.setX(np.nan)
            c.setY(np.nan)
            s.set(ckey, c)

        # Measure the "deblended" normal sources
        im <<= fullim
        measure.plotpat = 'joint-%(sourcenum)i.png'
        measure.run(exposure, sources)
        xx2,yy2,vx2,vy2 = [],[],[],[]
        for s in sources:
            xx2.append(s.getX())
            yy2.append(s.getY())
            vx2.append(s.getIxx())
            vy2.append(s.getIyy())

        # Measure the "deblended" no-wings sources
        im <<= fullim2
        measure.plotpat = 'joint2-%(sourcenum)i.png'
        measure.run(exposure, sources, noiseImage=noiseim)
        xx3,yy3,vx3,vy3 = [],[],[],[]
        for s in sources:
            xx3.append(s.getX())
            yy3.append(s.getY())
            vx3.append(s.getIxx())
            vy3.append(s.getIyy())

        print 'Normal:'
        print 'xx  ', xx0
        print '  vs', xx2
        print 'yy  ', yy0
        print '  vs', yy2
        print 'vx  ', vx0
        print '  vs', vx2
        print 'vy  ', vy0
        print '  vs', vy2

        print 'No wings:'
        print 'xx  ', xx1
        print '  vs', xx3
        print 'yy  ', yy1
        print '  vs', yy3
        print 'vx  ', vx1
        print '  vs', vx3
        print 'vy  ', vy1
        print '  vs', vy3

        # These "normal" tests are not very stringent.
        # 0.1-pixel centroids
        self.assertTrue(all([abs(v1-v2) < 0.1 for v1,v2 in zip(xx0,xx2)]))
        self.assertTrue(all([abs(v1-v2) < 0.1 for v1,v2 in zip(yy0,yy2)]))
        # 10% variances
        self.assertTrue(all([abs(v1-v2)/((v1+v2)/2.) < 0.1 for v1,v2 in zip(vx0,vx2)]))
        self.assertTrue(all([abs(v1-v2)/((v1+v2)/2.) < 0.1 for v1,v2 in zip(vy0,vy2)]))

        # The "no-wings" tests should be exact.
        self.assertTrue(xx1 == xx3)
        self.assertTrue(yy1 == yy3)
        self.assertTrue(vx1 == vx3)
        self.assertTrue(vy1 == vy3)

        # Reset sources
        for s in sources:
            sh = s.get(shkey)
            sh.setIxx(np.nan)
            sh.setIyy(np.nan)
            sh.setIxy(np.nan)
            s.set(shkey, sh)
            c = s.get(ckey)
            c.setX(np.nan)
            c.setY(np.nan)
            s.set(ckey, c)

        # Test that the parent/child order is unimportant.
        im <<= fullim2
        measure.doplot = False
        sources2 = sources.copy()
        perm = [2,1,0,3,4]
        for i,j in enumerate(perm):
            sources2[i] = sources[j]
            # I'm not convinced that HeavyFootprints get copied correctly...
            sources2[i].setFootprint(sources[j].getFootprint())
        measure.run(exposure, sources2, noiseImage=noiseim)
        # "measure.run" reorders the sources!
        xx3,yy3,vx3,vy3 = [],[],[],[]
        for s in sources:
            xx3.append(s.getX())
            yy3.append(s.getY())
            vx3.append(s.getIxx())
            vy3.append(s.getIyy())
        self.assertTrue(xx1 == xx3)
        self.assertTrue(yy1 == yy3)
        self.assertTrue(vx1 == vx3)
        self.assertTrue(vy1 == vy3)

        # Reset sources
        for s in sources:
            sh = s.get(shkey)
            sh.setIxx(np.nan)
            sh.setIyy(np.nan)
            sh.setIxy(np.nan)
            s.set(shkey, sh)
            c = s.get(ckey)
            c.setX(np.nan)
            c.setY(np.nan)
            s.set(ckey, c)

        # Test that it still works when the parent ID falls in the middle of
        # the child IDs.
        im <<= fullim2
        measure.doplot = False
        sources2 = sources.copy()
        parentid = 3
        ids = [parentid, 1,2,4,5]
        for i,s in enumerate(sources2):
            s.setId(ids[i])
            if i != 0:
                s.setParent(parentid)
            s.setFootprint(sources[i].getFootprint())
            
        measure.run(exposure, sources2, noiseImage=noiseim)
        # The sources get reordered!
        xx3,yy3,vx3,vy3 = [],[],[],[]
        xx3,yy3,vx3,vy3 = [0]*5,[0]*5,[0]*5,[0]*5
        for i,j in enumerate(ids):
            xx3[i] = sources2[j-1].getX()
            yy3[i] = sources2[j-1].getY()
            vx3[i] = sources2[j-1].getIxx()
            vy3[i] = sources2[j-1].getIyy()
        self.assertTrue(xx1 == xx3)
        self.assertTrue(yy1 == yy3)
        self.assertTrue(vx1 == vx3)
        self.assertTrue(vy1 == vy3)




    def test1(self):
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

        # The SDSS adaptive moments code seems sometimes to latch onto
        # an incorrect answer (maybe from a noise spike or something).
        # None of the flags seem to be set.  The result are variance
        # measurements a bit bigger than the PSF.  With different
        # noise draws the source values here will show this effect
        # (hence the loop in "test1" to try "runone" will different
        # noise draws).

        # The real point of this test case, though, is to show that
        # replacing other detections by noise results in better
        # measurements.  We do this by constructing a fake image
        # containing six rows.  In the top three rows, we have a
        # galaxy flanked by two stars that are far enough away that
        # they don't confuse the SDSS adaptive moments code.  In the
        # bottom three rows, they're close enough that the detections
        # don't merge, but the stars cause the variance of the galaxy
        # to be mis-estimated.  We want to show that with the
        # "doRemoveOtherSources" option, the measurements on the
        # bottom three improve.

        # If you love ASCII art (and who doesn't, really), the
        # synthetic image is going to look like this:
        #
        #    *     GGG     *
        #    *     GGG     *
        #    *     GGG     *
        #       *  GGG  *
        #       *  GGG  *
        #       *  GGG  *

        # We have three of each to work around the instability
        # mentioned above.

        x = 60
        y0 = 16
        ystep = 33
        for i in range(6):
            dx = [28,29,30, 35,36,37][i]
            y = y0 + i*ystep
            #                x y sx sy flux
            addGaussian(im, x, y, 10, 3, 2e5)
            addPsf(im, psf, x+dx, y, 1000)
            addPsf(im, psf, x-dx, y, 1000)

        #im.writeFits('im.fits')

        mi = afwImage.MaskedImageF(im)
        var = mi.getVariance()
        var.set(skystd**2)
        exposure = afwImage.makeExposure(mi)
        exposure.setPsf(psf)

        detconf = measAlg.SourceDetectionConfig()
        detconf.reEstimateBackground = False

        measconf = measAlg.SourceMeasurementConfig()
        measconf.doApplyApCorr = False

        #newalgs = [ 'shape.hsm.ksb', 'shape.hsm.bj', 'shape.hsm.linear' ]
        #measconf.algorithms = list(measconf.algorithms.names) + newalgs

        schema = afwTable.SourceTable.makeMinimalSchema()
        detect = measAlg.SourceDetectionTask(config=detconf, schema=schema)
        measure = measAlg.SourceMeasurementTask(config=measconf, schema=schema)

        print 'Running detection...'
        table = afwTable.SourceTable.make(schema)
        detected = detect.makeSourceCatalog(table, exposure)
        sources = detected.sources

        # We don't want the sources to be close enough that their
        # detection masks touch.
        self.assertEqual(len(sources), 18)

        # Run measurement with and without "doRemoveOtherSources"...
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

            # Now we want to find the galaxy variance measurements...
            # Sort, first vertically then horizontally
            # iy ~ row number
            iy = [int(round((y - y0) / float(ystep))) for y in yy]
            iy = np.array(iy)
            xx = np.array(xx)
            vx = np.array(vx)
            vy = np.array(vy)
            I = np.argsort(iy * 1000 + xx)
            vx = vx[I]
            vy = vy[I]
            # The "left" stars will be indices 0, 3, 6, ...
            # The galaxies will be 1, 4, 7, ...
            vx = vx[slice(1, 18, 3)]

            # Bottom three galaxies may be contaminated by the stars
            bad = vx[:3]
            # Top three should be clean
            good = vx[3:]

            # When SdssShape fails, we get variance ~ 11

            I = np.flatnonzero(bad > 50.)
            # Hope that we got at least one valid measurement
            self.assertTrue(len(I) > 0)
            bad = bad[I]
            I = np.flatnonzero(good > 50.)
            self.assertTrue(len(I) > 0)
            good = good[I]

            print 'bad:', bad
            print 'good:', good

            # Typical:
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

            # Set "doRemoveOtherSources" for the second time through the loop...
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


            
