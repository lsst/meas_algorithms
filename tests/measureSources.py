#!/usr/bin/env python
"""
Tests for measuring things

Run with:
   python measureSources.py
or
   python
   >>> import measureSources; measureSources.run()
"""
import sys
import itertools
import math
import unittest
import numpy
import lsst.utils.tests as tests
import lsst.pex.exceptions as pexExceptions
import lsst.pex.logging as pexLogging
import lsst.pex.policy as pexPolicy
import lsst.afw.detection as afwDetection
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.coord as afwCoord
import lsst.meas.algorithms as measAlg

try:
    type(verbose)
except NameError:
    verbose = 0
pexLogging.Trace_setVerbosity("afwDetection.Measure", verbose)

try:
    type(display)
except NameError:
    display = False

import lsst.afw.display.ds9 as ds9

FwhmPerSigma = 2*math.sqrt(2*math.log(2)) # FWHM for an N(0, 1) Gaussian

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class MeasureSourcesTestCase(unittest.TestCase):
    """A test case for Measure"""

    def setUp(self):
        pass
        
    def tearDown(self):
        pass

    def testNaiveMeasure(self):
        mi = afwImage.MaskedImageF(afwGeom.ExtentI(100, 200))
        mi.set(10)
        #
        # Create our measuring engine
        #
        exp = afwImage.makeExposure(mi)
        x0, y0 = 1234, 5678
        exp.setXY0(afwGeom.Point2I(x0, y0))
        
        control = measAlg.NaiveFluxControl()
        control.radius = 10.0
        schema = afwTable.SourceTable.makeMinimalSchema()
        mp = measAlg.MeasureSourcesBuilder().addAlgorithm(control).build(schema)
        table = afwTable.SourceTable.make(schema)
        source = table.makeRecord()
        mp.apply(source, exp, afwGeom.Point2D(30 + x0, 50 + y0))
        flux = 3170.0
        self.assertEqual(source.get(control.name), flux)

    def testApertureMeasure(self):
        mi = afwImage.MaskedImageF(afwGeom.ExtentI(100, 200))
        mi.set(10)
        #
        # Create our measuring engine
        #

        radii =  ( 1.0,   5.0,   10.0)  # radii to use
        fluxes = [50.0, 810.0, 3170.0]  # corresponding correct fluxes
        
        control = measAlg.ApertureFluxControl()
        control.radii = radii
        
        exp = afwImage.makeExposure(mi)
        x0, y0 = 1234, 5678
        exp.setXY0(afwGeom.Point2I(x0, y0))

        schema = afwTable.SourceTable.makeMinimalSchema()
        mp = measAlg.MeasureSourcesBuilder().addAlgorithm(control).build(schema)
        table = afwTable.SourceTable.make(schema)
        source = table.makeRecord()

        mp.apply(source, exp, afwGeom.Point2D(30 + x0, 50 + y0))
        measured = source[control.name]
        for i, f in enumerate(fluxes):
            self.assertEqual(f, measured[i])

    def testEllipticalGaussian(self):
        """Test measuring the properties of an elliptical Gaussian"""

        width, height = 200, 200
        xcen, ycen = 0.5*width, 0.5*height
        #
        # Make the object
        #
        gal = afwImage.ImageF(afwGeom.ExtentI(width, height))
        a, b, theta = float(10), float(5), 20
        flux = 1e4
        I0 = flux/(2*math.pi*a*b)

        c, s = math.cos(math.radians(theta)), math.sin(math.radians(theta))
        for y in range(height):
            for x in range(width):
                dx, dy = x - xcen, y - ycen
                u =  c*dx + s*dy
                v = -s*dx + c*dy
                val = I0*math.exp(-0.5*((u/a)**2 + (v/b)**2))
                if val < 0:
                    val = 0
                gal.set(x, y, val)

        objImg = afwImage.makeExposure(afwImage.makeMaskedImage(gal))
        objImg.getMaskedImage().getVariance().set(1.0)
        del gal
        objImg.setXY0(afwGeom.Point2I(1234, 5678))
        #
        # We need a PSF to be able to centroid well.  Cf. #2540
        #
        FWHM = 5
        ksize = 25                      # size of desired kernel
        objImg.setPsf(measAlg.DoubleGaussianPsf(ksize, ksize,
                                             FWHM/(2*math.sqrt(2*math.log(2))), 1, 0.1))
        

        if display:
            frame = 0
            ds9.mtv(objImg, frame=frame, title="Elliptical")

        self.assertAlmostEqual(1.0, afwMath.makeStatistics(objImg.getMaskedImage().getImage(),
                                                           afwMath.SUM).getValue()/flux)
        #
        # Test elliptical apertures
        #
        #
        msConfig = measAlg.SourceMeasurementConfig()
        msConfig.algorithms.names.add("flux.aperture.elliptical")
        radii = math.sqrt(a*b)*numpy.array([0.45, 1.0, 2.0, 3.0, 10.0,])

        msConfig.algorithms["flux.aperture.elliptical"].radii = radii
        schema = afwTable.SourceTable.makeMinimalSchema()
        ms = msConfig.makeMeasureSources(schema)
        
        table = afwTable.SourceTable.make(schema)
        msConfig.slots.setupTable(table)
        source = table.makeRecord()

        ss = afwDetection.FootprintSet(objImg.getMaskedImage(), afwDetection.Threshold(0.1))
        fp = ss.getFootprints()[0]
        source.setFootprint(fp)

        center =  fp.getPeaks()[0].getF()
        ms.apply(source, objImg, center)

        self.assertEqual(source.get("flux.aperture.elliptical.nProfile"), len(radii))

        r0 = 0.0
        if display:
            shape = source.getShape().clone()
            xy = afwGeom.ExtentD(source.getCentroid()) - afwGeom.ExtentD(objImg.getXY0())
            ds9.dot("x", xcen, ycen, ctype=ds9.RED)
            ds9.dot("+", *xy, frame=frame)
        with ds9.Buffering():
            for r, apFlux in zip(radii, source.get("flux.aperture.elliptical")):
                if display:                 # draw the inner and outer boundaries of the aperture
                    shape.scale(r/shape.getDeterminantRadius())
                    ds9.dot(shape, *xy, frame=frame)

                trueFlux = flux*(math.exp(-r0**2/(2*a*b)) - math.exp(-r**2/(2*a*b)))
                if verbose:
                    print "%5.2f %6.3f%%" % (r, 100*((trueFlux - apFlux)/flux))
                self.assertAlmostEqual(trueFlux/flux, apFlux/flux, 5)
                r0 = r
        #
        # Now measure some annuli "by hand" (we'll repeat this will EllipticalAperture algorithm soon)
        #

        for r1, r2 in [(0.0,    0.45*a),
                       (0.45*a, 1.0*a),
                       ( 1.0*a, 2.0*a),
                       ( 2.0*a, 3.0*a),
                       ( 3.0*a, 5.0*a),
                       ( 3.0*a, 10.0*a),
                       ]:
            control = measAlg.SincFluxControl()
            control.radius1 = r1
            control.radius2 = r2
            control.angle = math.radians(theta)
            control.ellipticity = 1 - b/a

            schema = afwTable.SourceTable.makeMinimalSchema()
            mp = measAlg.MeasureSourcesBuilder().addAlgorithm(control).build(schema)
            table = afwTable.SourceTable.make(schema)
            source = table.makeRecord()

            if display:                 # draw the inner and outer boundaries of the aperture
                Mxx = 1
                Myy = (b/a)**2

                mxx, mxy, myy = c**2*Mxx + s**2*Myy, c*s*(Mxx - Myy), s**2*Mxx + c**2*Myy
                for r in (r1, r2):
                    ds9.dot("@:%g,%g,%g" % (r**2*mxx, r**2*mxy, r**2*myy), xcen, ycen, frame=frame)

            mp.apply(source, objImg, center)

            self.assertAlmostEqual(math.exp(-0.5*(r1/a)**2) - math.exp(-0.5*(r2/a)**2),
                                   source.get(control.name)/flux, 5)

        control = measAlg.GaussianFluxControl()

        schema = afwTable.SourceTable.makeMinimalSchema()
        mp = measAlg.MeasureSourcesBuilder().addAlgorithm(control).build(schema)
        table = afwTable.SourceTable.make(schema)
        source = table.makeRecord()

        objImg.setPsf(None)             # no Psf
        mp.apply(source, objImg, center)
        # we haven't provided a PSF, so the built-in aperture correction won't work...but we'll get
        # a result anyway
        # Note that flags.psffactor==True sets flags=True IFF we attempt aperture corrections
        self.assertEqual(source.get(control.name + ".flags"), False)
        self.assertEqual(source.get(control.name + ".flags.psffactor"), True)
        gflux = source.get(control.name)
        err = gflux/flux - 1
        if abs(err) > 1.5e-5:
            self.assertEqual(gflux, flux, ("%g, %g: error is %g" % (gflux, flux, err)))

    def testPeakLikelihoodFlux(self):
        """Test measurement with PeakLikelihoodFlux
        """
        # make mp: a flux measurer
        measControl = measAlg.PeakLikelihoodFluxControl()
        schema = afwTable.SourceTable.makeMinimalSchema()
        mp = measAlg.MeasureSourcesBuilder().addAlgorithm(measControl).build(schema)
 
        # make and measure a series of exposures containing just one star, approximately centered
        bbox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.Extent2I(100, 101))
        kernelWidth = 35
        var = 100
        fwhm = 3.0
        sigma = fwhm/FwhmPerSigma
        convolutionControl = afwMath.ConvolutionControl()
        psf = measAlg.SingleGaussianPsf(kernelWidth, kernelWidth, sigma)
        psfKernel = psf.getLocalKernel()
        psfImage = psf.computeKernelImage()
        sumPsfSq = numpy.sum(psfImage.getArray()**2)
        psfSqArr = psfImage.getArray()**2
        for flux in (1000, 10000):
            ctrInd = afwGeom.Point2I(50, 51)
            ctrPos = afwGeom.Point2D(ctrInd)

            kernelBBox = psfImage.getBBox(afwImage.PARENT)
            kernelBBox.shift(afwGeom.Extent2I(ctrInd))

            # compute predicted flux error
            unshMImage = makeFakeImage(bbox, [ctrPos], [flux], fwhm, var)

            # filter image by PSF
            unshFiltMImage = afwImage.MaskedImageF(unshMImage.getBBox(afwImage.PARENT))
            afwMath.convolve(unshFiltMImage, unshMImage, psfKernel, convolutionControl)
            
            # compute predicted flux = value of image at peak / sum(PSF^2)
            # this is a sanity check of the algorithm, as much as anything
            predFlux = unshFiltMImage.getImage().get(ctrInd[0], ctrInd[1]) / sumPsfSq
            self.assertLess(abs(flux - predFlux), flux * 0.01)
            
            # compute predicted flux error based on filtered pixels
            # = sqrt(value of filtered variance at peak / sum(PSF^2)^2)
            predFluxErr = math.sqrt(unshFiltMImage.getVariance().get(ctrInd[0], ctrInd[1])) / sumPsfSq

            # compute predicted flux error based on unfiltered pixels
            # = sqrt(sum(unfiltered variance * PSF^2)) / sum(PSF^2)
            # and compare to that derived from filtered pixels;
            # again, this is a test of the algorithm
            varView = afwImage.ImageF(unshMImage.getVariance(), kernelBBox)
            varArr = varView.getArray()
            unfiltPredFluxErr = math.sqrt(numpy.sum(varArr*psfSqArr)) / sumPsfSq
            self.assertLess(abs(unfiltPredFluxErr - predFluxErr), predFluxErr * 0.01)
            
            for fracOffset in (afwGeom.Extent2D(0, 0), afwGeom.Extent2D(0.2, -0.3)):
                adjCenter = ctrPos + fracOffset
                if fracOffset == (0, 0):
                    maskedImage = unshMImage
                    filteredImage = unshFiltMImage
                else:
                    maskedImage = makeFakeImage(bbox, [adjCenter], [flux], fwhm, var)
                    # filter image by PSF
                    filteredImage = afwImage.MaskedImageF(maskedImage.getBBox(afwImage.PARENT))
                    afwMath.convolve(filteredImage, maskedImage, psfKernel, convolutionControl)

                exposure = afwImage.makeExposure(filteredImage)
                exposure.setPsf(psf)
                
                table = afwTable.SourceTable.make(schema)
                source = table.makeRecord()
                mp.apply(source, exposure, afwGeom.Point2D(*adjCenter))
                measFlux = source.get(measControl.name)
                measFluxErr = source.get(measControl.name + ".err")
                self.assertFalse(source.get(measControl.name + ".flags"))
                self.assertLess(abs(measFlux - flux), flux * 0.003)
                
                self.assertLess(abs(measFluxErr - predFluxErr), predFluxErr * 0.2)

                # try nearby points and verify that the flux is smaller;
                # this checks that the sub-pixel shift is performed in the correct direction
                for dx in (-0.2, 0, 0.2):
                    for dy in (-0.2, 0, 0.2):
                        if dx == dy == 0:
                            continue
                        offsetCtr = afwGeom.Point2D(adjCenter[0] + dx, adjCenter[1] + dy)
                        table = afwTable.SourceTable.make(schema)
                        source = table.makeRecord()
                        mp.apply(source, exposure, offsetCtr)
                        offsetFlux = source.get(measControl.name)
                        self.assertLess(offsetFlux, measFlux)
        
        # source so near edge of image that PSF does not overlap exposure should result in failure
        
        for edgePos in (
            (1, 50),
            (50, 1),
            (50, bbox.getHeight() - 1),
            (bbox.getWidth() - 1, 50),
        ):
            table = afwTable.SourceTable.make(schema)
            source = table.makeRecord()
            mp.apply(source, exposure, afwGeom.Point2D(*edgePos))
            self.assertTrue(source.get(measControl.name + ".flags"))
        
        # no PSF should result in failure: flags set
        noPsfExposure = afwImage.ExposureF(filteredImage)
        table = afwTable.SourceTable.make(schema)
        source = table.makeRecord()
        mp.apply(source, noPsfExposure, afwGeom.Point2D(*adjCenter))
        self.assertTrue(source.get(measControl.name + ".flags"))

    def testPixelFlags(self):
        width, height = 100, 100
        mi = afwImage.MaskedImageF(width, height)
        exp = afwImage.makeExposure(mi)
        mi.getImage().set(0)
        mask = mi.getMask()
        sat = mask.getPlaneBitMask('SAT')
        interp = mask.getPlaneBitMask('INTRP')
        edge = mask.getPlaneBitMask('EDGE')
        bad = mask.getPlaneBitMask('BAD')
        nodata = mask.getPlaneBitMask('NO_DATA')
        mask.set(0)
        mask.set(20, 20, sat)
        mask.set(60, 60, interp)
        mask.set(40, 20, bad)
        mask.set(20, 80, nodata)
        mask.Factory(mask, afwGeom.Box2I(afwGeom.Point2I(0,0), afwGeom.Extent2I(3, height))).set(edge)

        x0, y0 = 1234, 5678
        exp.setXY0(afwGeom.Point2I(x0, y0))

        control = measAlg.PixelFlagControl()
        schema = afwTable.SourceTable.makeMinimalSchema()
        mp = measAlg.MeasureSourcesBuilder().addAlgorithm(control).build(schema)
        table = afwTable.SourceTable.make(schema)

        allFlags = ["flags.pixel.edge",
                    "flags.pixel.bad",
                    "flags.pixel.saturated.center",
                    "flags.pixel.saturated.any",
                    "flags.pixel.interpolated.center",
                    "flags.pixel.interpolated.any",
                    ]
        for x, y, setFlags in [(1, 50, ["flags.pixel.edge"]),
                               (40, 20, ["flags.pixel.bad"]),
                               (20, 20, ["flags.pixel.saturated.center", "flags.pixel.saturated.any"]),
                               (20, 22, ["flags.pixel.saturated.any"]),
                               (60, 60, ["flags.pixel.interpolated.center", "flags.pixel.interpolated.any"]),
                               (60, 62, ["flags.pixel.interpolated.any"]),
                               (float("NAN"), 50, ["flags.pixel.edge"]),
                               (20, 80, ["flags.pixel.edge"]),
                               ]:
            source = table.makeRecord()
            foot = afwDetection.Footprint(afwGeom.Point2I(afwGeom.Point2D(x + x0, y + y0)), 5)
            source.setFootprint(foot)
            mp.apply(source, exp, afwGeom.Point2D(x + x0, y + y0))
            for flag in allFlags:
                value = source.get(flag)
                if flag in setFlags:
                    self.assertTrue(value, "Flag %s should be set for %f,%f" % (flag, x, y))
                else:
                    self.assertFalse(value, "Flag %s should not be set for %f,%f" % (flag, x, y))


class ForcedMeasureSourcesTestCase(unittest.TestCase):
    """A test case for making forced measurements"""

    def setUp(self):
        # We create an image that has a ramp (unique values for each pixel),
        # with a single high pixel that allows for centering
        self.width, self.height = 50, 50
        self.xcen, self.ycen = self.width//2, self.height//2
        self.image = afwImage.ImageF(afwGeom.ExtentI(self.width, self.height))
        for y in range(self.height):
            for x in range(self.width):
                self.image.set(x, y, self.width * y + x)
        self.image.set(self.xcen, self.ycen, 1234567.89)
        self.exp = afwImage.makeExposure(afwImage.makeMaskedImage(self.image))
        self.exp.getMaskedImage().getVariance().set(1.0)
        scale = 0.2 / 3600
        wcs = afwImage.makeWcs(afwCoord.Coord(0 * afwGeom.degrees, 0 * afwGeom.degrees),
                               afwGeom.Point2D(self.xcen, self.ycen), scale, 0, 0, scale)
        self.exp.setWcs(wcs)

        if display:
            frame = 1
            ds9.mtv(self.exp, frame=frame, title="Single pixel")

        # We will use a GaussianCentroid to tweak the center (it should not, for forced measurement)
        # and a NaiveFlux to measure the single pixel.  We'll start offset from the high pixel,
        # so that a forced measurement should yield a flux of zero, while a measurement that was allowed to
        # center should yield a flux of unity.
        # Note that previous versions used NaiveCentroid, which was so nonrobust that it failed to get
        # right answer when the input value was round-tripped through Wcs and modified by ~1E-8.
        gaussianCentroid = measAlg.GaussianCentroidControl()
        naiveFlux = measAlg.NaiveFluxControl()
        naiveFlux.radius = 0.5
        self.x, self.y = self.xcen - 1, self.ycen - 1

        self.foot = afwDetection.Footprint(afwGeom.Point2I(self.x, self.y), 2)
        peak = afwDetection.Peak(self.x, self.y)
        self.foot.getPeaks().push_back(peak)

        schema = afwTable.SourceTable.makeMinimalSchema()
        msb = measAlg.MeasureSourcesBuilder()
        msb.addAlgorithm(naiveFlux)
        msb.setCentroider(gaussianCentroid)
        self.measurer = msb.build(schema)
        self.table = afwTable.SourceTable.make(schema)
        self.table.defineCentroid("centroid.gaussian")

        schemaF = afwTable.SourceTable.makeMinimalSchema()
        msbF = measAlg.MeasureSourcesBuilder("", True)
        msbF.addAlgorithm(naiveFlux)
        msbF.setCentroider(gaussianCentroid)
        self.measurerF = msbF.build(schemaF)
        self.tableF = afwTable.SourceTable.make(schemaF)
        self.tableF.defineCentroid("centroid.gaussian")

    def tearDown(self):
        del self.image
        del self.exp
        del self.measurer
        del self.measurerF
        del self.table
        del self.tableF
        del self.foot

    def makeSource(self):
        return self.table.makeRecord()

    def makeSourceF(self):
        return self.tableF.makeRecord()

    def checkForced(self, source, forced):
        """Check whether the forced photometry was done with centering or not"""
        self.assertEqual(source.get("flux.naive"),
                         self.image.get(self.x, self.y) if forced else self.image.get(self.xcen, self.ycen))

    def testExplicit(self):
        sys.stderr.write("Explicit\n"); sys.stderr.flush()
        # Explicit center, without refinement
        source = self.makeSource()
        source.setFootprint(self.foot)
        self.measurer.apply(source, self.exp, afwGeom.Point2D(self.x, self.y), False)
        self.checkForced(source, True)

        # Explicit center, with refinement
        source = self.makeSource()
        source.setFootprint(self.foot)
        self.measurer.apply(source, self.exp, afwGeom.Point2D(self.x, self.y), True)
        self.checkForced(source, False)

    def testWithPeak(self):
        sys.stderr.write("WithPeak\n"); sys.stderr.flush()
        # Start with peak, don't refine
        source = self.makeSource()
        source.setFootprint(self.foot)
        self.measurer.applyWithPeak(source, self.exp, False)
        self.checkForced(source, True)

        # Normal use (single frame measurement): center up on peak
        source = self.makeSource()
        source.setFootprint(self.foot)
        self.measurer.applyWithPeak(source, self.exp, True)
        self.checkForced(source, False)

    def testWithPixel(self):
        sys.stderr.write("WithPixel\n"); sys.stderr.flush()
        # Center defined by previous centroid
        source = self.makeSource()
        source.set(self.table.getCentroidKey(), afwGeom.Point2D(self.x, self.y))
        self.measurer.applyWithPixel(source, self.exp, False)
        self.checkForced(source, True)

        # Start with previous centroid, refine
        source = self.makeSource()
        source.set(self.table.getCentroidKey(), afwGeom.Point2D(self.x, self.y))
        self.measurer.applyWithPixel(source, self.exp, True)
        self.checkForced(source, False)

    def testWithCoord(self):
        sys.stderr.write("WithCoord\n"); sys.stderr.flush()
        # Center defined by previous coordinates
        source = self.makeSource()
        source.setCoord(self.exp.getWcs().pixelToSky(afwGeom.Point2D(self.x, self.y)))
        self.measurer.applyWithCoord(source, self.exp, False)
        self.checkForced(source, True)

        # Center defined by previous coordinates, with refinement
        source = self.makeSource()
        source.setCoord(self.exp.getWcs().pixelToSky(afwGeom.Point2D(self.x, self.y)))
        self.measurer.applyWithCoord(source, self.exp, True)
        self.checkForced(source, False)

    def testWithReference(self):
        # Center defined by reference source
        wcs = self.exp.getWcs()
        wcs.flipImage(True, True, self.exp.getDimensions())
        
        ref = self.makeSourceF()
        ref.setFootprint(self.foot)
        ref.setCoord(self.exp.getWcs().pixelToSky(afwGeom.Point2D(self.x, self.y)))
        
        source = self.makeSourceF()
        self.measurerF.applyForced(source, self.exp, ref, wcs, False)
        self.checkForced(source, True)

        # Center defined by reference source
        wcs = self.exp.getWcs()
        wcs.flipImage(True, True, self.exp.getDimensions())
        
        ref = self.makeSourceF()
        ref.setFootprint(self.foot)
        ref.setCoord(self.exp.getWcs().pixelToSky(afwGeom.Point2D(self.x, self.y)))
        
        source = self.makeSourceF()
        self.measurerF.applyForced(source, self.exp, ref, wcs, True)
        self.checkForced(source, False)

def addStar(image, center, flux, fwhm):
    """Add a perfect single Gaussian star to an image
    
    @warning uses Python to iterate over all pixels (because there is no C++
    function that computes a Gaussian offset by a non-integral amount).
    
    @param[in,out] image: Image to which to add star
    @param[in] center: position of center of star on image (pair of float)
    @param[in] flux: flux of Gaussian star, in counts
    @param[in] fwhm: FWHM of Gaussian star, in pixels
    """
    sigma = fwhm/FwhmPerSigma
    func = afwMath.GaussianFunction2D(sigma, sigma, 0)
    starImage = afwImage.ImageF(image.getBBox(afwImage.PARENT))
    # The flux in the region of the image will not be exactly the desired flux because the Gaussian
    # does not extend to infinity, so keep track of the actual flux and correct for it
    actFlux = 0
    # No function exists that has a fractional x and y offset, so set the image the slow way
    for i in range(image.getWidth()):
        x = center[0] - i
        for j in range(image.getHeight()):
            y = center[1] - j
            pixVal = flux * func(x, y)
            actFlux += pixVal
            starImage[i, j] += pixVal
    starImage *= flux / actFlux
    
    image += starImage

def makeFakeImage(bbox, centerList, fluxList, fwhm, var):
    """Make a fake image containing a set of stars variance = image + var
    
    (It is trivial to add Poisson noise, which would be more accurate,
    but hard to make a unit test  that can reliably determine whether such an image passes a test)

    @param[in] bbox: bounding box for image
    @param[in] centerList: list of positions of center of star on image (pairs of float)
    @param[in] fluxList: flux of each star, in counts
    @param[in] fwhm: FWHM of Gaussian star, in pixels
    @param[in] var: value of variance plane (counts)
    """
    if len(centerList) != len(fluxList):
        raise RuntimeError("len(centerList) != len(fluxList)")
    maskedImage = afwImage.MaskedImageF(bbox)
    image = maskedImage.getImage()
    for center, flux in itertools.izip(centerList, fluxList):
        addStar(image, center=center, flux=flux, fwhm=fwhm)
    variance = maskedImage.getVariance()
    variance[:] = image
    variance += var
    return maskedImage


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    #suites += unittest.makeSuite(MeasureSourcesTestCase)
    suites += unittest.makeSuite(ForcedMeasureSourcesTestCase)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
