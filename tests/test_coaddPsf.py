#
# LSST Data Management System
#
# Copyright 2008-2017  AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#
import unittest

import lsst.geom
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.meas.algorithms as measAlg
import lsst.pex.exceptions as pexExceptions
import lsst.utils.tests


def getPsfMoments(psf, point):
    #    import os, pdb; print "PID =", os.getpid(); pdb.set_trace()
    image = psf.computeImage(point)
    array = image.getArray()
    sumx2 = 0.0
    sumy2 = 0.0
    sumy = 0.0
    sumx = 0.0
    sum = 0.0
    for x in range(image.getWidth()):
        for y in range(image.getHeight()):
            f = array[y][x]
            sumx2 += x*x*f
            sumy2 += y*y*f
            sumx += x*f
            sumy += y*f
            sum += f
    xbar = sumx/sum
    ybar = sumy/sum
    mxx = sumx2 - 2*xbar*sumx + xbar*xbar*sum
    myy = sumy2 - 2*ybar*sumy + ybar*ybar*sum
    return sum, xbar, ybar, mxx, myy, image.getX0(), image.getY0()


def getPsfSecondMoments(psf, point):
    sum, xbar, ybar, mxx, myy, x0, y0 = getPsfMoments(psf, point)
    return mxx, myy


def makeBiaxialGaussianPsf(sizex, sizey, sigma1, sigma2, theta):
    kernel = afwMath.AnalyticKernel(sizex, sizey, afwMath.GaussianFunction2D(sigma1, sigma2, theta))
    return measAlg.KernelPsf(kernel)

# This is a mock method for coadding the moments of the component Psfs at a point
# Check that the coaddpsf passed in is really using the correct components and weighting them properly
# The components in this case are all single gaussians, and we will just add the moments
# If useValidPolygon = True then the exposures are expected to have validPolygons defined, otherwise
# it will set the whoe region as valid


def getCoaddSecondMoments(coaddpsf, point, useValidPolygon=False):
    count = coaddpsf.getComponentCount()
    coaddWcs = coaddpsf.getCoaddWcs()
    weight_sum = 0.0
    m1_sum = 0.0
    m2_sum = 0.0
    for i in range(count):
        wcs = coaddpsf.getWcs(i)
        psf = coaddpsf.getPsf(i)
        bbox = lsst.geom.Box2D(coaddpsf.getBBox(i))
        if useValidPolygon:
            validPolygon = coaddpsf.getValidPolygon(i)
        else:
            validPolygon = afwGeom.Polygon(bbox)

        point_rel = wcs.skyToPixel(coaddWcs.pixelToSky(lsst.geom.Point2D(point)))
        if bbox.contains(point_rel) and validPolygon.contains(point_rel):
            weight = coaddpsf.getWeight(i)
            m0, xbar, ybar, mxx, myy, x0, y0 = getPsfMoments(psf, point)  # , extent)
            m1_sum += mxx*weight
            m2_sum += myy*weight
            weight_sum += weight
    if weight_sum == 0.0:
        return 0, 0
    else:
        return m1_sum/weight_sum, m2_sum/weight_sum


class CoaddPsfTest(lsst.utils.tests.TestCase):

    def setUp(self):
        scale = 5.55555555e-05*lsst.geom.degrees
        self.cdMatrix = afwGeom.makeCdMatrix(scale=scale, flipX=True)
        self.crpix = lsst.geom.PointD(1000, 1000)
        self.crval = lsst.geom.SpherePoint(0.0, 0.0, lsst.geom.degrees)
        self.wcsref = afwGeom.makeSkyWcs(crpix=self.crpix, crval=self.crval, cdMatrix=self.cdMatrix)

        schema = afwTable.ExposureTable.makeMinimalSchema()
        self.weightKey = schema.addField("weight", type="D", doc="Coadd weight")
        self.mycatalog = afwTable.ExposureCatalog(schema)

    def tearDown(self):
        del self.crpix
        del self.crval
        del self.wcsref
        del self.weightKey
        del self.mycatalog

    #   This is a test which checks to see that all of the ExposureCatalog rows are correctly
    #   ingested by the CoaddPsf constructor, and that they can be read back in the right order
    #   and with the right values
    #   The weightname mechanism is also tested.  Whatever input column name is used should be
    #   mapped to "weight"

    def testCreate(self):
        """Check that we can create a CoaddPsf with 9 elements."""
        print("CreatePsfTest")

        # also test that the weight field name is correctly observed
        schema = afwTable.ExposureTable.makeMinimalSchema()
        schema.addField("customweightname", type="D", doc="Coadd weight")
        mycatalog = afwTable.ExposureCatalog(schema)

        # Each of the 9 has its peculiar Psf, Wcs, weight, and bounding box.
        for i in range(1, 10, 1):
            record = mycatalog.getTable().makeRecord()
            psf = measAlg.DoubleGaussianPsf(100, 100, i, 1.00, 0.0)
            record.setPsf(psf)
            crpix = lsst.geom.PointD(i*1000.0, i*1000.0)
            wcs = afwGeom.makeSkyWcs(crpix=crpix, crval=self.crval, cdMatrix=self.cdMatrix)

            record.setWcs(wcs)
            record['customweightname'] = 1.0 * (i+1)
            record['id'] = i
            bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(i*1000, i*1000))
            record.setBBox(bbox)
            mycatalog.append(record)

        # create the coaddpsf
        mypsf = measAlg.CoaddPsf(mycatalog, self.wcsref, 'customweightname')

        # check to be sure that we got the right number of components, in the right order
        self.assertEqual(mypsf.getComponentCount(), 9)
        for i in range(1, 10, 1):
            wcs = mypsf.getWcs(i-1)
            psf = mypsf.getPsf(i-1)
            bbox = mypsf.getBBox(i-1)
            weight = mypsf.getWeight(i-1)
            id = mypsf.getId(i-1)
            self.assertEqual(i, id)
            self.assertEqual(weight, 1.0*(i+1))
            self.assertEqual(bbox.getBeginX(), 0)
            self.assertEqual(bbox.getBeginY(), 0)
            self.assertEqual(bbox.getEndX(), 1000 * i)
            self.assertEqual(bbox.getEndY(), 1000 * i)
            self.assertAlmostEqual(wcs.getPixelOrigin().getX(), (1000.0 * i))
            self.assertAlmostEqual(wcs.getPixelOrigin().getY(), (1000.0 * i))
            m0, xbar, ybar, mxx, myy, x0, y0 = getPsfMoments(psf, lsst.geom.Point2D(0, 0))
            self.assertAlmostEqual(i*i, mxx, delta=0.01)
            self.assertAlmostEqual(i*i, myy, delta=0.01)

    def testFractionalPixel(self):
        """Check that we can create a CoaddPsf with 10 elements."""
        print("FractionalPixelTest")
        cdMatrix = afwGeom.makeCdMatrix(
            scale=5.55555555e-05*lsst.geom.degrees,
            orientation=90*lsst.geom.degrees,
        )
        wcs = afwGeom.makeSkyWcs(crpix=self.crpix, crval=self.crval, cdMatrix=cdMatrix)

        # make a single record with an oblong Psf
        record = self.mycatalog.getTable().makeRecord()
        psf = makeBiaxialGaussianPsf(100, 100, 6.0, 6.0, 0.0)
        record.setPsf(psf)
        record.setWcs(wcs)
        record['weight'] = 1.0
        record['id'] = 1
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(2000, 2000))
        record.setBBox(bbox)
        self.mycatalog.append(record)
        mypsf = measAlg.CoaddPsf(self.mycatalog, self.wcsref)
        psf.computeImage(lsst.geom.PointD(0.25, 0.75))
        psf.computeImage(lsst.geom.PointD(0.25, 0.75))
        psf.computeImage(lsst.geom.PointD(1000, 1000))
        m0, xbar, ybar, mxx, myy, x0, y0 = getPsfMoments(psf, lsst.geom.Point2D(0.25, 0.75))
        cm0, cxbar, cybar, cmxx, cmyy, cx0, cy0 = getPsfMoments(mypsf, lsst.geom.Point2D(0.25, 0.75))
        self.assertAlmostEqual(x0+xbar, cx0+cxbar, delta=0.01)
        self.assertAlmostEqual(y0+ybar, cy0+cybar, delta=0.01)

    def testRotatePsf(self):
        """Check that we can create a CoaddPsf with 10 elements."""
        print("RotatePsfTest")
        cdMatrix = afwGeom.makeCdMatrix(
            scale=5.55555555e-05*lsst.geom.degrees,
            orientation=90*lsst.geom.degrees,
            flipX=True,
        )
        wcs = afwGeom.makeSkyWcs(crpix=self.crpix, crval=self.crval, cdMatrix=cdMatrix)

        # make a single record with an oblong Psf
        record = self.mycatalog.getTable().makeRecord()
        psf = makeBiaxialGaussianPsf(100, 100, 1.0, 6.0, 0.0)
        record.setPsf(psf)
        record.setWcs(wcs)
        record['weight'] = 1.0
        record['id'] = 1
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(2000, 2000))
        record.setBBox(bbox)
        self.mycatalog.append(record)
        mypsf = measAlg.CoaddPsf(self.mycatalog, self.wcsref)
        m0, xbar, ybar, mxx, myy, x0, y0 = getPsfMoments(psf, lsst.geom.Point2D(0.25, 0.75))
        cm0, cxbar, cybar, cmxx, cmyy, cx0, cy0 = getPsfMoments(mypsf, lsst.geom.Point2D(0.25, 0.75))
        self.assertAlmostEqual(mxx, cmyy, delta=0.01)
        self.assertAlmostEqual(myy, cmxx, delta=0.01)

    def testDefaultSize(self):
        """Test of both default size and specified size."""
        print("DefaultSizeTest")
        sigma0 = 5
        # set the peak of the outer guassian to 0 so this is really a single gaussian.

        psf = measAlg.DoubleGaussianPsf(60, 60, 1.5*sigma0, 1, 0.0)

        # Now make the catalog
        record = self.mycatalog.getTable().makeRecord()
        psf = measAlg.DoubleGaussianPsf(100, 100, 10.0, 1.00, 1.0)
        record.setPsf(psf)
        wcs = afwGeom.makeSkyWcs(crpix=self.crpix, crval=self.crval, cdMatrix=self.cdMatrix)
        record.setWcs(wcs)
        record['weight'] = 1.0
        record['id'] = 1
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(2000, 2000))
        record.setBBox(bbox)
        self.mycatalog.append(record)

        mypsf = measAlg.CoaddPsf(self.mycatalog, self.wcsref)  # , 'weight')

        m1coadd, m2coadd = getCoaddSecondMoments(mypsf, lsst.geom.Point2D(0, 0))
        m1, m2 = getPsfSecondMoments(mypsf, lsst.geom.Point2D(1000, 1000))
        self.assertAlmostEqual(m1, m1coadd, delta=.01)
        self.assertAlmostEqual(m2, m2coadd, delta=.01)

    def testSimpleGaussian(self):
        """Check that we can measure a single Gaussian's attributes."""
        print("SimpleGaussianTest")
        sigma0 = 5
        # set the peak of the outer guassian to 0 so this is really a single gaussian.

        psf = measAlg.DoubleGaussianPsf(60, 60, 1.5*sigma0, 1, 0.0)

        sigma = [5, 6, 7, 8]  # 5 pixels is the same as a sigma of 1 arcsec.

        # lay down a simple pattern of four ccds, set in a pattern of 1000 pixels around the center
        offsets = [(1999, 1999), (1999, 0), (0, 0), (0, 1999)]

#       Imagine a ccd in each of positions +-1000 pixels from the center
        for i in range(4):
            record = self.mycatalog.getTable().makeRecord()
            psf = measAlg.DoubleGaussianPsf(100, 100, sigma[i], 1.00, 1.0)
            record.setPsf(psf)
            crpix = lsst.geom.PointD(offsets[i][0], offsets[i][1])
            wcs = afwGeom.makeSkyWcs(crpix=crpix, crval=self.crval, cdMatrix=self.cdMatrix)
            record.setWcs(wcs)
            record['weight'] = 1.0
            record['id'] = i
            bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(2000, 2000))
            record.setBBox(bbox)
            self.mycatalog.append(record)

        mypsf = measAlg.CoaddPsf(self.mycatalog, self.wcsref)  # , 'weight')
        m1coadd, m2coadd = getCoaddSecondMoments(mypsf, lsst.geom.Point2D(1000, 1000))

        m1, m2 = getPsfSecondMoments(mypsf, lsst.geom.Point2D(1000, 1000))
        self.assertAlmostEqual(m1, m1coadd, delta=.01)

        m1, m2 = getPsfSecondMoments(mypsf, lsst.geom.Point2D(1000, 1001))
        m1coadd, m2coadd = getCoaddSecondMoments(mypsf, lsst.geom.Point2D(1000, 1001))
        self.assertAlmostEqual(m1, m1coadd, delta=0.01)

#   This test checks to be sure that the weights are being applied correctly in doComputeImage
#   Since the 2nd moments are linear in the function value, we can simply weight the moments
#   and be sure that the resulting moments are correct

    def testWeight(self):
        """Check that we can measure a single Gaussian's attributes."""
        print("WeightTest")
        sigma0 = 5
        # set the peak of the outer guassian to 0 so this is really a single gaussian.

        psf = measAlg.DoubleGaussianPsf(60, 60, 1.5*sigma0, 1, 0.0)

        sigma = [5, 6, 7, 8]  # 5 pixels is the same as a sigma of 1 arcsec.

        # lay down a simple pattern of four ccds, set in a pattern of 1000 pixels around the center
        offsets = [(1999, 1999), (1999, 0), (0, 0), (0, 1999)]

#       Imagine a ccd in each of positions +-1000 pixels from the center
        for i in range(4):
            record = self.mycatalog.getTable().makeRecord()
            psf = measAlg.DoubleGaussianPsf(100, 100, sigma[i], 1.00, 0.0)
            record.setPsf(psf)
            crpix = lsst.geom.PointD(offsets[i][0], offsets[i][1])
            wcs = afwGeom.makeSkyWcs(crpix=crpix, crval=self.crval, cdMatrix=self.cdMatrix)

            # print out the coorinates of this supposed 2000x2000 ccd in wcsref coordinates
            record.setWcs(wcs)
            record['weight'] = 1.0 * (i+1)
            record['id'] = i
            bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(2000, 2000))
            record.setBBox(bbox)
            self.mycatalog.append(record)

        mypsf = measAlg.CoaddPsf(self.mycatalog, self.wcsref)  # , 'weight')

        m1, m2 = getPsfSecondMoments(mypsf, lsst.geom.Point2D(1000, 1000))
        m1coadd, m2coadd = getCoaddSecondMoments(mypsf, lsst.geom.Point2D(1000, 1000))
        self.assertAlmostEqual(m1, m1coadd, delta=0.01)

        m1, m2 = getPsfSecondMoments(mypsf, lsst.geom.Point2D(1000, 1001))
        m1coadd, m2coadd = getCoaddSecondMoments(mypsf, lsst.geom.Point2D(1000, 1001))
        self.assertAlmostEqual(m1, m1coadd, delta=0.01)

        m1, m2 = getCoaddSecondMoments(mypsf, lsst.geom.Point2D(1001, 1000))
        m1coadd, m2coadd = getCoaddSecondMoments(mypsf, lsst.geom.Point2D(1001, 1000))
        self.assertAlmostEqual(m1, m1coadd, delta=0.01)

    def testTicket2872(self):
        """Test that CoaddPsf.getAveragePosition() is always a position at which
        we can call computeImage().
        """
        scale = 0.2*lsst.geom.arcseconds
        cdMatrix = afwGeom.makeCdMatrix(scale=scale)
        wcs = afwGeom.makeSkyWcs(
            crpix=lsst.geom.Point2D(50, 50),
            crval=lsst.geom.SpherePoint(45.0, 45.0, lsst.geom.degrees),
            cdMatrix=cdMatrix,
        )
        kernel = measAlg.DoubleGaussianPsf(7, 7, 2.0).getKernel()
        psf1 = measAlg.KernelPsf(kernel, lsst.geom.Point2D(0, 50))
        psf2 = measAlg.KernelPsf(kernel, lsst.geom.Point2D(100, 50))
        record1 = self.mycatalog.addNew()
        record1.setPsf(psf1)
        record1.setWcs(wcs)
        record1.setD(self.weightKey, 1.0)
        record1.setBBox(lsst.geom.Box2I(lsst.geom.Point2I(-40, 0), lsst.geom.Point2I(40, 100)))
        record2 = self.mycatalog.addNew()
        record2.setPsf(psf2)
        record2.setWcs(wcs)
        record2.setD(self.weightKey, 1.0)
        record2.setBBox(lsst.geom.Box2I(lsst.geom.Point2I(60, 0), lsst.geom.Point2I(140, 100)))
        coaddPsf = measAlg.CoaddPsf(self.mycatalog, wcs)
        naiveAvgPos = lsst.geom.Point2D(50, 50)
        with self.assertRaises(pexExceptions.InvalidParameterError):
            coaddPsf.computeKernelImage(naiveAvgPos)
        # important test is that this doesn't throw:
        coaddPsf.computeKernelImage()

    def testValidPolygonPsf(self):
        """Demonstrate that we can use the validPolygon on Exposures in the CoaddPsf."""
        # Create 9 separate records, each with its own peculiar Psf, Wcs,
        # weight, bounding box, and valid region.
        for i in range(1, 10):
            record = self.mycatalog.getTable().makeRecord()
            record.setPsf(measAlg.DoubleGaussianPsf(100, 100, i, 1.00, 0.0))
            crpix = lsst.geom.PointD(1000-10.0*i, 1000.0-10.0*i)
            wcs = afwGeom.makeSkyWcs(crpix=crpix, crval=self.crval, cdMatrix=self.cdMatrix)
            record.setWcs(wcs)
            record['weight'] = 1.0*(i + 1)
            record['id'] = i
            record.setBBox(lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(1000, 1000)))
            validPolygon = afwGeom.Polygon(lsst.geom.Box2D(lsst.geom.Point2D(0, 0),
                                                           lsst.geom.Extent2D(i*100, i*100)))
            record.setValidPolygon(validPolygon)
            self.mycatalog.append(record)

        # Create the CoaddPsf and check at three different points to ensure that the validPolygon is working
        mypsf = measAlg.CoaddPsf(self.mycatalog, self.wcsref, 'weight')

        for position in [lsst.geom.Point2D(50, 50), lsst.geom.Point2D(500, 500), lsst.geom.Point2D(850, 850)]:
            m1coadd, m2coadd = getCoaddSecondMoments(mypsf, position, True)
            m1, m2 = getPsfSecondMoments(mypsf, position)
            self.assertAlmostEqual(m1, m1coadd, delta=0.01)
            self.assertAlmostEqual(m2, m2coadd, delta=0.01)

    def testGoodPix(self):
        """Demonstrate that we can goodPix information in the CoaddPsf."""
        bboxSize = lsst.geom.Extent2I(2000, 2000)
        schema = afwTable.ExposureTable.makeMinimalSchema()
        schema.addField("weight", type="D", doc="Coadd weight")
        schema.addField("goodpix", type="I", doc="Number of good pixels")
        mycatalog = afwTable.ExposureCatalog(schema)

        # Create several records, each with its own peculiar center and numGoodPixels.
        # Each PSF has the same shape and size, and the position offsets are small
        # relative to the FWHM, in order to make it easy to predict the resulting
        # weighted mean position.
        xwsum = 0
        ywsum = 0
        wsum = 0
        for i, (xOff, yOff, numGoodPix) in enumerate((
            (30.0, -20.0, 25),
            (32.0, -21.0, 10),
            (28.0, -19.0, 30),
        )):
            xwsum -= xOff * numGoodPix
            ywsum -= yOff * numGoodPix
            wsum += numGoodPix
            record = mycatalog.getTable().makeRecord()
            record.setPsf(measAlg.DoubleGaussianPsf(25, 25, 10, 1.00, 0.0))
            offPix = self.crpix + lsst.geom.Extent2D(xOff, yOff)
            wcs = afwGeom.makeSkyWcs(crpix=offPix, crval=self.crval, cdMatrix=self.cdMatrix)
            record.setWcs(wcs)
            record['weight'] = 1.0
            record['id'] = i
            record['goodpix'] = numGoodPix
            record.setBBox(lsst.geom.Box2I(lsst.geom.Point2I(0, 0), bboxSize))
            mycatalog.append(record)

        mypsf = measAlg.CoaddPsf(mycatalog, self.wcsref, 'weight')
        predPos = lsst.geom.Point2D(xwsum/wsum, ywsum/wsum)
        self.assertPairsAlmostEqual(predPos, mypsf.getAveragePosition())

    def testBBox(self):
        """Check that computeBBox returns same BBox as realized Kernel Image

        and resized raises a Not Implemented Error"""
        sigma0 = 5
        size = [50, 60, 70, 80]

        for i in range(4):
            record = self.mycatalog.getTable().makeRecord()
            psf = measAlg.DoubleGaussianPsf(size[i], size[i], sigma0, 1.00, 0.0)
            record.setPsf(psf)
            wcs = afwGeom.makeSkyWcs(crpix=self.crpix, crval=self.crval, cdMatrix=self.cdMatrix)
            record.setWcs(wcs)
            record['weight'] = 1.0 * (i + 1)
            record['id'] = i
            bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(2000, 2000))
            record.setBBox(bbox)
            self.mycatalog.append(record)

        mypsf = measAlg.CoaddPsf(self.mycatalog, self.wcsref, 'weight')

        self.assertEqual(mypsf.computeKernelImage().getBBox(), mypsf.computeBBox())

        with self.assertRaises(pexExceptions.LogicError):
            mypsf.resized(100, 100)

    def testLargeTransform(self):
        """Test that images with bad astrometry are identified"""
        multiplier = 1000.0  # CD matrix multiplier for bad input
        badId = 1  # ID of bad input
        for ii in range(3):
            record = self.mycatalog.addNew()
            record.setPsf(measAlg.DoubleGaussianPsf(50, 50, 5.0, 1.00, 0.0))
            cdMatrix = self.cdMatrix
            if ii == badId:
                # This image has bad astrometry:
                cdMatrix *= multiplier
            record['id'] = ii
            record['weight'] = 1.0
            record.setWcs(afwGeom.makeSkyWcs(crpix=self.crpix, crval=self.crval, cdMatrix=cdMatrix))
            record.setBBox(lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(2000, 2000)))

        coaddPsf = measAlg.CoaddPsf(self.mycatalog, self.wcsref)
        with self.assertRaises(pexExceptions.RangeError) as cm:
            coaddPsf.computeKernelImage()
        self.assertIn("id=%d" % (badId,), str(cm.exception))


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
