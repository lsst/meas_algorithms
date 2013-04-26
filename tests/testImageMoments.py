#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2013 LSST Corporation.
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
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import unittest
import numpy
import os

import lsst.utils.tests
import lsst.afw.geom
import lsst.afw.geom.ellipses as el
import lsst.afw.image
from lsst.meas.algorithms import ImageMoments

class ImageMomentsTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        self.ellipseCores = [
            el.Quadrupole(27.0, 22.0, -5.0),
            el.Quadrupole(23.0, 28.0, 2.0),
            ]
        self.centers = [
            lsst.afw.geom.Point2D(2.0, 3.0),
            lsst.afw.geom.Point2D(-1.0, 2.5),
            ]
        for ellipseCore in self.ellipseCores:
            ellipseCore.scale(2)
        self.bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-500, -500), lsst.afw.geom.Point2I(50, 50))
        self.xg, self.yg = numpy.meshgrid(
            numpy.arange(self.bbox.getBeginX(), self.bbox.getEndX(), dtype=float),
            numpy.arange(self.bbox.getBeginY(), self.bbox.getEndY(), dtype=float)
            )
        self.mom = ImageMoments(5.0)
    
    @staticmethod
    def computeNumericJacobian(func, x, epsilon=1E-6):
        y = func(x)
        result = numpy.zeros(y.shape + x.shape, dtype=float)
        for ix in range(x.size):
            h = epsilon * (x[ix] + epsilon)
            x[ix] += h
            dy = func(x)
            x[ix] -= 2.0*h
            dy -= func(x)
            x[ix] += h
            dy /= 2.0*h
            result[:,ix] = dy
        return result

    def evaluateGaussian(self, ellipse):
        gt = ellipse.getGridTransform()
        xt = gt[gt.XX] * self.xg + gt[gt.XY] * self.yg + gt[gt.X]
        yt = gt[gt.YX] * self.xg + gt[gt.YY] * self.yg + gt[gt.Y]
        return numpy.exp(-0.5 * (xt**2 + yt**2))

    def checkMoments(self, dEllipseCore, dCenter, wEllipseCore, wCenter):
        dEllipse = el.Ellipse(dEllipseCore, dCenter)
        image = lsst.afw.image.ImageD(self.bbox)
        image.getArray()[:,:] = self.evaluateGaussian(dEllipse)
        result = self.mom.measureEllipse(image, wEllipseCore, wCenter)
        return result

    def testGeneral(self):
        """Test every step of measuring moments, but without looking at the uncertainty estimates.
        """
        for dEllipseCore in self.ellipseCores:
            for dCenter in self.centers:
                dEllipse = el.Ellipse(dEllipseCore, dCenter)
                dArray = self.evaluateGaussian(dEllipse)
                dImage = lsst.afw.image.ImageD(self.bbox)
                dImage.getArray()[:,:] = dArray
                for wEllipseCore in self.ellipseCores:
                    for wCenter in self.centers:
                        wEllipse = el.Ellipse(wEllipseCore, wCenter)
                        wArray = self.evaluateGaussian(wEllipse)
                        product = dArray * wArray
                        wSum = wArray.sum()
                        product /= wSum

                        # test raw moments calculation against naive NumPy implementation
                        xr = self.xg - wCenter.getX()
                        yr = self.yg - wCenter.getY()
                        q0 = numpy.sum(product)
                        qx = numpy.sum(product * xr)
                        qy = numpy.sum(product * yr)
                        qxx = numpy.sum(product * xr**2)
                        qyy = numpy.sum(product * yr**2)
                        qxy = numpy.sum(product * xr * yr)
                        q = numpy.array([q0, qxx, qyy, qxy, qx, qy])
                        qMeas = self.mom.measureRaw(dImage, wEllipseCore, wCenter)
                        self.assertClose(qMeas.moments, q, rtol=1E-5, atol=1E-9)

                        # test conversion from raw moments quadrupole/centroid moments
                        m, dm_dq = self.mom.convertRawMoments(qMeas.moments)
                        self.assertEqual(dm_dq.shape, (5,6))
                        mx = qx / q0
                        my = qy / q0
                        mxx = numpy.sum(product * (xr - mx)**2) / q0
                        myy = numpy.sum(product * (yr - my)**2) / q0
                        mxy = numpy.sum(product * (xr - mx) * (yr - my)) / q0
                        m1 = numpy.array([mxx, myy, mxy, mx, my])
                        self.assertClose(m1[3:], m[3:], rtol=1E-5, atol=1E-9)
                        self.assertClose(m1[:3], m[:3], rtol=1E-6, atol=1E-8)
                        # test that doing a single pass over the data isn't hurting us
                        m2 = numpy.array([qxx/q0 - mx**2, qyy/q0 - my**2, qxy/q0 - mx*my, mx, my])
                        self.assertClose(m1, m2, rtol=1E-8, atol=1E-13)

                        c, dc_dm = self.mom.correctWeightedMoments(wEllipseCore, m)
                        self.assertEqual(dc_dm.shape, (5,5))
                        self.assertClose(c[:3], dEllipseCore.getParameterVector(), rtol=5E-7, atol=1E-6)
                        self.assertClose(c[3:], numpy.array(dCenter - wCenter), rtol=1E-8, atol=2E-7)

                        cMeas = self.mom.measureEllipse(dImage, wEllipseCore, wCenter)
                        self.assertClose(cMeas.quadrupole.getParameterVector(),
                                         dEllipseCore.getParameterVector(),
                                         rtol=5E-7, atol=1E-6)
                        self.assertClose(numpy.array(cMeas.centroid),
                                         numpy.array(dCenter - wCenter),
                                         rtol=1E-8, atol=2E-7)
                        

    def testOptimalWeight(self):
        """Test with stricter precision requirements in the case where the weighted-moment correction
        is trivial because the weight and the data are the same.
        """
        for ellipseCore in self.ellipseCores:
            for center in self.centers:
                result = self.checkMoments(ellipseCore, center, ellipseCore, center)
                self.assertClose(result.quadrupole.getParameterVector(),
                                 ellipseCore.getParameterVector(),
                                 rtol=1E-7, atol=1E-9)
                self.assertClose(numpy.array(result.centroid), 0.0, atol=2E-8)

    def testCovarianceWrappers(self):
        """Test Swig wrappers for the covariances in the result objects.
        """
        for dEllipseCore in self.ellipseCores:
            for dCenter in self.centers:
                dEllipse = el.Ellipse(dEllipseCore, dCenter)
                dArray = self.evaluateGaussian(dEllipse)
                dImage = lsst.afw.image.ImageD(self.bbox)
                dMaskedImage = lsst.afw.image.MaskedImageD(self.bbox)
                dImage.getArray()[:,:] = dArray
                dMaskedImage.getImage().getArray()[:,:] = dArray
                for wEllipseCore in self.ellipseCores:
                    for wCenter in self.centers:
                        raw1 = self.mom.measureRaw(dImage, wEllipseCore, wCenter)
                        raw2 = self.mom.measureRaw(dMaskedImage, wEllipseCore, wCenter)
                        self.assertClose(raw1.moments, raw2.moments)
                        self.assertIsNone(raw1.covariance)
                        self.assertEqual(raw2.covariance.shape, (6,6))
                        self.assertClose(raw2.covariance, 0.0)
                        ell1 = self.mom.measureEllipse(dImage, wEllipseCore, wCenter)
                        ell2 = self.mom.measureEllipse(dMaskedImage, wEllipseCore, wCenter)
                        self.assertClose(ell1.moments, ell2.moments)
                        self.assertIsNone(ell1.covariance)
                        self.assertEqual(ell2.covariance.shape, (5,5))
                        self.assertClose(ell2.covariance, 0.0)



    def testDerivatives(self):
        """Test derivatives of convertRawMoments and correctWeightedMoments.

        This is essentially a test of the uncertainty calibration - because the uncertainties
        are non-Gaussian, this is as close as we can really get to testing the correctness of
        the uncertainties without a whole lot of Monte Carlo.
        Happily, the bits we can't test are pretty small and bulletproof, especially compared
        to the horribly complicated derivatives we test here.
        """
        for dEllipseCore in self.ellipseCores:
            for dCenter in self.centers:
                dEllipse = el.Ellipse(dEllipseCore, dCenter)
                dArray = self.evaluateGaussian(dEllipse)
                dImage = lsst.afw.image.ImageD(self.bbox)
                dImage.getArray()[:,:] = dArray
                for wEllipseCore in self.ellipseCores:
                    for wCenter in self.centers:
                        qMeas = self.mom.measureRaw(dImage, wEllipseCore, wCenter)
                        # Test derivative of convertRawMoments
                        def func1(q):
                            m, _ = ImageMoments.convertRawMoments(q)
                            return m
                        jNumeric = self.computeNumericJacobian(func1, qMeas.moments)
                        mMeas, jAnalytic = ImageMoments.convertRawMoments(qMeas.moments)
                        self.assertClose(jAnalytic, jNumeric, rtol=1E-4, atol=1E-7)

                        # Test derivative of correctWeightedMoments
                        def func2(m):
                            c, _ = ImageMoments.correctWeightedMoments(wEllipseCore, m)
                            return c
                        jNumeric = self.computeNumericJacobian(func2, mMeas)
                        cMeas, jAnalytic = ImageMoments.correctWeightedMoments(wEllipseCore, mMeas)
                        self.assertClose(jAnalytic, jNumeric, rtol=6E-3, atol=1E-7)

def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(ImageMomentsTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
