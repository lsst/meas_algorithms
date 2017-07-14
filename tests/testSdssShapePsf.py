#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2017 AURA/LSST.
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

from __future__ import absolute_import, division, print_function
import unittest

import numpy as np

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg
import lsst.meas.base as measBase
import lsst.meas.base.tests as measBaseTests
import lsst.utils.tests


class SdssShapePsfTestCase(measBaseTests.AlgorithmTestCase, lsst.utils.tests.TestCase):
    """Test case to ensure base_SdssShape_psf is being measured at source position

    Note: this test lives here in meas_algorithms rather than meas_base (where SdssShape
    lives) due to the need to apply a spatially varying PSF model such that the PSF is
    different at each position in the image.  This varying PSF model is built with
    meas_algorithms' PcaPsf (which is not accessible from meas_base).
    """
    def setUp(self):
        self.bbox = afwGeom.Box2I(afwGeom.Point2I(-20, -30), afwGeom.Extent2I(240, 160))
        self.dataset = measBaseTests.TestDataset(self.bbox)
        # first two sources are points
        self.pointCentroid1 = afwGeom.Point2D(50.1, 49.8)
        self.pointCentroid2 = afwGeom.Point2D(-11.6, -1.7)
        self.dataset.addSource(flux=1E5, centroid=self.pointCentroid1)
        self.dataset.addSource(flux=2E5, centroid=self.pointCentroid2)
        # third source is extended
        self.extendedCentroid = afwGeom.Point2D(149.9, 50.3)
        self.dataset.addSource(flux=1E5, centroid=self.extendedCentroid,
                               shape=afwGeom.ellipses.Quadrupole(8, 9, 3))
        self.config = self.makeSingleFrameMeasurementConfig("base_SdssShape")

    def tearDown(self):
        del self.bbox
        del self.dataset
        del self.pointCentroid1
        del self.pointCentroid2
        del self.extendedCentroid
        del self.config

    def _computeVaryingPsf(self):
        """Compute a varying PSF as a linear combination of PCA (== Karhunen-Loeve) basis functions

        We simply desire a PSF that is not constant across the image, so the precise choice of
        parameters (e.g., sigmas, setSpatialParameters) are not crucial.
        """
        kernelSize = 31
        sigma1 = 1.75
        sigma2 = 2.0*sigma1
        basisKernelList = []
        for sigma in (sigma1, sigma2):
            basisKernel = afwMath.AnalyticKernel(kernelSize, kernelSize,
                                                 afwMath.GaussianFunction2D(sigma, sigma))
            basisImage = afwImage.ImageD(basisKernel.getDimensions())
            basisKernel.computeImage(basisImage, True)
            basisImage /= np.sum(basisImage.getArray())
            if sigma == sigma1:
                basisImage0 = basisImage
            else:
                basisImage -= basisImage0
            basisKernelList.append(afwMath.FixedKernel(basisImage))

        order = 1
        spFunc = afwMath.PolynomialFunction2D(order)
        exactKernel = afwMath.LinearCombinationKernel(basisKernelList, spFunc)
        exactKernel.setSpatialParameters([[1.0, 0, 0], [0.0, 0.5E-2, 0.2E-2]])
        exactPsf = measAlg.PcaPsf(exactKernel)

        return exactPsf

    def _runMeasurementTask(self, psf=None):
        task = self.makeSingleFrameMeasurementTask("base_SdssShape", config=self.config)
        exposure, catalog = self.dataset.realize(10.0, task.schema)
        if psf:
            exposure.setPsf(psf)
        task.run(catalog, exposure)
        return exposure, catalog

    def _checkPsfShape(self, result, psfResult, psfTruth):
        self.assertFloatsAlmostEqual(psfResult.getIxx(), psfTruth.getIxx(), rtol=1E-4)
        self.assertFloatsAlmostEqual(psfResult.getIyy(), psfTruth.getIyy(), rtol=1E-4)
        self.assertFloatsAlmostEqual(psfResult.getIxy(), psfTruth.getIxy(), rtol=1E-4)
        self.assertFalse(result.getFlag(measBase.SdssShapeAlgorithm.PSF_SHAPE_BAD.number))

    def testMeasureGoodPsf(self):
        """Test that we measure shapes and record the PSF shape correctly

        To ensure this, apply a varying PSF to the image such that different positions
        can be distinguished by their different PSF model shapes.
        """
        # Apply varying PSF model to the exposure
        varyingPsf = self._computeVaryingPsf()
        exposure, catalog = self._runMeasurementTask(psf=varyingPsf)
        key = measBase.SdssShapeResultKey(catalog.schema["base_SdssShape"])
        # First make sure we did indeed get a varying PSF model across the exposure
        psf = exposure.getPsf()
        # Compare truth PSF at positions of two point sources
        self.assertFloatsNotEqual(psf.computeShape(self.pointCentroid1).getIxx(),
                                  psf.computeShape(self.pointCentroid2).getIxx(), rtol=1E-1)
        self.assertFloatsNotEqual(psf.computeShape(self.pointCentroid1).getIyy(),
                                  psf.computeShape(self.pointCentroid2).getIyy(), rtol=1E-1)
        # Compare truth PSF at average position vs. truth PSF at extended source position
        self.assertFloatsNotEqual(psf.computeShape(self.extendedCentroid).getIxx(),
                                  psf.computeShape().getIxx(), rtol=1E-1)
        self.assertFloatsNotEqual(psf.computeShape(self.extendedCentroid).getIyy(),
                                  psf.computeShape().getIyy(), rtol=1E-1)
        # Now check the base_SdssShape_psf entries against the PSF truth values
        for record in catalog:
            psfTruth = psf.computeShape(afwGeom.Point2D(record.getX(), record.getY()))
            result = record.get(key)
            psfResult = key.getPsfShape(record)
            self._checkPsfShape(result, psfResult, psfTruth)

    def testResizedPcaPsf(self):
        """Test that  PcaPsf can resize itself.

        This test resides here because PcaPsfs do not have their own test module"""
        psf = self._computeVaryingPsf()
        dim = psf.computeBBox().getDimensions()
        for pad in [0, 4, -2]:
            resizedPsf = psf.resized(dim.getX() + pad, dim.getY() + pad)
            self.assertEqual(resizedPsf.computeBBox().getDimensions(),
                             afwGeom.Extent2I(dim.getX() + pad, dim.getY() + pad))
            if psf.getKernel().isSpatiallyVarying():
                self.assertEqual(resizedPsf.getKernel().getSpatialParameters(),
                                 psf.getKernel().getSpatialParameters())
            else:
                self.assertEqual(resizedPsf.getKernel().getKernelParameters(),
                                 psf.getKernel().getKernelParameters())
            self._compareKernelImages(resizedPsf, psf)

    def _compareKernelImages(self, psf1, psf2):
        """Test that overlapping portions of kernel images are identical
        """
        # warning: computeKernelImage modifies kernel parameters if spatially varying
        im1 = psf1.computeKernelImage()
        im2 = psf2.computeKernelImage()
        bboxIntersection = im1.getBBox()
        bboxIntersection.clip(im2.getBBox())
        im1Intersection = afwImage.ImageD(im1, bboxIntersection)
        im2Intersection = afwImage.ImageD(im2, bboxIntersection)
        scale1 = im1.getArray().sum() / im1Intersection.getArray().sum()
        scale2 = im2.getArray().sum() / im2Intersection.getArray().sum()
        im1Arr = scale1 * im1Intersection.getArray()
        im2Arr = scale2 * im2Intersection.getArray()
        self.assertTrue(np.allclose(im1Arr, im2Arr),
                        "kernel images %s, %s do not match" % (im1Arr, im2Arr))


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
