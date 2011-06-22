#!/usr/bin/env python
import math

import lsst.pex.policy as pexPolicy
import lsst.afw.detection as afwDet
import lsst.afw.image as afwImage
import lsst.meas.algorithms as measAlg

model = "SingleGaussian"
size = 15
fwhm = 7.1428571432942825

mi = afwImage.MaskedImageF("mi.fits")

psf = afwDet.createPsf(model, size, size, fwhm/(2*math.sqrt(2*math.log(2))))

policy = pexPolicy.Policy()
policy.set("cond3_fac2", 0.6)
policy.set("minSigma", 6.0)
policy.set("cond3_fac", 2.5)
policy.set("min_DN", 150.0)
policy.set("nCrPixelMax", 200000)
policy.set("niteration", 3)

bg = -0.0030176613945513964
keepCRs = True

measAlg.findCosmicRays(mi, psf, bg, policy, keepCRs)
