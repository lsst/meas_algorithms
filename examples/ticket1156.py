#!/usr/bin/env python

import glob, math, os, sys, math
import eups

from lsst.pex.logging import Log
import lsst.pex.policy as policy
import lsst.afw.detection as afwDetection
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.meas.algorithms as algorithms
import lsst.meas.algorithms.Psf as Psf
import lsst.meas.algorithms.defects as defects

def main():

    logger = Log.getDefaultLog()
    logger.setThreshold(Log.FATAL)
    logger.setThresholdFor("meas.algorithms.measureSources", Log.FATAL)

    # load the fits file, make gausssian psf
    fileName = os.path.join(eups.productDir("afwdata"), "CFHT", "D4", "cal-53535-i-797722_%d" % 1)
    mi = afwImage.MaskedImageF(fileName, 0, None, afwImage.BBox(afwImage.PointI(824, 140), 512, 512))
    exposure = afwImage.makeExposure(mi)
    psf        = algorithms.createPSF("DoubleGaussian", 15, 15, 2.0)

    # Fix defects
    pafFile = os.path.join(eups.productDir("meas_algorithms"), "policy", "BadPixels.paf")
    algorithms.interpolateOverDefects(mi, psf, defects.policyToBadRegionList(pafFile))

    # Subtract background
    backobj = afwMath.makeBackground(mi.getImage(),afwMath.BackgroundControl("NATURAL_SPLINE", 3, 3))
    img = mi.getImage(); img -= backobj.getImageF(); del img

    # Smooth 
    cnvImage = mi.Factory(mi.getDimensions())
    cnvImage.setXY0(afwImage.PointI(mi.getX0(), mi.getY0()))
    psf.convolve(cnvImage, mi, True, cnvImage.getMask().getMaskPlane("EDGE"))

    # detect
    threshold = afwDetection.Threshold(3, afwDetection.Threshold.STDEV)
    fs = afwDetection.FootprintSetF(cnvImage, threshold, "DETECTED")
    footprints = fs.getFootprints()

    # measure
    pafFile = os.path.join(eups.productDir("meas_pipeline"), "policy", "MeasureSources.paf")
    moPolicy = policy.Policy.createPolicy(pafFile).getPolicy("measureObjects")
    measureSources = algorithms.makeMeasureSources(exposure, moPolicy, psf)

    for i in range(len(footprints)):
        source = afwDetection.Source()
        measureSources.apply(source, footprints[i])
            
if __name__ == "__main__":
    main()
